from typing import Dict, List, Union
import pandas as pd
from syndesign.utils import revcom

def get_all_sub_combo(fullseq: str, targetseq: str, edit_start: int, flank: int, strand: str) -> List[List[str]]:
    """
    Generate all possible SNV edits for a given window.
    """
    inputseqs = []
    for i in range(len(targetseq)):
        for alt in 'ACGT':
            if targetseq[i] == alt:
                continue

            wtseq = (
                fullseq[edit_start - flank + i:edit_start + i] +
                targetseq[i] +
                fullseq[edit_start + i + 1:edit_start + i + 1 + flank]
            )
            edseq = (
                fullseq[edit_start - flank + i:edit_start + i] +
                alt +
                fullseq[edit_start + i + 1:edit_start + i + 1 + flank]
            )

            if strand == '-':
                wtseq = revcom(wtseq)
                edseq = revcom(edseq)
                refnuc = revcom(targetseq[i])
                altnuc = revcom(alt)
            else:
                refnuc = targetseq[i]
                altnuc = alt

            inputseqs.append([wtseq, edseq, edit_start + i, refnuc, altnuc])

    return inputseqs


def find_PAMs_for_input_v2(
    exon_no: int,
    chrID: str,
    exon_start: int,
    exon_end: int,
    seq: str,
    fasta,
    pe_system: str,
    flank: int,
    dict_out: Dict[int, List[List[str]]]
) -> Dict[int, List[List[str]]]:
    """
    Scan input sequence for PAMs and generate candidate edit windows.
    """
    max_rttlen = 40
    guidelen = 20

    # PAM regexes
    if 'NRCH' in pe_system:
        dict_pam_re = {
            '+': '[ACGT][ACGT]G[ACGT]|[ACGT][CG]A[ACGT]|[ACGT][AG]CC|[ATCG]ATG',
            '-': '[ACGT]C[ACGT][ACGT]|[ACGT]T[CG][ACGT]|G[GT]T[ACGT]|ATT[ACGT]|CAT[ACGT]|GGC[ACGT]|GTA[ACGT]'
        }
    else:
        dict_pam_re = {'+': '[ACGT]GG', '-': 'CC[ACGT]'}

    import regex

    for strand in ['+', '-']:
        pam_re = dict_pam_re[strand]
        for match in regex.finditer(pam_re, seq, overlapped=True):
            i_start = match.start()
            i_end = match.end()

            if strand == '+':
                nickpos = i_start - 3
                nickpos = max(0, nickpos)
                winsize = min(len(seq), nickpos + max_rttlen)
                alt_window = [exon_start + nickpos, exon_start + winsize]
            else:
                nickpos = i_end + 3
                nickpos = min(len(seq), nickpos)
                winsize = max(0, nickpos - max_rttlen)
                alt_window = [exon_start + winsize, exon_start + nickpos]

            # Extract sequence from genome
            targetseq = fasta.fetch(chrID, alt_window[0] - 1, alt_window[1] - 1).upper()

            # Generate all SNV candidates for window
            inputseqs = get_all_sub_combo(fullseq=seq, targetseq=targetseq, edit_start=alt_window[0], flank=flank, strand=strand)

            for wt, ed, editpos, refnuc, altnuc in inputseqs:
                guidekey = f'{exon_no}.{editpos}|{chrID}:{editpos}|{refnuc}>{altnuc}'
                if editpos not in dict_out:
                    dict_out[editpos] = []
                dict_out[editpos].append([guidekey, exon_no, wt, ed])

    return dict_out


def make_dp_input_v2(
    fasta,
    pe_system: str,
    chrID: str,
    strand: str,
    coord: Union[str, Dict[str, str]],
    flank: int,
    target: int
) -> pd.DataFrame:
    """
    Construct input dataframe for a given target exon or full transcript.
    """
    dict_out = {}
    if target == 0:
        # All exons
        list_targets = []
        for exon_no, poskey in coord.items():
            s, e = map(int, poskey.split('-'))
            seq = fasta.fetch(chrID, s - 1, e).upper()
            list_targets.append([int(exon_no), s, e, seq])

        for exon_no, s, e, seq in list_targets:
            find_PAMs_for_input_v2(exon_no, chrID, s, e, seq, fasta, pe_system, flank, dict_out)
    else:
        s, e = map(int, coord.split('-'))
        seq = fasta.fetch(chrID, s - 1, e).upper()
        find_PAMs_for_input_v2(target, chrID, s, e, seq, fasta, pe_system, flank, dict_out)

    list_out = []
    for editpos, guides in dict_out.items():
        list_out += guides[:3]  # limit to 3 guides per edit

    header = ['ID', 'target', 'wtseq', 'edseq']
    df = pd.DataFrame(list_out, columns=header)

    if df.empty:
        raise RuntimeError(f'No PAM sequence was detected in the exon == {target}')

    return df
