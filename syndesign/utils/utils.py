
"""
SynDesign - Configuration and Utility Section
Author: PJM
Description: Basic constants, environment setup, and shared utility functions
"""

# === Imports ===
import os, sys, re, gzip
import pickle
import torch
import regex
import pandas as pd
import numpy as np
from collections import defaultdict
from scipy import stats
from sklearn.preprocessing import minmax_scale
from itertools import cycle
import multiprocessing as mp
import matplotlib as mpl
import matplotlib.pyplot as plt

# === Set non-interactive matplotlib backend ===
mpl.use('Agg')

# === Environment Configuration ===
BASE_DIR = os.getcwd()
os.environ['MPLCONFIGDIR'] = f'{BASE_DIR}/images/matplt_tmp'

# === Hardware Detection ===
NUM_CPUs = os.cpu_count()
NUM_GPUs = torch.cuda.device_count()

# === Constants ===
ALT_INDEX = 60  # Number of nucleotides flanking the variant
CORES = 1

# === Utility Functions ===

dict_sAA_TABLE = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                  'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                  'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                  'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                  'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                  'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                  'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                  'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                  'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                  'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
                  'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
                  }


# def END: reverse_complement

def translate(seq):
    # Ensure the sequence length is divisible by 3 (codon length)
    if len(seq) % 3 != 0:
        raise ValueError("Input DNA sequence length must be a multiple of 3.")

    AAseq = ''
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        AAseq += dict_sAA_TABLE.get(codon, 'X')  # 'X' for unknown codons

    return AAseq


# def END: translate

def reverse_complement(seq):
    dict_bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'U': 'U', 'n': '',
                  '.': '.', '*': '*', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', '>': '>', 'X': 'X'}
    list_out = [dict_bases[base] for base in list(seq)]
    return ''.join(list_out)[::-1]


def lower_list(input: list): return [v.lower() for v in input]


def lower_dict(input: dict):
    return dict((k.lower(), v.lower) for k, v in input.items())


re_nonchr = re.compile('[^a-zA-Z]')


class Fasta:
    def __init__(self, fasta):

        # V-S Check
        if not os.path.isfile(fasta):
            sys.exit('(): File does not exist')

        self.fafile = open(fasta, 'r')
        self.chrlist = []
        self.chrlen = []
        self.seekpos = []
        self.len1 = []
        self.len2 = []

        # V-S Check
        if not os.path.isfile('%s.fai' % fasta):
            sys.exit('.fai file does not exist')

        faifile = open('%s.fai' % fasta, 'r')
        for line in faifile:
            columns = line.strip('\n').split()  # Goes backwards, -1 skips the new line character

            self.chrlist.append(columns[0])
            self.chrlen.append(int(columns[1]))
            self.seekpos.append(int(columns[2]))
            self.len1.append(int(columns[3]))
            self.len2.append(int(columns[4]))
        # loop END: sLINE
        faifile.close()
        self.type = []

    # def END: __init_

    def fetch(self, chrom, start=None, end=None, strand='+'):

        assert chrom in self.chrlist, chrom

        index = self.chrlist.index(chrom)

        if start == None: nFrom = 0
        if end == None:   nTo = self.chrlen[index]
        # if nTo >= self.nChromLen[nChrom]: nTo = self.nChromLen[nChrom]-1

        assert (0 <= start) and (start < end) and (end <= self.chrlen[index])

        nBlank = self.len2[index] - self.len1[index]

        start = int(start + (start / self.len1[index]) * nBlank)  # Start Fetch Position

        end = int(end + (end / self.len1[index]) * nBlank)  # End Fetch Position

        self.fafile.seek(self.seekpos[index] + start)  # Get Sequence

        seq = re.sub(re_nonchr, '', self.fafile.read(end - start))

        if strand == '+':
            return seq
        elif strand == '-':
            return reverse_complement(seq)
        else:
            sys.exit('Invalid Strand')
    # def END: fetch


# class END: Fasta

class RefAA:
    def __init__(self, genesym):
        ref = '%s/ref' % os.getcwd()
        genome = 'GRCh38'
        fa = '%s/%s.fa' % (ref, genome)
        fasta = Fasta(fa)
        idfile = '%s/forwebtool/%s/%s' % (ref, 'GeneSym', genesym)
        inf = open(idfile, 'r')
        self.gi, self.genesym, self.chrom, self.strand, self.exons = inf.readline().strip('\n').split('\t')
        inf.close()

        self.exoncnt = len(self.exons.split(','))
        self.coord = dict([pos.split('|') for pos in self.exons.split(',')])
        self.dict_refAA = self.fetch_orf(fasta)

        # def END: __init_

    def fetch_orf(self, fasta):

        fullseq = ''
        dict_index = {}
        exonsum = 0
        frame_numbers = cycle([0, 1, 2])
        for exon_no, poskey in self.coord.items():
            s, e = [int(pos) for pos in poskey.split('-')]
            seq = fasta.fetch(self.chrom, s - 1, e, self.strand).upper()

            for i in range(len(seq)):
                adjusted_i = i + exonsum
                current_frame = next(frame_numbers)

                if adjusted_i not in dict_index:
                    dict_index[adjusted_i] = ''

                if self.strand == '+':
                    genomic = s + i
                    refnuc = seq[i]

                else:
                    genomic = e - i
                    refnuc = reverse_complement(seq[i])
                # Extract codon and codon index information
                dict_index[adjusted_i] = {
                    'exon_no': exon_no,
                    'exonseq': seq,
                    'refnuc': refnuc,
                    'index': i,
                    'gpos': genomic,
                    'frame': current_frame,
                    'strand': self.strand,
                }
            # loop end: i
            exonsum += e - s + 1
            fullseq += seq
        # loop END:

        for i in range(0, len(fullseq), 3):
            codon = fullseq[i:i + 3]
            codon_i = i // 3
            aa = dict_sAA_TABLE[codon]
            for i in range(i, i + 3):
                dict_index[i]['codon'] = codon
                dict_index[i]['anticodon'] = reverse_complement(codon)
                dict_index[i]['codon_i'] = codon_i + 1  # start codon = 1
                dict_index[i]['aa'] = aa
            # loop END: i
        # loop END: i

        dict_refAA = {}
        for i, info in dict_index.items():

            gpos = info['gpos']
            if gpos not in dict_refAA:
                dict_refAA[gpos] = {'exon_no': info['exon_no'],
                                    'exonseq': info['exonseq'],
                                    'refnuc': info['refnuc'],
                                    'index': info['index'],
                                    'frame': info['frame'],
                                    'codon': info['codon'],
                                    'anticodon': info['anticodon'],
                                    'codon_i': info['codon_i'],
                                    'aa': info['aa'],
                                    'strand': info['strand']}

            # loop END: i, info

        return dict_refAA
    # def END: fetch_orf


# class END: RefGenome


class cVCFData:
    def __init__(self):
        self.sPatID = ''
        self.sChrID = ''
        self.nPos = 0
        self.sDBSNP_ID = ''
        self.sRefNuc = ''
        self.sAltNuc = ''
        self.fQual = 0.0
        self.sFilter = ''
        self.sInfo = ''
        self.sFormat = ''
        self.list_sMisc = []

        # Extra
        self.nClusterID = 0

    # def END: __int__

    def parse_vcf_file_v2(self, infile):
        if not os.path.isfile(infile):
            sys.exit('File Not Found %s' % infile)

        dict_out = {}

        if infile.endswith('.gz'):
            inf = gzip.open(infile, 'rt')
        else:
            inf = open(infile, 'r')

        for line in inf:

            # File Format
            # Column Number:     | 0       | 1        | 2          | 3       | 4
            # Column Description:| sChrID  | nPos     | sDBSNP_ID  | sRefNuc | sAltNuc
            # Column Example:    | 1       | 32906558 | rs79483201 | T       | A
            # Column Number:     | 5       | 6        | 7          | 8              | 9./..
            # Column Description:| fQual   | sFilter  | sInfo      | sFormat        | sSampleIDs
            # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to sFormat

            if line.startswith('#'): continue  # SKIP Information Headers
            list_col = line.strip('\n').split('\t')

            self.chrom = 'chr%s' % list_col[0]

            try:
                self.pos = int(list_col[1])
            except ValueError:
                continue

            self.varID = int(list_col[2])
            self.ref = list_col[3]
            self.alt = list_col[4]
            self.qual = float(list_col[5]) if list_col[5] != '.' else list_col[5]
            self.filter = list_col[6]
            self.addinfo = list_col[7]

            dict_addinfo = dict([info.split('=') for info in self.addinfo.split(';') if len(info.split('=')) == 2])

            self.vartype = dict_addinfo['CLNVC'].lower()
            try:
                self.geneinfo = dict_addinfo['GENEINFO'].upper()
            except KeyError:
                self.geneinfo = 'N/A'

            if self.varID not in dict_out:
                dict_out[self.varID] = ''
            dict_out[self.varID] = self

        # loop END: sReadLine
        inf.close()

        return dict_out

    # def END: parse_vcf_file

    def parse_vcf_stdout2(self, stdout):
        list_sOutput = []
        for sReadLine in stdout:
            # File Format
            # Column Number:     | 0       | 1        | 2          | 3       | 4
            # Column Description:| sChrID  | nPos     | sDBSNP_ID  | sRefNuc | sAltNuc
            # Column Example:    | chr13   | 32906558 | rs79483201 | T       | A
            # Column Number:     | 5       | 6        | 7          | 8              | 9./..
            # Column Description:| fQual   | sFilter  | sInfo      | sFormat        | sSampleIDs
            # Column Example:    | 5645.6  | PASS     | .          | GT:AD:DP:GQ:PL | Scores corresponding to sFormat
            ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Ph
            sReadLine = str(sReadLine, 'UTF-8')
            list_sColumn = sReadLine.strip('\n').split('\t')

            if list_sColumn[0] == 'MT': continue
            if list_sColumn[0].startswith('<GL00'): continue

            dict_sChrKey = {'X': '23', 'Y': '24'}
            self.sChrID = 'chr%s' % list_sColumn[0]

            if list_sColumn[0] in ['X', 'Y']:
                self.nChrID = int(dict_sChrKey[list_sColumn[0]])
            else:
                self.nChrID = int(list_sColumn[0])
            self.nPos = int(list_sColumn[1])
            self.sDBSNP_ID = list_sColumn[2]
            self.sRefNuc = list_sColumn[3]
            self.sAltNuc = list_sColumn[4]
            self.fQual = float(list_sColumn[5]) if list_sColumn[5] != '.' else list_sColumn[5]
            self.sFilter = list_sColumn[6]
            self.sInfo = list_sColumn[7]
            dict_sInfo = dict([sInfo.split('=') for sInfo in self.sInfo.split(';') if len(sInfo.split('=')) == 2])

            try:
                self.fAlleleFreq = float(dict_sInfo['AF_raw'])
            except ValueError:
                self.fAlleleFreq = np.mean([float(f) for f in dict_sInfo['AF_raw'].split(',')])

            list_sOutput.append(self)
        # loop END: sReadLine

        return list_sOutput
    # def END: parse_vcf_stdout2


# class END: cVCFData


class SplitFastq:
    def __init__(
            self,
            file: str,
            n_split: int,
            out_name: str,
            out_path: str = './',
            silence: bool = False,
    ):
        """fastq file을 원하는 수 만큼 균등하게 나눠주는 함수.

        Args:
            file (str): fastq 파일 경로
            n_split (int): 몇 등분 할 것인지 적는 칸
            out_name (str): 나눴을 때 저장되는 파일들의 prefix
            out_path (str, optional): Output이 저장 될 경로. Defaults to './'.
            silence (bool, optional): Logging을 위한 print 되는 메시지를 끄는 용도. Defaults to False.
        """

        output_format = 'fastq'
        lineset = 4

        self.names = []
        self.dir = '%s/%s_subsplits' % (os.path.abspath(out_path), out_name)
        os.makedirs(self.dir, exist_ok=True)

        with open(file, 'r') as f:
            lines = f.readlines()
            total = len(lines)
            rec_cnt = total / lineset

            list_nBins = [[int(rec_cnt * (i + 0) / n_split), int(rec_cnt * (i + 1) / n_split)] for i in range(n_split)]
            self.meta = {}
            cnt = 0

            for nStart, nEnd in list_nBins:
                if silence == False: print('[Info] Make data subsplits: %s - %s' % (nStart, nEnd))

                sSplit_file_name = '%s_%s.%s' % (out_name, cnt, output_format)
                with open('%s/%s' % (self.dir, sSplit_file_name), 'w') as outfile:
                    for l in lines[nStart * lineset:nEnd * lineset]: outfile.write(l)

                self.names.append(sSplit_file_name)

                self.meta[sSplit_file_name] = {
                    'start': nStart,
                    'end': nEnd,
                    'count': nEnd - nStart
                }
                cnt += 1


# class END: SplitFastq


def load_HGNC_reference(infile, output_by_genesym=0):
    dict_out = {}
    inf = open(infile, 'r')

    for line in inf:
        ## File Format ##
        # HGNC:ID	Approved symbol
        # HGNC:5	A1BG

        if line.startswith('HGNC:ID'): continue

        list_col = line.strip('\n').split('\t')

        HGNC_ID = list_col[0].split(':')[1]
        genesym = list_col[1].upper()

        if output_by_genesym:
            if genesym not in dict_out:
                dict_out[genesym] = ''
            dict_out[genesym] = HGNC_ID
        else:
            if HGNC_ID not in dict_out:
                dict_out[HGNC_ID] = ''
            dict_out[HGNC_ID] = genesym
        # if END:
    # loop END: line
    return dict_out


# def END: load_HGNC_reference
def safe_open(path, mode='r', encoding='utf-8'):
    """Open file or gzip transparently based on extension."""
    if path.endswith('.gz'):
        return gzip.open(path, mode + 't', encoding=encoding)
    return open(path, mode, encoding=encoding)

def load_indexed_pickle(index_file: str, variant_id: int):
    """
    Load specific data entry from an indexed pickle file based on variant ID.
    Parameters:
        index_file (str): Path to the index.txt file.
        variant_id (int): COSMIC or ClinVar ID to lookup.
    Returns:
        Any: Object mapped to the variant ID.
    """
    with open(index_file, 'r') as inf:
        for line in inf:
            idx_range, pickle_path = line.strip().split('\t')
            start, end = map(int, idx_range.split('-'))
            if start <= variant_id <= end:
                with open(pickle_path, 'rb') as pf:
                    return pickle.load(pf)[variant_id]
    raise KeyError(f"Variant ID {variant_id} not found in index.")

def extract_mutation(alt_type_str):
    """
    Extract a SNV mutation from an HGVS-style string like 'c.35G>T'.
    Parameters:
        alt_type_str (str): Mutation string in 'Ref>Alt' format.
    Returns:
        str or None: 'G>T' or None if not found.
    """
    match = re.search(r'([ACGT]+)>([ACGT]+)', alt_type_str)
    return f"{match.group(1)}>{match.group(2)}" if match else None

def revcom(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    Parameters:
        seq (str): Input DNA sequence.
    Returns:
        str: Reverse complement.
    """
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(complement)[::-1]
