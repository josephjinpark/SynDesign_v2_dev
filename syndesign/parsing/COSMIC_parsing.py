from typing import Dict, Union
import os
import sys
import gzip
import pickle

# === COSMIC Record Class ===

class COSMICRecord:
    """Container for COSMIC mutation entry."""
    def __init__(self):
        self.sGeneName = ''
        self.sAccID = ''
        self.nCDSLen = 0
        self.sHGCNID = ''
        self.sSample = ''
        self.sSampleID = ''
        self.sTumorID = ''
        self.sPriSite = ''
        self.sSiteSub1 = ''
        self.sSiteSub2 = ''
        self.sSiteSub3 = ''
        self.sPriHist = ''
        self.sHistSub1 = ''
        self.sHistSub2 = ''
        self.sHistSub3 = ''
        self.bGenomeWide = ''
        self.sMutaID = ''
        self.sAltType = ''
        self.sRef = ''
        self.sAlt = ''
        self.sAAType = ''
        self.sMutaDescri = ''
        self.sMutaZygo = ''
        self.bLOH = ''
        self.sGRCh = ''
        self.sGenicPos = ''
        self.nChrID = ''
        self.sChrID = ''
        self.sPos = ''
        self.sStrand = ''
        self.bSNP = ''
        self.sDelete = ''


def parse_cosmic(file_path: str) -> Dict[int, list]:
    """Parses COSMIC mutations and returns a dictionary of SNV entries keyed by COSMIC ID."""
    from deepprime_utils_refactored import safe_open, extract_mutation

    dict_out = {}

    with safe_open(file_path, 'r') as inf:
        for line in inf:
            if line.startswith('Gene'):
                continue  # Skip header

            cols = line.strip().split('\t')
            record = COSMICRecord()
            try:
                record.sGeneName = cols[0].upper()
                record.sAccID = cols[1]
                record.nCDSLen = int(cols[2])
                record.sPriSite = cols[7]
                record.sPriHist = cols[11]
                record.bGenomeWide = (cols[15] == 'y')
                record.sMutaID = cols[16]
                record.sAltType = cols[17]
                record.sAAType = cols[18]
                record.sMutaDescri = cols[19]
                record.bLOH = (cols[21] == 'y')
                record.sGRCh = cols[22]
                record.sGenicPos = cols[23]
                record.sStrand = cols[24]
                record.bSNP = (cols[25] == 'y')
            except IndexError:
                continue  # Malformed line

            if not record.sGenicPos:
                continue

            try:
                record.nChrID = record.sGenicPos.split(':')[0]
                chrom_map = {'24': 'chrX', '25': 'chrY'}
                record.sChrID = chrom_map.get(record.nChrID, f'chr{record.nChrID}')

                pos_parts = set(record.sGenicPos.split(':')[1].split('-'))
                record.sPos = next(iter(pos_parts))
            except Exception:
                continue

            altnotation = extract_mutation(record.sAltType)
            if not altnotation or len(altnotation) != 3:
                continue  # skip non-SNVs

            try:
                mutaid = int(record.sMutaID.replace('COSM', ''))
            except ValueError:
                continue

            dict_out[mutaid] = [
                record.sGeneName, altnotation, record.sChrID, record.sStrand, record.sGenicPos
            ]

    return dict_out


class VCFRecord:
    """Container for a single VCF entry."""
    def __init__(self):
        self.chrom = ''
        self.pos = 0
        self.varID = 0
        self.ref = ''
        self.alt = ''
        self.qual = ''
        self.filter = ''
        self.addinfo = ''
        self.vartype = ''
        self.geneinfo = 'N/A'


def parse_vcf_file(file_path: str) -> Dict[int, VCFRecord]:
    """Parses ClinVar VCF and returns a dictionary of records keyed by variant ID."""
    from deepprime_utils_refactored import safe_open

    dict_out = {}
    with safe_open(file_path, 'r') as inf:
        for line in inf:
            if line.startswith('#'):
                continue

            cols = line.strip().split('\t')
            try:
                record = VCFRecord()
                record.chrom = f'chr{cols[0]}'
                record.pos = int(cols[1])
                record.varID = int(cols[2].replace('VCV', '').split('.')[0])
                record.ref = cols[3]
                record.alt = cols[4]
                record.qual = float(cols[5]) if cols[5] != '.' else None
                record.filter = cols[6]
                record.addinfo = cols[7]

                info_dict = dict(
                    kv.split('=') for kv in record.addinfo.split(';') if '=' in kv
                )
                record.vartype = info_dict.get('CLNVC', '').lower()
                record.geneinfo = info_dict.get('GENEINFO', 'N/A').upper()

                dict_out[record.varID] = record
            except (IndexError, ValueError):
                continue  # Skip malformed lines

    return dict_out
