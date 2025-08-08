import os
import sys
import pandas as pd

# Enable local package discovery
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '.')))

from syndesign.pipeline.get_pegrnas import make_dp_input_v2
from syndesign.pipeline.dp_batch import run_deepprime_parallel
from syndesign.pipeline.generate_output import create_html_output_table


def main():
    # Parameters would typically come from a config file or command-line arguments
    input_fasta     = '/path/to/genome.fa'
    input_id        = 'GENE1'
    chrID           = 'chr1'
    strand          = '+'
    coord           = {'1': '1000-1100', '2': '1200-1300'}  # Example exon coords
    pe_system       = 'PE2max-HCT116'
    flank           = 60
    target_exon     = 0  # all exons
    
    results_dir = f'./results/{input_id}/{pe_system}'
    tempdir = f'./temp/{input_id}/{pe_system}'
    model_dir = './models/DeepPrime'
    base_name = input_id

    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(tempdir, exist_ok=True)

    # Split PE string
    pe, celltype = pe_system.split('-', 1)

    # Options dict passed throughout
    opts = {
        'pe': pe,
        'celltype': celltype,
        'edittype': 'sub',
        'editsize': 1,
        'rtt_max': 40,
        'pbs_max': 17,
        'pbs_min': 7,
        'target': target_exon,
        'inputtype': 'GeneSym',
        'NUM_GPUs': 1
    }

    from syndesign.utils import Fasta
    fasta = Fasta(input_fasta)

    # 1. Generate input sequences
    df_input = make_dp_input_v2(fasta, pe_system, chrID, strand, coord, flank, target_exon)

    # 2. Run DeepPrime model scoring
    df_scores = run_deepprime_parallel(tempdir, model_dir, results_dir, df_input, opts)

    # 3. Output plots + HTML-compatible tables
    create_html_output_table(df_scores, pe=pe, top_n=50, results_dir=results_dir, base_name=base_name)


if __name__ == '__main__':
    main()
