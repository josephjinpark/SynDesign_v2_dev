import os
import sys
import pandas as pd
import pickle
import numpy as np
import multiprocessing as mp
from sklearn.preprocessing import minmax_scale
from scipy import stats

def run_deepprime_batch(params):
    """Worker function for multiprocessing."""
    index, target_dir, df, bin_indices, gpu_id, opts = params

    from genet import predict as prd

    pe = opts['pe']
    celltype = opts['celltype']
    edittype = opts['edittype']
    editsize = int(opts['editsize'])
    rtt_max = opts['rtt_max']
    pbs_max = opts['pbs_max']
    pbs_min = opts['pbs_min']
    inputtype = opts['inputtype']
    target = opts['target']

    for i in bin_indices:
        row = df.iloc[i]
        ID, wt, ed = row['ID'], row['wtseq'], row['edseq']
        pam = 'NRCH' if pe.startswith('NRCH') else 'NGG'

        # Output path
        if inputtype in ['GeneSym', 'NMID', 'EnsemblID', 'HGNC'] and target == 0:
            exon = ID.split('.')[0]
            outdir = os.path.join(target_dir, exon)
            os.makedirs(outdir, exist_ok=True)
            output_path = os.path.join(outdir, f'{i}.csv.pkl')
        else:
            output_path = os.path.join(target_dir, f'{i}.csv.pkl')

        dp = prd.DeepPrime(
            sID=ID, Ref_seq=wt, ED_seq=ed,
            edit_type=edittype, edit_len=editsize,
            pam=pam,
            pbs_min=pbs_min, pbs_max=pbs_max,
            rtt_min=0, rtt_max=rtt_max,
            silence=True, gpu=int(gpu_id)
        )

        if dp.pegRNAcnt == 0:
            continue

        df_scores = dp.predict(pe_system=pe, cell_type=celltype, show_features='syn_pe')
        df_scores = df_scores.sort_values(by=f'{pe}_score', ascending=False)
        df_scores['wtseq'] = wt

        with open(output_path, 'wb') as pf:
            pickle.dump(df_scores, pf)


def run_deepprime_parallel(tempdir, model_dir, results_dir, df, opts):
    """
    Runs DeepPrime prediction in parallel and collects result tables.
    """
    from genet import models

    pe = opts['pe']
    target = opts['target']
    inputtype = opts['inputtype']

    bin_count = 3 if inputtype in ['GeneSym', 'NMID', 'EnsemblID', 'HGNC', 'Position'] else 1
    indices = df.index.tolist()
    bins = np.array_split(indices, bin_count)

    if target == 0:
        target_dir = os.path.join(tempdir, 'full_run')
    else:
        target_dir = os.path.join(tempdir, str(target))
    os.makedirs(target_dir, exist_ok=True)

    # Launch parallel processes
    params = [
        (i, target_dir, df, bin.tolist(), i % opts.get('NUM_GPUs', 1), opts)
        for i, bin in enumerate(bins)
    ]

    mp.set_start_method('spawn', force=True)
    with mp.Pool(bin_count) as pool:
        pool.map(run_deepprime_batch, params)

    # Collect all results
    full_dfs = []
    for root, _, files in os.walk(target_dir):
        for f in files:
            if f.endswith('.pkl'):
                with open(os.path.join(root, f), 'rb') as pf:
                    full_dfs.append(pickle.load(pf))

    if not full_dfs:
        raise RuntimeError("No PAM / pegRNAs Found. Try PE variants with expanded PAM repertoire.")

    # Combine all results
    df_all = pd.concat(full_dfs, ignore_index=True)
    df_sorted = df_all.sort_values(by=f'{pe}_score', ascending=False)

    # Compute z-scores
    if inputtype in ['GeneSym', 'NMID', 'EnsemblID', 'HGNC']:
        scores = df_sorted[f'{pe}_score']
        mean, std = scores.mean(), scores.std()
    else:
        with open(f'{model_dir}/dp_mean.csv') as f:
            mean = float(f.readline().strip())
        with open(f'{model_dir}/dp_std.csv') as f:
            std = float(f.readline().strip())

    df_sorted['Zscore'] = (df_sorted[f'{pe}_score'] - mean) / std

    # Output
    os.makedirs(results_dir, exist_ok=True)
    output_file = os.path.join(results_dir, f'{target}.csv')
    df_sorted.to_csv(output_file, index=False)

    return df_sorted
