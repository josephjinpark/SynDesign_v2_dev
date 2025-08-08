import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale
from scipy import stats


def generate_zscore_plot(zscore, all_zscores, label, output_path):
    """
    Create a scatter plot showing the Z-score percentile rank.
    """
    percents = minmax_scale(all_zscores)
    yvals = [p * 100 for p in percents]

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.scatter(all_zscores, yvals, c='black', s=25, label='All pegRNAs')
    ax.axvline(zscore, color='red', label='Selected pegRNA')

    ax.set_title(f'Z-score Distribution: {label}')
    ax.set_xlabel('Z-Score')
    ax.set_ylabel('Percentile Rank')
    ax.set_ylim([0, 100])
    ax.legend()

    fig.savefig(output_path)
    plt.close(fig)


def create_html_output_table(df_sorted: pd.DataFrame, pe: str, top_n: int, results_dir: str, base_name: str):
    """
    Save top-N pegRNA predictions and plots to output directory.
    """
    html_cols = [
        'ID', 'GuideSeq', 'PBSlen', 'RTlen', 'RHA_len', f'{pe}_score', 'Zscore', 'percentile', 'plotpath',
        'wtseq', 'Edited74_On', 'Edited_wNote'
    ]
    plots_dir = os.path.join(results_dir, 'matplot')
    os.makedirs(plots_dir, exist_ok=True)

    df_top = df_sorted.head(top_n).copy()
    df_top['percentile'] = df_top['Zscore'].apply(lambda z: f"{int(stats.percentileofscore(df_sorted['Zscore'], z))}th")
    df_top['plotpath'] = df_top['ID'].apply(lambda i: f'{plots_dir}/{i}.png')

    # Generate rank plots
    for _, row in df_top.iterrows():
        generate_zscore_plot(
            zscore=row['Zscore'],
            all_zscores=df_sorted['Zscore'],
            label=pe,
            output_path=row['plotpath']
        )

    # Save output tables
    html_table = df_top[html_cols]
    html_table.to_csv(os.path.join(results_dir, f'{base_name}_html.csv'), index=False)
    df_sorted.to_csv(os.path.join(results_dir, f'{base_name}_fullout.csv'), index=False)

    return html_table
