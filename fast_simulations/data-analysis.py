import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def plot_normalized_detection_vs_mu_for_shapes(df, fixed_n,
                                               fixed_n_volume=None,
                                               fixed_n_targets=None,
                                               fixed_D=None,
                                               fixed_probability=None,
                                               save_path=None):
    base_filter = (df['n_walkers'] == fixed_n)
    if fixed_n_volume is not None:
        base_filter &= (df['n_volume'] == fixed_n_volume)
    if fixed_n_targets is not None:
        base_filter &= (df['n_targets'] == fixed_n_targets)
    if fixed_D is not None:
        base_filter &= (df['D'] == fixed_D)
    if fixed_probability is not None:
        base_filter &= (df['probability'] == fixed_probability)

    subset = df[base_filter].copy()
    if subset.empty:
        print("Nessun dato corrispondente ai filtri specificati.")
        return

    group_cols = ['mu', 'TargetShape', 'n_walkers']
    if fixed_n_volume is not None: group_cols.append('n_volume')
    if fixed_n_targets is not None: group_cols.append('n_targets')
    if fixed_D is not None: group_cols.append('D')
    if fixed_probability is not None: group_cols.append('probability')

    stats = subset.groupby(group_cols)['detection_time'].mean().reset_index()
    stats = stats.rename(columns={'detection_time': 'detection_time_mean'})

    merge_keys = [col for col in group_cols if col != 'TargetShape']
    ball_ref = stats[stats['TargetShape'] == 'Ball'][merge_keys + ['detection_time_mean']]
    ball_ref = ball_ref.rename(columns={'detection_time_mean': 'ball_detection_time'})

    merged = stats.merge(ball_ref, on=merge_keys, how='left')
    merged['normalized_mean'] = merged['detection_time_mean'] / merged['ball_detection_time']
    merged = merged.dropna(subset=['normalized_mean'])

    plt.figure(figsize=(10, 6))
    for shape, data in merged.groupby('TargetShape'):
        data = data.sort_values('mu')
        plt.plot(data['mu'], data['normalized_mean'], label=shape, marker='o')

    plt.xticks(sorted(merged['mu'].unique()))

    titolo = f'Normalized Detection Time vs μ (n_walkers={fixed_n}'
    if fixed_n_volume is not None:
        side = round(fixed_n_volume ** (1/3))
        titolo += f', side={side}'
    if fixed_n_targets is not None: titolo += f', n_targets={fixed_n_targets}'
    if fixed_D is not None: titolo += f', D={fixed_D}'
    if fixed_probability is not None: titolo += f', p={fixed_probability}'
    titolo += ')'

    plt.axhline(1.0, color='gray', linestyle='--', label='Ball = 1')
    plt.title(titolo)
    plt.xlabel('μ')
    plt.ylabel('Detection Time / Ball Detection Time')
    plt.legend(title='Shape')
    plt.grid(True)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()

def plot_detection_vs_mu_for_shapes(df, fixed_n,
                                    fixed_n_volume=None,
                                    fixed_n_targets=None,
                                    fixed_D=None,
                                    fixed_probability=None,
                                    save_path=None):
    base_filter = (df['n_walkers'] == fixed_n)
    if fixed_n_volume is not None:
        base_filter &= (df['n_volume'] == fixed_n_volume)
    if fixed_n_targets is not None:
        base_filter &= (df['n_targets'] == fixed_n_targets)
    if fixed_D is not None:
        base_filter &= (df['D'] == fixed_D)
    if fixed_probability is not None:
        base_filter &= (df['probability'] == fixed_probability)

    subset = df[base_filter].copy()
    if subset.empty:
        print("Nessun dato corrispondente ai filtri specificati.")
        return

    group_cols = ['mu', 'TargetShape', 'n_walkers']
    if fixed_n_volume is not None: group_cols.append('n_volume')
    if fixed_n_targets is not None: group_cols.append('n_targets')
    if fixed_D is not None: group_cols.append('D')
    if fixed_probability is not None: group_cols.append('probability')

    stats = subset.groupby(group_cols)['detection_time'].agg(['mean', 'std']).reset_index()
    stats = stats.rename(columns={'mean': 'detection_time_mean', 'std': 'detection_time_std'})

    plt.figure(figsize=(10, 6))
    for shape, data in stats.groupby('TargetShape'):
        data = data.sort_values('mu')
        if shape == 'Line':
            # Plot con barre di errore per Line
            plt.plot(data['mu'], data['detection_time_mean'],
                         label=shape, marker='o', linestyle='-', alpha=0.8)
        else:
            # Plot con banda di deviazione standard per gli altri
            plt.plot(data['mu'], data['detection_time_mean'], label=shape, marker='o')
            plt.fill_between(data['mu'],
                             data['detection_time_mean'] - data['detection_time_std'],
                             data['detection_time_mean'] + data['detection_time_std'],
                             alpha=0.2)

        # Aggiungi valore del primo punto vicino al marker (solo Ball, Line, Disk)
        if shape in ['Ball', 'Line', 'Disk']:
            first_point = data.iloc[0]
            x = first_point['mu']
            y = first_point['detection_time_mean']
            plt.text(x, y, f'{y*1e-6:.2f}M', fontsize=9, fontweight='bold',
                     ha='left', va='bottom', color='black')

    plt.xticks(sorted(stats['mu'].unique()))

    # Imposta tick y personalizzati in milioni
    y_min = stats['detection_time_mean'].min()
    y_max = stats['detection_time_mean'].max()
    y_min_tick = np.floor(y_min / 1e6) * 1e6
    y_max_tick = np.ceil(y_max / 1e6) * 1e6
    yticks = np.arange(y_min_tick, y_max_tick + 2e6, 2e6)
    plt.yticks(yticks)

    # Formatter asse y in milioni
    def millions_formatter(x, pos):
        return f'{x*1e-6:.1f}M'
    plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(millions_formatter))
    plt.ylabel('Detection Time (×10⁶)')

    titolo = f'Detection Time vs μ (n_walkers={fixed_n}'
    if fixed_n_volume is not None:
        side = round(fixed_n_volume ** (1/3))
        titolo += f', side={side}'
    if fixed_n_targets is not None: titolo += f', n_targets={fixed_n_targets}'
    if fixed_D is not None: titolo += f', D={fixed_D}'
    if fixed_probability is not None: titolo += f', p={fixed_probability}'
    titolo += ')'

    plt.title(titolo)
    plt.xlabel('μ')
    plt.legend(title='Shape')
    plt.grid(True)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()

def plot_detection_vs_walkers(df, diametri, distance_values, fixed_mu, fixed_n_volume, fixed_probability, save_path=None):
    # Filtro i dati
    subset = df[
        (df['TargetShape'] == 'Ball') &
        (df['D'].isin(diametri)) &
        (df['fixed_target_dist'].isin(distance_values)) &
        (df['n_targets'] == 1) &
        (df['mu'] == fixed_mu) &
        (df['n_volume'] == fixed_n_volume) &
        (df['probability'] == fixed_probability)
    ]

    # Valori effettivi di n_walkers
    x_ticks = sorted(subset['n_walkers'].unique())

    # FacetGrid: y non condivisa
    g = sns.FacetGrid(subset, col='D', row='fixed_target_dist', margin_titles=True, height=3.5, sharey=False)
    g.map_dataframe(sns.lineplot, x='n_walkers', y='detection_time', marker='o')

    # Label e titoli
    g.set_axis_labels("Number of Walkers", "Detection Time")
    g.set_titles(row_template='Distance = {row_name}', col_template='D = {col_name}')
    plt.subplots_adjust(top=0.9)

    # Imposta tick asse x
    for ax in g.axes.flat:
        ax.set_xticks(x_ticks)

    side = round(fixed_n_volume ** (1/3))
    g.fig.suptitle(f'Detection Time vs Number of Walkers — shape=Ball, μ={fixed_mu}, side={side}, p={fixed_probability}, n_targets=1')

    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()

def plot_detection_vs_walkers_reverse(df, diametri, n_walkers_values, fixed_mu, fixed_n_volume, fixed_probability, save_path=None):
    # Filtra dati
    subset = df[
        (df['TargetShape'] == 'Ball') &
        (df['D'].isin(diametri)) &
        (df['n_walkers'].isin(n_walkers_values)) &
        (df['n_targets'] == 1) &
        (df['mu'] == fixed_mu) &
        (df['n_volume'] == fixed_n_volume) &
        (df['probability'] == fixed_probability)
    ]

    # Tick distanza x asse
    x_ticks = sorted(subset['fixed_target_dist'].unique())

    # FacetGrid con righe = D, colonne = n_walkers
    g = sns.FacetGrid(subset, row='D', col='n_walkers', margin_titles=True, height=3.5, sharey=False)
    g.map_dataframe(sns.lineplot, x='fixed_target_dist', y='detection_time', marker='o')

    # Labels e titoli
    g.set_axis_labels("Distance from Target", "Detection Time")
    g.set_titles(row_template='D = {row_name}', col_template='n_walkers = {col_name}')
    plt.subplots_adjust(top=0.9)

    # Setta ticks asse x solo per valori presenti
    for ax in g.axes.flat:
        ax.set_xticks(x_ticks)

    side = round(fixed_n_volume ** (1/3))
    g.fig.suptitle(f'Detection Time vs Distance — shape=Ball, μ={fixed_mu}, side={side}, p={fixed_probability}, n_targets=1')

    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()

def plot_detection_vs_mu_by_probability(df,
                                       fixed_D,
                                       fixed_n_walkers,
                                       fixed_n_volume,
                                       fixed_n_targets,
                                       save_path=None):
    valid_ps = [0.1, 1]
    valid_shapes = ['Ball', 'Line']
    
    base_filter = (df['D'] == fixed_D) & \
                  (df['n_walkers'] == fixed_n_walkers) & \
                  (df['n_volume'] == fixed_n_volume) & \
                  (df['n_targets'] == fixed_n_targets) & \
                  (df['probability'].isin(valid_ps)) & \
                  (df['TargetShape'].isin(valid_shapes))
    
    subset = df[base_filter]
    
    # Raggruppamento e aggregazione
    grouped = subset.groupby(['mu', 'TargetShape', 'probability'], as_index=False)['detection_time'].agg(['mean', 'std']).reset_index()
    grouped.rename(columns={'mean': 'detection_time_mean', 'std': 'detection_time_std'}, inplace=True)

    # Plotting
    plt.figure(figsize=(10, 6))
    for shape in valid_shapes:
        for p in valid_ps:
            data = grouped[(grouped['TargetShape'] == shape) & (grouped['probability'] == p)].sort_values('mu')
            if data.empty:
                continue
            label = f"{shape}, p={p}"
            plt.plot(data['mu'], data['detection_time_mean'], marker='o', label=label)
            plt.fill_between(data['mu'],
                             data['detection_time_mean'] - data['detection_time_std'],
                             data['detection_time_mean'] + data['detection_time_std'],
                             alpha=0.2)
    
    plt.xlabel('μ')
    plt.ylabel('Detection Time')
    side = round(fixed_n_volume ** (1/3))
    plt.title(f'Detection Time vs μ — D={fixed_D}, side={side}, n_agents={fixed_n_walkers}, n_targets={fixed_n_targets}')
    plt.legend(title='Shape, p')
    plt.grid(True)
    plt.xticks(sorted(grouped['mu'].unique()))
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()

def plot_detection_vs_k_targets(df,
                                 fixed_mu,
                                 fixed_n_walkers,
                                 fixed_n_volume,
                                 fixed_probability,
                                 fixed_shape,
                                 D_values,
                                 save_path=None):
    plt.figure(figsize=(10, 6))
    
    # Ciclo sui diversi valori di D
    for D in D_values:
        k_values = df['n_targets'].unique()
        detection_means = []
        k_clean = []
        
        for k in [1,2,4,8, 16]:
            d_k = D / k
            mask = (
                (df['mu'] == fixed_mu) &
                (df['n_walkers'] == fixed_n_walkers) &
                (df['n_volume'] == fixed_n_volume) &
                (df['probability'] == fixed_probability) &
                (df['TargetShape'] == fixed_shape) &
                (df['n_targets'] == k) &
                (df['D'] == d_k)
            )
            subset = df[mask]
            if not subset.empty:
                detection_time_mean = subset['detection_time'].mean()
                detection_means.append(detection_time_mean)
                k_clean.append(k)
        
        if detection_means:
            plt.plot(k_clean, detection_means, marker='o', label=f'D = {D}')
    
    plt.xlabel('Number of Targets (k)')
    plt.ylabel('Detection Time')
    side = round(fixed_n_volume ** (1/3))
    plt.title(f'Detection Time vs Number of Targets — μ={fixed_mu}, shape={fixed_shape}, side={side}, p={fixed_probability}, n={fixed_n_walkers}')
    plt.legend(title='Original D')
    plt.grid(True)
    plt.xticks(sorted(df['n_targets'].unique()))
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()

def plot_detection_vs_k_targets_with_line(df,
                                 fixed_mu,
                                 fixed_n_walkers,
                                 fixed_n_volume,
                                 fixed_probability,
                                 fixed_shape,
                                 D_values,
                                 save_path=None):
    plt.figure(figsize=(10, 6))
    
    # Generatore di colori
    colors = plt.cm.get_cmap('tab10', len(D_values))

    for i, D in enumerate(D_values):
        k_values = [1, 2, 4, 8, 16]
        detection_means = []
        k_clean = []

        for k in k_values:
            d_k = D / k
            mask = (
                (df['mu'] == fixed_mu) &
                (df['n_walkers'] == fixed_n_walkers) &
                (df['n_volume'] == fixed_n_volume) &
                (df['probability'] == fixed_probability) &
                (df['TargetShape'] == fixed_shape) &
                (df['n_targets'] == k) &
                (df['D'] == d_k)
            )
            subset = df[mask]
            if not subset.empty:
                detection_time_mean = subset['detection_time'].mean()
                detection_means.append(detection_time_mean)
                k_clean.append(k)

        if detection_means:
            color = colors(i)
            plt.plot(k_clean, detection_means, marker='o', label=f'{fixed_shape}, D = {D}', color=color)

            # Linea di riferimento con stesso colore
            mask_line = (
                (df['mu'] == fixed_mu) &
                (df['n_walkers'] == fixed_n_walkers) &
                (df['n_volume'] == fixed_n_volume) &
                (df['probability'] == fixed_probability) &
                (df['TargetShape'] == 'Line') &
                (df['n_targets'] == 1) &
                (df['D'] == D)
            )
            line_subset = df[mask_line]
            if not line_subset.empty:
                line_time = line_subset['detection_time'].mean()
                plt.hlines(line_time, min(k_values), max(k_values),
                           colors=[color], linestyles='--', linewidth=1.5,
                           label=None)  # niente legenda per la linea

    plt.xlabel('Number of Targets (k)')
    plt.ylabel('Detection Time')
    side = round(fixed_n_volume ** (1/3))
    plt.title(f'Detection Time vs Number of Targets — μ={fixed_mu}, shape={fixed_shape}, side={side}, p={fixed_probability}, n={fixed_n_walkers}')
    plt.legend(title='Shape & D')
    plt.grid(True)
    plt.xticks(k_values)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()

def plot_detection_vs_diameter(df,
                                                fixed_mu,
                                                fixed_n_volume,
                                                fixed_probability,
                                                fixed_n_targets,
                                                fixed_n_walkers,
                                                save_path=None):
    shapes = df['TargetShape'].unique()
    plt.figure(figsize=(10, 6))

    plotted_d_values = set()

    # Plotta detection time medio per ogni shape
    for shape in shapes:
        subset = df[
            (df['mu'] == fixed_mu) &
            (df['n_volume'] == fixed_n_volume) &
            (df['probability'] == fixed_probability) &
            (df['n_targets'] == fixed_n_targets) &
            (df['n_walkers'] == fixed_n_walkers) &
            (df['TargetShape'] == shape) &
            (df['D'] != 0)
        ]

        if subset.empty:
            continue

        grouped = subset.groupby('D')['detection_time'].mean().reset_index()
        grouped = grouped.sort_values('D')
        plotted_d_values.update(grouped['D'].values)

        plt.plot(grouped['D'], grouped['detection_time'], marker='o', label=shape)

    # Calcola e aggiungi curve teoriche
    d_values = sorted(plotted_d_values)
    volume_over_d = [fixed_n_volume / d for d in d_values]
    volume_over_d2 = [fixed_n_volume / (d**2) for d in d_values]

    plt.plot(d_values, volume_over_d, linestyle='--', color='gray', label='Volume / D')
    plt.plot(d_values, volume_over_d2, linestyle=':', color='black', label='Volume / D²')

    # Etichette e stile
    plt.xlabel('Diametro D')
    plt.ylabel('Detection Time')
    side = round(fixed_n_volume ** (1/3))
    plt.title(f'Detection Time vs D — μ={fixed_mu}, side={side}, n={fixed_n_walkers}, p={fixed_probability}, n_targets={fixed_n_targets}')
    plt.xticks(d_values)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
    plt.close()



if __name__ == "__main__":
    # Supponiamo che i tuoi file si trovino in una cartella e abbiano estensione .csv
    file_paths = glob.glob('./results/*.csv')  # Sostituisci con il path corretto

    # Carica tutti i CSV e concatenali in un unico DataFrame
    df_list = [pd.read_csv(file) for file in file_paths]
    df = pd.concat(df_list, ignore_index=True)
    output_dir = './plots'
    #plot_normalized_detection_vs_mu_for_shapes(df, fixed_n_targets=1, fixed_n=1, fixed_D=16, fixed_n_volume=16777216, save_path=f'{output_dir}/normalized_detection_vs_mu_shapes_n1.png', fixed_probability=1)
    #plot_detection_vs_mu_for_shapes(df, fixed_n_targets=1, fixed_n=1, fixed_D=16, fixed_n_volume=16777216, save_path=f'{output_dir}/unnormalized_detection_vs_mu_shapes_n1.png', fixed_probability=1)
    # ripeti dopo aver rieseguito le simulazioni
    #plot_detection_vs_walkers(df, diametri=[2, 4, 8, 16, 32, 64], distance_values=[2, 4, 8, 16, 32, 64],save_path=f'{output_dir}/detection_distance.png', fixed_n_volume=2097152, fixed_probability=1, fixed_mu = 2)
    plot_detection_vs_walkers_reverse(df, diametri=[0,1,2, 4, 8, 16, 32, 64], n_walkers_values=[1, 2, 4, 8, 16, 32],save_path=f'{output_dir}/detection_distance_reverse2.png', fixed_n_volume=2097152, fixed_probability=1, fixed_mu = 2)
    #plot_detection_vs_mu_by_probability(df, fixed_D=16, fixed_n_walkers=1, fixed_n_volume=2097152, fixed_n_targets=1, save_path=f'{output_dir}/detection_vs_mu_by_shape_and_probability.png')
    #plot_detection_vs_k_targets(df, fixed_mu=2, fixed_n_walkers=1,fixed_n_volume=16777216,fixed_probability=1,fixed_shape="Ball",D_values=[8,16, 32, 64, 128],save_path=f'{output_dir}/multi_target.png')
    #plot_detection_vs_k_targets_with_line(df, fixed_mu=2, fixed_n_walkers=1,fixed_n_volume=16777216,fixed_probability=1,fixed_shape="Ball",D_values=[8,16, 32, 64, 128],save_path=f'{output_dir}/multi_target_with_line.png')
    #plot_detection_vs_diameter(df,fixed_mu=2,fixed_n_volume=16777216,fixed_probability=1,fixed_n_targets=1, fixed_n_walkers=1, save_path=f'{output_dir}/detection_vs_diameter_fixed_n.png')