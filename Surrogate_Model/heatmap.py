import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("LH_samples/latin_sample_combined.csv")

def plot_raw_lhs_heatmap(df, metric_col_name):
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]
    z = df[metric_col_name]

    plt.figure(figsize=(10, 8))
    plt.tricontourf(x, y, z, levels=50, cmap='plasma')
    plt.scatter(x, y, color='white', s=2, alpha=0.3) 
    plt.colorbar(label=metric_col_name)
    plt.title(f'Raw Data Interpolation: {metric_col_name}')
    plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_metric_distribution(csv_path, metric_index):
    # Load data
    df = pd.read_csv(csv_path)
    
    metric_data = df.iloc[:, 2 + metric_index]
    metric_name = df.columns[2 + metric_index]

    plt.figure(figsize=(10, 6))
    sns.histplot(metric_data, kde=True, color='royalblue', bins=30)
    
    plt.title(f'Distribution of {metric_name}', fontsize=14)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.3)
    plt.show()

# is_phase_lock,width,CV1_at_target_freq,CV2_at_tartget_freq,centre_of_phase_locking,left_depth,righ_depth,depth
for i in range(8):
    plot_metric_distribution('LH_samples/latin_sample_combined.csv', metric_index=i)

# plot_raw_lhs_heatmap(df, "depth")