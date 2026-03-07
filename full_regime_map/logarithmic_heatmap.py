import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap

def load_data(csv_filename):
    df = pd.read_csv(csv_filename)
    V_ext = df['V_ext'].values
    R_ext = df['R_ext'].values
    x_max = df['x_max'].values
    x_min = df['x_min'].values
    spike_count = df['spike_count'].values
    burst_freq = df ['burst_freq'].values
    avg_ibi = df['avg_ibi'].values
    intra_burst_freq = df['intra_burst_freq'].values

def block_heatmap(X,Y, Z, Z_lims, thresholds=(None, None), hatching=None):
    
    xi = np.logspace(np.log10(X.min()+1e-5), np.log10(X.max()), 200)
    yi = np.linspace(Y.min(), Y.max(), 500)
    X_grid, Y_grid = np.meshgrid(xi, yi)

    block_map = np.zeros(len(Z[0]))
    block_id = 0


    for i in range(0,len(thresholds[0][:])):
        if thresholds[1][i] == '>':
            thresholded_Z = Z[i] > thresholds[0][i]
        else:
            thresholded_Z = Z[i] < thresholds[0][i]    
        thresholded_Z.astype(int)

        block_map[thresholded_Z] = block_id + 1
        block_id+=1
                
    Z_grid = griddata((X, Y), block_map, (X_grid, Y_grid), method='nearest')

    fig, ax = plt.subplots(figsize=(10, 6))

    colours = [
        '#008000', 
        "#B34100", 
        '#FF00FF', 
        '#FFD700',
        '#FF0000',         
        '#000000',
        "#00C3FF", 
        '#0000FF', 
        ]

    custom_cmap = ListedColormap(colours)
    mesh = ax.pcolormesh(X_grid, Y_grid, Z_grid, vmin=Z_lims[0], vmax=Z_lims[1], shading='auto', cmap=custom_cmap)

    ax.set_xscale('log')
    ax.set_xlim([1.0, None])

    ax.set_xlabel('$log_{10}(V_{ext})$')
    ax.set_ylabel('$R_{ext}$')
    ax.set_title('Codimensional Regime Map: Up Sweep')

    plt.tight_layout()

    # add hatching to spiking regimes
    for i in range(len(thresholds[0])):
        if thresholds[1][i] == '>':
            mask = Z[i] > thresholds[0][i]
        else:
            mask = Z[i] < thresholds[0][i]

        thresholded_Z = mask.astype(int)

        hatch = hatching[i]
        if hatch:
            Z_grid = griddata((X, Y), thresholded_Z, (X_grid, Y_grid), method='nearest')
            if np.any(mask):
                plt.contourf(X_grid, Y_grid, Z_grid, levels=[0.5, 1.5], 
                             colors='none', hatches=[hatch], alpha=0)
    
                

    plt.show()

def create_heatmap(X,Y, Z, Z_lims):
    
    xi = np.logspace(np.log10(X.min()+1e-5), np.log10(X.max()), 200)
    yi = np.linspace(Y.min(), Y.max(), 200)
    X_grid, Y_grid = np.meshgrid(xi, yi)
        
    Z_grid = griddata((X, Y), Z, (X_grid, Y_grid), method='nearest')

    fig, ax = plt.subplots(figsize=(10, 6))

    mesh = ax.pcolormesh(X_grid, Y_grid, Z_grid, vmin=Z_lims[0], vmax=Z_lims[1], shading='auto', cmap='viridis')

    ax.set_xscale('log')

    plt.colorbar(mesh, label='x_{max}')
    ax.set_xlim([1.0, None])
    ax.set_xlabel('$V_{ext}$')
    ax.set_ylabel('$R_{ext}$')
    ax.set_title('Bi-parametric Heatmap: Up Sweep')

    plt.tight_layout()
    #plt.savefig('heatmap_output_down.png', dpi=300)
    plt.show()

df = pd.read_csv('Master_Aggregated_Grid_interpolated.csv')
#df = pd.read_csv('spiking_extension_up.csv', skiprows=1)

V_ext = df['V_ext'].values
R_ext = df['R_ext'].values
x_max = df['x_max'].values
x_min = df['x_min'].values
spike_count = df['spike_count'].values
burst_freq = df ['burst_freq'].values
avg_ibi = df['avg_ibi'].values
intra_burst_freq = df['intra_burst_freq'].values

block_heatmap(V_ext, R_ext, (x_max, x_max, x_max, x_max, spike_count, intra_burst_freq),
                Z_lims=(None, None),
                thresholds=([[0.0, 0.75, -0.75, -10.0, 4.0, 0.5],
                             ['>', '>', '<', '<', '>', '>',  '>']]),
                hatching=(False, False, False, False, '\|', 'O\|\|'))


"""
block_heatmap(V_ext, R_ext, (x_max, x_max),
              Z_lims=(None, None),
              thresholds=[[0.0, -0.5],
                          ['>',  '<']],
              replace_with=(1, 2)) 
"""

create_heatmap(V_ext, R_ext, x_max, (-1.0, 1.0))

#create_heatmap(V_ext, R_ext, spike_count/500, (0.0, 0.4))