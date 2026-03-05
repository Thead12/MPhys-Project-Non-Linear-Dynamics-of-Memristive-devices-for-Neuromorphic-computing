import pandas as pd
import numpy as np
import os
from scipy.interpolate import griddata

FILE_LIST = [
    'spike_count_up.csv',
    'burst_mapping_v2_up.csv',
    'burst_extension_up.csv',
    'spiking_extension_up.csv'    
]

def aggregate_data(files):
    processed_data = []
    
    for f in files:
        if not os.path.exists(f):
            print(f"Skipping: {f} (File not found)")
            continue
            
        # Load the file
        df = pd.read_csv(f, skiprows=1)
        print(df.columns.values.tolist())

        processed_data.append(df)

    if not processed_data:
        print("No valid data found to aggregate.")
        return
    
    aggregate_map = pd.concat(processed_data, )
    aggregate_map = aggregate_map.replace(np.nan, 0)
    aggregate_map = aggregate_map.replace('-nan(ind)', 0)
    
    all_v = aggregate_map['V_ext'].values
    all_r = aggregate_map['R_ext'].values    

    #v_grid = np.logspace(np.log10(all_v.min()+1), np.log10(all_v.max()), 200)

    v_grid = np.linspace(all_v.min(), all_v.max(), 200)
    r_grid = np.linspace(all_r.min(), all_r.max(), 200)
    V_MESH, R_MESH = np.meshgrid(v_grid, r_grid)

    interpolated_map = pd.DataFrame({'V_ext': V_MESH.flatten(), 'R_ext': R_MESH.flatten()})

    for column in aggregate_map.iloc[:, 2:]:
        print('Column: ', column)

        points = (all_v, all_r)
        values = aggregate_map[column]
        xi = (V_MESH, R_MESH)

        z_interp = griddata(points=points, values=values, xi=xi, method='nearest', fill_value=0)
        flat_z = z_interp.flatten()
        
        interpolated_map[column] = flat_z


    print(interpolated_map.head)

    interpolated_map.to_csv("Master_Aggregated_Grid_interpolated.csv", index=False)
    print(f"Successfully aggregated {len(processed_data)} files into 'Master_Aggregated_Grid.csv'")

if __name__ == "__main__":
    aggregate_data(FILE_LIST)