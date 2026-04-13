import SelectivityLib as sl
import numpy as np
import csv
import time

from scipy.stats import qmc


def SimulationOutput_to_list(result):
    list_output = []
    list_output.append(result.is_phase_lock)
    list_output.append(result.width)
    list_output.append(result.CV1_at_target_freq)
    list_output.append(result.CV2_at_target_freq)
    list_output.append(result.centre_of_phase_locking)
    list_output.append(result.left_depth)
    list_output.append(result.right_depth)
    list_output.append(result.depth)

    return list_output
        
p = sl.Params()
p.V_AC = 3.0
p.target_freq = 0.064
p.dt = 1e-3
p.warmUp_cycles = 20.0
p.num_cycles = 100.0
p.frequency_step = 0.0005

sim = sl.Simulator(p, verbose=False)

# linear grid
# v_space = np.linspace(5.0, 7.0, 3)
# r_space = np.linspace(400, 600, 3)
# batch_inputs = []
# for v in v_space:
#     for r in r_space:
#         batch_inputs.append(sl.SimulationInput(v, r))

batch_inputs = []
# latin hypercube sampling
l_bounds = [3.5, 150.0]
h_bounds = [8.0, 1000.0]
sampler = qmc.LatinHypercube(d=2)
samples = sampler.random(n=256)
samples = qmc.scale(samples, l_bounds, h_bounds)
for sample in samples:
    batch_inputs.append(sl.SimulationInput(sample[0], sample[1]))

print(f"Starting parallel execution of {len(batch_inputs)} simulations...")

start_time = time.perf_counter()
results = sim.run_batch(batch_inputs)
end_time = time.perf_counter()

duration = end_time - start_time
print(f"Batch simulation completed in {duration:.2f} seconds.")

header = ["V_DC","R_ext","is_phase_lock","width","CV1_at_target_freq","CV2_at_tartget_freq","centre_of_phase_locking","left_depth","righ_depth","depth"]
print("Saving to CSV file latin_sample.csv")
with open('latin_sample.csv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(header)
    for i, res in enumerate(results):
        v_val = batch_inputs[i].V_DC
        r_val = batch_inputs[i].R_ext

        list_res = SimulationOutput_to_list(res)

        writer.writerow([v_val,r_val] + list_res)

print("Results saved to latin_sample.csv") 