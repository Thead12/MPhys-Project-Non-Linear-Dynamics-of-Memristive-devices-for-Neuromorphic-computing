import matplotlib.pyplot as plt
import csv
import textwrap


parameters = []
t = []
x = []
T_prime = []
V = []
V_prime = []

with open("single_neuron_det_results_20251016_184032.csv") as csvfile:
    plots = csv.reader(csvfile, delimiter = ',')

    for row in plots:
        if row[0] == 't':
            parameters = row
            continue
        else:
            t.append(float(row[0]))
            x.append(float(row[1]))
            T_prime.append(float(row[2]))
            V.append(float(row[3]))
            V_prime.append(float(row[4]))


params = ', '.join(parameters[5:])

plt.figure()
plt.plot(t, x)
plt.title('x vs t')
plt.xlabel('t')
plt.ylabel('x')
ax = plt.gca()

wrapped = '\n'.join(textwrap.wrap(params, width=60))
plt.text(0.6, 0.8, f'Parameters: {wrapped}', fontsize=8, ha='center',
                 transform=ax.transAxes, bbox=dict(boxstyle='round,pad=0.3', alpha=0.9))
plt.legend()
plt.show()
