import matplotlib.pyplot as plt
import csv
import textwrap
import numpy as np
import re


parameters = []
t = []
x = []
T_prime = []
V = []
V_prime = []

def format_param_value(num_str, sig=3):
    try:
        x=float(num_str)
    except ValueError:
        return num_str
    if abs(x-round(x)) < 1e-12:
        return str(int(round(x)))
    if 1e-4 <= abs(x) < 1e4:
        return np.format_float_positional(x, precision=sig, unique=False, fractional=False, trim='k')
    return np.format_float_scientific(x, precision=sig-1, unique=False, trim='k')

with open("data\Vext_2.csv") as csvfile:
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

params = re.sub(r'=(\-?\d+\.?\d*(?:[eE][+-]?\d+)?)',
                lambda m: '=' + format_param_value(m.group(1)),
                params)
print(params)

plt.figure()
plt.plot(t, x)
plt.title('x vs t: V_ext=2')
plt.xlabel('t')
plt.ylabel('x')
ax = plt.gca()

wrapped = '\n'.join(textwrap.wrap(params, width=40))
plt.text(0.7, 0.2, f'Parameters: {wrapped}', fontsize=8, ha='center',
                 transform=ax.transAxes, bbox=dict(boxstyle='round,pad=0.3', alpha=0.7))
plt.legend()
plt.show()
