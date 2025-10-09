import matplotlib.pyplot as plt
import numpy as np

i0 = 1


def U_Akther(x):
    #print('x+1', x+1)
    #print('x-1', x-1)
    return (0.05/((x+1)**2)) - 1/(x+1) - 0.5/(x-1) + (1.2*x+0.168)**100


def U_Savelev_stoch(x, i0=i0):
    pot = -( (x - 0.1)**2 + 0.1 ) - 180 * np.exp(-((x - 0.8)**2 / 0.01)) + 0.2 * np.sqrt(10) * i0
    return pot

def U_Savelev_det(x, i0=i0):
    pot = -( (x - 0.1)**2 + 0.1 ) - 180 * np.exp(-((x - 0.8)**2 / 0.01)) + 0.2 * np.sqrt(10) * i0 - 100 * (x ** (100 - 1))
    return pot 

def tiltedPotential(U, q=0.68, V=1, x=0):
    return U - q*V*(x+1)


x1 = -1
x2 = 1

V_ext = 5
q = 0.68

x = np.linspace(start=-1.5, stop=1.5, num=500)

U_Akther = U_Akther(x)
U_Savelev_stoch = U_Savelev_stoch(x)
U_Savelev_det = U_Savelev_det(x)

# Akther Potential plot
plt.figure()
plt.plot(x, U_Akther)
plt.title('Ref.2 Potential')
plt.xlabel('x')
plt.ylabel('U(x)')
plt.xlim(x1, 0.75)
plt.ylim(-5, 2)
plt.show()


"""
# Akther potential plot with V sweep
fig, ax = plt.subplots()
for V in range(0, 35, 5):
    plt.plot(x, tiltedPotential(U_Akther, q=0.68, V=V, x=x), label=f'V={V}')
plt.title('Ref.2 Potential with Varying Voltage bias')
plt.xlabel('x')
plt.ylabel('U(x)+qV(x+1)')
plt.legend()
plt.xlim(x1, x2)
plt.ylim(-35, 1)

#inset axes
axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47], xlim=(-1, -0.75), ylim=(-8, -3))
ax.indicate_inset_zoom(axins, edgecolor="black")
plt.show()
"""

# Akther potential plot with V sweep
fig, ax = plt.subplots()
for V in range(0, 35, 5):
    ax.plot(x, tiltedPotential(U_Akther, q=0.68, V=V, x=x), label=f'V={V}')
ax.set_title('Ref.2 Potential with Varying Voltage bias')
ax.set_xlabel('x')
ax.set_ylabel('U(x)+qV(x+1)')
ax.legend()
ax.set_xlim(x1, 0.75)
ax.set_ylim(-35, 1)

# Inset axes zoomed-in section
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Region to zoom in
x_inset_min, x_inset_max = -1, -0.75
y_inset_min, y_inset_max = -8, -3

axins = inset_axes(ax, width="80%", height="80%", loc='lower left', 
                   bbox_to_anchor=(0.2, 0.05, 0.47, 0.47), bbox_transform=ax.transAxes)
for V in range(0, 35, 5):
    axins.plot(x, tiltedPotential(U_Akther, q=0.68, V=V, x=x))
axins.set_xlim(x_inset_min, x_inset_max)
axins.set_ylim(y_inset_min, y_inset_max)
axins.set_xticks([])
axins.set_yticks([])

mark_inset(ax, axins, loc1=2, loc2=1, fc="none", ec="black")

plt.show()


# Savel'ev potential plot
plt.figure()
plt.plot(x, U_Savelev_stoch, label='Stochastic Potential')
plt.plot(x, U_Savelev_det+30, label='Deterministic Potential with offset + 30')
plt.title('Ref.1 Potentials')
plt.xlabel('x')
plt.ylabel('U(x)')
plt.legend()
plt.xlim(-1.25, 1.25)
plt.ylim(-190, 60)
plt.show()

# Savel'ev potential plot with V sweep
plt.figure()
for V in range(0, 6, 1):
    plt.plot(x, tiltedPotential(U_Savelev_det, q, V, x), label=f'V={V}')
plt.title('Ref.1 Deterministic Potential with Various Voltage bias\'')
plt.xlabel('x')
plt.ylabel('U(x)+qV(x+1)')
plt.legend()
plt.xlim(x1, 0.75)
plt.ylim(-3, 0.55)    
plt.show()


