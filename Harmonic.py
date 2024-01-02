import matplotlib.pyplot as plt
import numpy as np
# from scipy.constants import hbar

from matplotlib import cm
from matplotlib.ticker import LinearLocator

h = 1
hbar = h/(2*np.pi)

def Wigner0(x, p, a) :
    return 2/h * np.exp( -a*a*p*p/(hbar*hbar) - x*x/(a*a) )

def Wigner1(x, p, a) :
    return 2/h * (-1 + 2*(a*p/hbar)*(a*p/hbar) + 2*(x/a)*(x/a)) * np.exp(-a*a*p*p/(hbar*hbar) - x*x/(a*a))
    

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
x = np.arange(-2, 2, 0.01)
p = np.arange(-0.2, 0.2, 0.01)
x, p = np.meshgrid(x, p)

a = 1
W = Wigner1(x, p, a)
print(W.shape)

# Plot the surface.
surf = ax.plot_surface(x, p, W, cmap=cm.coolwarm,
                    linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-2, 2)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')
ax.set_xlabel("x")
ax.set_ylabel("p")
ax.set_zlabel("W")
ax.set_title("a = "+str(a))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()