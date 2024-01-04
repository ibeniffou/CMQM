import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from matplotlib.ticker import LinearLocator

pi = np.pi
h = 1
m = 1/(2*pi)
w = 1
T = 2*pi/w

hbar = h/(2*pi)
a = np.sqrt(hbar/(m*w))
b = 4

def Wigner(x, p) :
    W = np.zeros((len(x), len(p)))
    for i in range(len(x)) :
        for j in range(len(p)) :
            W[i][j] = 1/(h*(1+np.exp(-b*b/(a*a)))) * np.exp(-a*a*p[j]*p[j]/(hbar*hbar)) * (np.exp(-(x[i]-b)*(x[i]-b)/(a*a)) + np.exp(-(x[i]+b)*(x[i]+b)/(a*a)) + 2*np.exp(-x[i]*x[i]/(a*a))*np.cos(2*b*p[j]/hbar))
    return W


x = np.linspace(-5.5,5.5,1000)
p = np.linspace(-0.45,0.45,1000)

W = Wigner(x,p); W = W.T


X, P = np.meshgrid(x, p)

# Plot the surface.
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

surf = ax.plot_surface(X, P, W, cmap=cm.coolwarm,
                    linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-2, 2)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')
ax.set_xlabel("x")
ax.set_ylabel("p")
ax.set_zlabel("W")

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()