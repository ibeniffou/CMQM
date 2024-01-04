import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from matplotlib.ticker import LinearLocator


# Wigner function of the harmonic oscillator ground state shifted in x by b 

pi = np.pi
h = 1
m = 1/(2*pi)
w = 1
T = 2*pi/w
hbar = h/(2*pi)
a = np.sqrt(hbar/(m*w))
b = 5

def Wigner(x, p, t) :
    W = np.zeros((len(x), len(p)))
    for i in range(len(x)) :
        for j in range(len(p)) :
            W[i][j] = 2/h * np.exp(-a*a/(hbar*hbar)*(p[j]*np.cos(w*t) + m*w*x[i]*np.sin(w*t))**2 - 1/(a*a)*(x[i]*np.cos(w*t) - p[j]/(m*w)*np.sin(w*t) - b)**2)
    return W

# def Wigner(x, p, t) :
#     W = np.zeros((len(x), len(p)))
#     for i in range(len(x)) :
#         for j in range(len(p)) :
#             W[i][j] = 2/h * np.exp(-a*a/(hbar*hbar)*(p[j] + b*hbar/(a*a)*np.sin(w*t))**2 - 1/(a*a)*(x[i]-b*np.cos(w*t))**2)
#     return W

x = np.linspace(-8,8,1000)
p = np.linspace(-1.5,1.5,1000)

W = Wigner(x,p,0) + Wigner(x,p,T/4) + Wigner(x,p,T/2); W = W.T


X, P = np.meshgrid(x, p)

# Plot the surface.
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

surf = ax.plot_surface(X, P, W, cmap=cm.coolwarm,
                    linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(0, 2)
ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.02f}')
ax.set_xlabel("x")
ax.set_ylabel("p")
ax.set_zlabel("W")

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()