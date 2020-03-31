'''
======================
3D surface (color map)
======================

Demonstrates plotting a 3D surface colored with the coolwarm color map.
The surface is made opaque by using antialiased=False.

Also demonstrates using the LinearLocator and custom formatting for the
z axis tick labels.
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter

class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex
    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)
    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))


fig1 = plt.figure(1,figsize=(8, 4))
ax = plt.subplot(1,1,1)

# Make data.
X = np.arange(0, 5*np.pi+0.1, 0.1)
Y = np.arange(0, (np.pi+0.1), 0.1)
X, Y = np.meshgrid(X, Y)
Z = (np.sin(X)**2 + np.sin(Y)**2 - 2*np.sin(X)*np.sin(Y))

plane=ax.pcolor(X, Y, Z)
ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
ax.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
ax.yaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
fig1.colorbar(plane, ax=ax)


fig2 = plt.figure(2)
ax1 = fig2.gca(projection='3d')

# Plot the surface.
surf = ax1.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)

# Customize the z axis.
#ax1.set_zlim(-1.01, 2.01)
ax1.zaxis.set_major_locator(LinearLocator(10))
ax1.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax1.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax1.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
ax1.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
ax1.yaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
ax1.yaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
ax1.yaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
# Add a color bar which maps values to colors.
fig2.colorbar(surf, shrink=0.5, aspect=5)

plt.show()


