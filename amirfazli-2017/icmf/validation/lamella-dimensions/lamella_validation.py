import numpy as np
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.ticker import FormatStrFormatter

import subprocess

from pylab import *

import os
import os.path


font_normal = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontsize'   : 16 }

font_italic = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontstyle'  : 'italic',
                'fontsize'   : 14 }

t = np.array([0.87, 1.08, 1.73, 1.95, 2.38, 2.81])

H_exp = np.array([1.23, 1.35, 1.69, 1.82, 2.13, 2.51])
H_exp_error = np.array([0.1227, 0.1345, 0.1685, 0.1818, 0.2129, 0.2513])
H_num = np.array([1.13, 1.25, 1.55, 1.77, 2.24, 2.63])

Db_exp = np.array([1.51, 1.63, 1.97, 2.04, 2.17, 2.41])
Db_exp_error = np.array([0.1508, 0.1626, 0.1966, 0.2040, 0.2173, 0.2409])
Db_num = np.array([1.28, 1.45, 1.99, 2.19, 2.39, 2.53])


# markers: 'k-', 'k--', 'k-.', 'k:'

legend = ('Simulation', 'Banitabaei (2017)')

figure = plt.figure(figsize=(6., 6.), dpi=300)

#plt.grid()
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'$H^*$', font_italic, rotation='vertical')
plt.plot(t, H_num, 'k-')
plt.errorbar(t, H_exp, H_exp_error, fmt='.', color='black', elinewidth=0.5, capsize=5, capthick=0.5)
plt.legend(legend, loc='upper left')

plt.tight_layout()

# TODO CHECK IMAGE NAME
plt.savefig('lamella-height.png')
figure.clear()
plt.close(figure)



figure = plt.figure(figsize=(6., 6.), dpi=300)

#plt.grid()
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'$D_b^*$', font_italic, rotation='vertical')
plt.plot(t, Db_num, 'k-')
plt.errorbar(t, Db_exp, Db_exp_error, fmt='.', color='black', elinewidth=0.5, capsize=5, capthick=0.5)
plt.legend(legend, loc='upper left')

plt.tight_layout()

# TODO CHECK IMAGE NAME
plt.savefig('lamella-diameter.png')
figure.clear()
plt.close(figure)

