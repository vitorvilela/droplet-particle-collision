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


figure = plt.figure(figsize=(6., 6.), dpi=300)
legend = ('$Re=10^3$ $We=10^2$', '$Re=10^3$ $We=10^3$', '$Re=10^4$ $We=10^2$', '$Re=10^4$ $We=10^3$')
markers = ('k-.', 'k-', 'k:', 'k--')

variablesList = []
i = -1
with open('time-area.csv') as csvFile:
  for line in csvFile:
    if line[0] != '#' and line[0] != '%':
      print(line)
      stringVariables = line.split(',')
      floatVariables = map(float, stringVariables)
      variablesList.append(floatVariables)
    if line[0] == '%':  
      print('line[0] == %')
      print(line)
      i += 1
      dimensionlessTime, dimensionlessArea = zip(*variablesList)
      plt.plot(dimensionlessTime, dimensionlessArea, markers[i])
      variablesList.clear()

plt.legend(legend, loc='lower right')
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'$A^*$', font_italic, rotation='vertical')

plt.savefig('time-area.png')
figure.clear()
plt.close(figure)
