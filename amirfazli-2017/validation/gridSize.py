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

# GRID SIZE

font_normal = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontsize'   : 24 }

font_italic = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontstyle'  : 'italic',
                'fontsize'   : 20 }

# TODO CHECK PATHS OF NUMERICAL RESULTS FILES
pwd = os.getcwd()
paths = (pwd+'/c5-d3-l10/simulation.csv', pwd+'/c5-d3-l11/simulation.csv', pwd+'/c5-d3-l12/simulation.csv')
markers = ('k:', 'k--', 'k-')

figure = plt.figure(figsize=(6., 6.), dpi=300)
legend = ('$l=10$', '$l=11$', '$l=12$')

for path, marker in zip(paths, markers):  
  if os.path.exists(path):    
    Log = []
    with open(path) as logFile:
      for line in logFile:
        log = line.split(',')
        if log[0] != 't':
          floatLog = map(float, log)
          if floatLog[1] < 2.5:
            Log.append(floatLog)
    t, td, dt, mgpi, mgui, grid_tn, perf_t, perf_speed = zip(*Log)
    plt.semilogy(td[10:], grid_tn[10:], marker, linewidth = 3)

#plt.grid()
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'Grid size', font_italic, rotation='vertical')
plt.legend(legend, loc='upper left', fontsize=20)

plt.tight_layout()

# TODO CHECK IMAGE FILE NAME
plt.savefig('gridSize.png', transparent=False)
figure.clear()
plt.close(figure)

