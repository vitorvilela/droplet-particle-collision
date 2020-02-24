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
                'fontsize'   : 16 }

font_italic = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontstyle'  : 'italic',
                'fontsize'   : 14 }

# TODO CHECK PATHS OF NUMERICAL RESULTS FILES
pwd = os.getcwd()
#paths = (pwd+'/second-results/simulation.csv', pwd+'/c3-d3-l10/simulation.csv', pwd+'/c3-d3-l11/simulation.csv')
#paths = ('/home/vitor/Dropbox/Doctorate/defesa/cases/mitra-2018/simulations/validation/second-results/simulation.csv')
markers = ('k-')
marker = 'k-'
#markers = ('k--', 'k-.', 'k-')

figure = plt.figure(figsize=(8., 8.), dpi=300)
#legend = ('$l=9$', '$l=10$', '$l=11$')
legend = '$l=8$'

#for path, marker in zip(paths, markers):  
  #if os.path.exists(path):    
Log = []
    #with open(path) as logFile:
with open('/home/vitor/Dropbox/Doctorate/defesa/cases/mitra-2018/simulations/validation/second-results/simulation.csv') as logFile:
  for line in logFile:
    log = line.split(',')
    if log[0] != 't':
      floatLog = map(float, log)
      if floatLog[1] < 2.5:
        Log.append(floatLog)
t, td, dt, mgpi, mgui, grid_tn, perf_t, perf_speed = zip(*Log)
plt.semilogy(td[10:], grid_tn[10:], marker)

plt.grid()
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'Grid size', font_italic, rotation='vertical')
#plt.legend(legend, loc='best')

# TODO CHECK IMAGE FILE NAME
plt.savefig('c3-gridSize.png')
figure.clear()
plt.close(figure)

