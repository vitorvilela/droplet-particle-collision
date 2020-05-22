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

# GRID ANALYSIS - DIMENSIONLESS FILM THICKNESS ON PARTICLE

font_normal = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontsize'   : 16 }

font_italic = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontstyle'  : 'italic',
                'fontsize'   : 14 }

# TODO CHECK DIMENSION
dimension = 3

# TODO CHECK BENCHMARK PATH
bench_path = '/home/vitor/Dropbox/Doctorate/defesa/cases/amirfazli-2017/imposed/conservative/benchmarks/'

Benchmark = []
with open(bench_path+'remaining-liquid.csv') as benchmarkFile:
  for line in benchmarkFile:
    benchmark = line.split(',')
    Benchmark.append(map(float, benchmark))
dimensionlessTimeBench, dimensionlessFilmthicknessBench, = zip(*Benchmark)

figure = plt.figure(figsize=(7., 7.), dpi=300)
plt.semilogy(dimensionlessTimeBench, dimensionlessFilmthicknessBench, 'ko')  

# TODO CHECK PATHS OF NUMERICAL RESULTS FILES
pwd = os.getcwd()
#paths = (pwd+'/c3-d3-l9/validation.csv', pwd+'/c3-d3-l10/validation.csv', pwd+'/c3-d3-l11/validation.csv')
paths = (pwd+'/c5-d3-l9/validation.csv', pwd+'/c5-d3-l10/validation.csv', pwd+'/c5-d3-l11/validation.csv', pwd+'/c5-d3-l12/validation.csv')
markers = ('k--', 'k-.', 'k-')
#markers = ('k:', 'k--', 'k-.', 'k-')
legend = ('Banitabaei (2017)', '$D_d/\Delta=30$', '$D_d/\Delta=60$', '$D_d/\Delta=120$')
#legend = ('Banitabaei (2017)', '$D_d/\Delta=30$', '$D_d/\Delta=60$', '$D_d/\Delta=120$', '$D_d/\Delta=240$')

for path, marker in zip(paths, markers):  
  if os.path.exists(path):    
    Numerical = []
    with open(path) as numericalFile:
      for line in numericalFile:
        numerical = line.split(',')
        if numerical[0] != 't':
          floatNumerical = map(float, numerical)
          if floatNumerical[1]<2.5:
            Numerical.append(floatNumerical)
    if dimension == 2:        
      t_n, td_n, vol_n, x_n, y_n, vx_n, vy_n, rl_n, lh_n, lw_n = zip(*Numerical)
    else:
      t_n, td_n, vol_n, x_n, y_n, z_n, vx_n, vy_n, vz_n, rl_n, lh_n, lw_n = zip(*Numerical)

    plt.semilogy(td_n, rl_n, marker)
        
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'$h^*$', font_normal, rotation='vertical')
plt.legend(legend, loc='upper right')

plt.tight_layout()

# TODO CHECK IMAGE FILE NAME
plt.savefig('grid-remaining-liquid.png')
figure.clear()
plt.close(figure)
