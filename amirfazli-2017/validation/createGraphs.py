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
                'fontsize'   : 24 }

font_italic = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontstyle'  : 'italic',
                'fontsize'   : 20 }


# TODO CHECK IT !!!
bench_path = '../../benchmarks/'

pwd = os.getcwd()
subcase = pwd.split('/')[-1]
dimension = int(subcase[4])



# GRID ANALYSIS - DIMENSIONLESS FILM THICKNESS ON PARTICLE

Benchmark = []
with open(bench_path+'remaining-liquid.csv') as benchmarkFile:
  for line in benchmarkFile:
    benchmark = line.split(',')
    Benchmark.append(map(float, benchmark))
t_b, var_b, = zip(*Benchmark)


figure = plt.figure(figsize=(8., 8.), dpi=300)
plt.semilogy(t_b, var_b, 'ko')  


# TODO CHECK IT !!!
paths = ['validation.csv']
marks = ['k-']
legend = ('Banitabaei (2017)', '$D_p/\Delta=120$')

for path, mark in zip(paths, marks):  
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

    plt.semilogy(td_n, rl_n, mark)
    
    
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'$h^*$', font_normal, rotation='vertical')
plt.legend(legend, loc='upper right')

plt.savefig(subcase+'-gridAnalysis-filmThickness.png')
figure.clear()
plt.close(figure)

# GRID ANALYSIS - DIMENSIONLESS FILM THICKNESS ON PARTICLE




# SIMULATION TRACKING - TIME x GRID 

Log = []
with open('simulation.csv') as logFile:
  for line in logFile:
    log = line.split(',')
    if log[0] != 't':
      floatLog = map(float, log)
      Log.append(floatLog)
t, td, dt, mgpi, mgui, grid_tn, perf_t, perf_speed = zip(*Log)


figure = plt.figure(figsize=(8., 8.), dpi=300)

plt.semilogy(td[10:], grid_tn[10:], 'k-')
plt.grid()
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'Grid size', font_italic, rotation='vertical')

plt.savefig(subcase+'-simulationTracking-grid.png')
figure.clear()
plt.close(figure)

# SIMULATION TRACKING - TIME x GRID 




# ENERGY ANALYSIS - FLOW ENERGY AND INTERFACIAL AREA

SIGMA = 0.0728

Numerical = []
with open('energy.csv') as numericalFile:
  for line in numericalFile:
    numerical = line.split(',')
    if numerical[0] != 't':
      floatNumerical = map(float, numerical)
      Numerical.append(floatNumerical)
t, td, g, kl, kg, st, a, ad = zip(*Numerical)

# I am not sure if the very first [0] has a valid (conservative) value, so it in set to [1]
# ekl0 = kl[1] 
# ekg0 = kg[1]
# est0 = st[1]
# eg0 = g[1]

e0 = st[1] + kl[1]

E0 = []
SumEg = []
Feed = []
Flow = []
Sink = []
Ae = []

sum_eg = 0.

for i in range(1, len(g)):
  
  sum_eg += g[i]
  
  E0.append(e0)
  SumEg.append(sum_eg)
  
  Feed.append(e0 + sum_eg)
  Flow.append(kl[i] + kg[i] + st[i])
  
  Sink.append(abs(Feed[i-1]-Flow[i-1]))
  Ae.append((kl[i] + kg[i]) / SIGMA )
  
  
  
legend = ('$\sum {E_g}$', '$Ek_l$', '$Ek_g$', '$E\sigma$')

figure = plt.figure(figsize=(9., 9.), dpi=300)
plt.subplots_adjust(hspace=0.4, wspace=0.4)

plt.subplot(221)
plt.semilogy(td[1:], SumEg, 'k-')
plt.ylabel(legend[0]+'    [ J ]', font_normal)
plt.xlabel(r'$t^*$', font_normal)

plt.subplot(222)
plt.semilogy(td[1:], kl[1:], 'k-')
plt.ylabel(legend[1]+'    [ J ]', font_normal)
plt.xlabel(r'$t^*$', font_normal)

plt.subplot(223)
plt.semilogy(td[1:], kg[1:], 'k-')
plt.ylabel(legend[2]+'    [ J ]', font_normal)
plt.xlabel(r'$t^*$', font_normal)

plt.subplot(224)
plt.semilogy(td[1:], st[1:], 'k-')
plt.ylabel(legend[3]+'    [ J ]', font_normal)
plt.xlabel(r'$t^*$', font_normal)

plt.savefig(subcase+'-energies.png')
figure.clear()
plt.close(figure)




legend = ('Feed', 'Flow', 'Sink')

figure = plt.figure(figsize=(8., 8.), dpi=300)

plt.semilogy(td[1:], Feed, 'k-') #, color='#000000')
plt.semilogy(td[1:], Flow, 'k--') #, color='#778899')
plt.semilogy(td[1:], Sink, 'k-.') #, color='#808080')
#plt.semilogy(td[1:], Flow, '-', color='#DCDCDC')


plt.legend(legend, loc='lower right')
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'Energy    [ J ]', font_italic, rotation='vertical')

plt.savefig(subcase+'-energy.png')
figure.clear()
plt.close(figure)


figure = plt.figure(figsize=(8., 7.), dpi=300)

plt.plot(td, ad, 'k-')
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'$A^*$', font_normal, rotation='vertical')
plt.savefig(subcase+'-normalizedArea.png')

figure.clear()
plt.close(figure)





figure = plt.figure(figsize=(8., 7.), dpi=300)

plt.plot(td, a, 'k-')
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'Area    [ $m^2$ ]', font_italic, rotation='vertical')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
plt.savefig(subcase+'-interfacialArea.png')

figure.clear()
plt.close(figure)




a_diff = [abs(a_vof-a_energy) for a_vof, a_energy in zip(a[1:], Ae)]

legend = ('$A_{VOF}$', '$A_{Energy}$', '$\delta$')

figure = plt.figure(figsize=(8., 8.), dpi=300)

plt.semilogy(td[1:], a[1:], '-', color='#000000')
plt.semilogy(td[1:], Ae, '-', color='#778899')
plt.semilogy(td[1:], a_diff, '-', color='#808080')
#plt.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))

plt.legend(legend, loc='upper right')
plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'Area    [ $m^2$ ]', font_normal, rotation='vertical')

plt.savefig(subcase+'-area.png')
figure.clear()
plt.close(figure)


# ENERGY ANALYSIS - FLOW ENERGY AND INTERFACIAL AREA





# VALIDATION - DROPLET WIDTH, HEIGHT AND FILM TICKNESS ON PARTICLE

Numerical = []
with open('validation.csv') as numericalFile:
  for line in numericalFile:
    numerical = line.split(',')
    if numerical[0] != 't':
      floatNumerical = map(float, numerical)
      if floatNumerical[1]<3.:
        Numerical.append(floatNumerical)
if dimension == 2:        
  t_n, td_n, vol_n, x_n, y_n, vx_n, vy_n, rl_n, lh_n, lw_n = zip(*Numerical)
else:
  t_n, td_n, vol_n, x_n, y_n, z_n, vx_n, vy_n, vz_n, rl_n, lh_n, lw_n = zip(*Numerical)


Benchmark = []
with open(bench_path+'lamella-width-we1146.csv') as benchmarkFile:
  for line in benchmarkFile:
    benchmark = line.split(',')
    floatBenchmark = map(float, benchmark)
    if floatBenchmark[0]<3.:
      Benchmark.append(floatBenchmark)
tb, bench = zip(*Benchmark)


figure = plt.figure(figsize=(8., 8.), dpi=300)

plt.plot(tb, bench, 'ko')
plt.plot(td_n, lw_n, 'k-')

plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'$D^*$', font_normal, rotation='vertical')

plt.savefig(subcase+'-lamellaWidth.png')
figure.clear()
plt.close(figure)


Benchmark = []
with open(bench_path+'lamella-height-we1146.csv') as benchmarkFile:
  for line in benchmarkFile:
    benchmark = line.split(',')
    floatBenchmark = map(float, benchmark)
    if floatBenchmark[0]<3.:
      Benchmark.append(floatBenchmark)
tb, bench = zip(*Benchmark)


figure = plt.figure(figsize=(8., 8.), dpi=300)

plt.plot(tb, bench, 'ko')
plt.plot(td_n, lh_n, 'k-')

plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'$H^*$', font_normal, rotation='vertical')

plt.savefig(subcase+'-lamellaHeight.png')
figure.clear()
plt.close(figure)


Benchmark = []
with open(bench_path+'remaining-liquid.csv') as benchmarkFile:
  for line in benchmarkFile:
    benchmark = line.split(',')
    floatBenchmark = map(float, benchmark)
    if floatBenchmark[0]<3.: 
      Benchmark.append(floatBenchmark)
tb, bench = zip(*Benchmark)

figure = plt.figure(figsize=(8., 8.), dpi=300)

plt.plot(tb, bench, 'ko')
plt.plot(td_n, rl_n, 'k-')

plt.xlabel(r'$t^*$', font_normal)
plt.ylabel(r'$h^*$', font_normal, rotation='vertical')

plt.savefig(subcase+'-filmThickness.png')
figure.clear()
plt.close(figure)


# VALIDATION - DROPLET WIDTH, HEIGHT AND FILM TICKNESS ON PARTICLE




# Options:

#ax.yaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
#ax.xaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
#plt.xlim(X.min()*1.1, X.max()*1.1)
#plt.ylim(C.min()*1.1, C.max()*1.1)
#plt.xticks( [-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
#plt.yticks([-1, 0, +1])
#fig, ax = plt.subplots()
#plt.loglog(t, grid_tn, basex=2)
#plt.legend(('Banitabaei (2017)', 'Current work'), loc='upper right')
