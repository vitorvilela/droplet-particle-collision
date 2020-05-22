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

# ENERGY ANALYSIS - FEED, FLOW AND LOST ENERGY

font_normal = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontsize'   : 16 }

font_italic = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontstyle'  : 'italic',
                'fontsize'   : 14 }

# TODO CHECK PATHS OF NUMERICAL RESULTS FILES
pwd = os.getcwd()
#paths = (pwd+'/c5-d3-l10/energy.csv', pwd+'/c5-d3-l11/energy.csv', pwd+'/c5-d3-l12/energy.csv')
#markers = ('k:', 'k--', 'k-.', 'k-')




path = 'energy.csv'
marker = 'k'


figure1 = plt.figure(1, figsize=(7., 7.), dpi=300)
figure10 = plt.figure(10, figsize=(7., 7.), dpi=300)
figure11 = plt.figure(11, figsize=(7., 7.), dpi=300)

figure2 = plt.figure(2, figsize=(8., 8.), dpi=300)
figure3 = plt.figure(3, figsize=(8., 8.), dpi=300)
figure4 = plt.figure(4, figsize=(8., 8.), dpi=300)
figure5 = plt.figure(5, figsize=(8., 8.), dpi=300)



SIGMA = 0.0728
#SIGMA = 0.0825

#for path, marker in zip(paths, markers):  
  #if os.path.exists(path): 
    
    #subcase = path.split('/')[-2]
    
Numerical = []
with open(path) as numericalFile:
  for line in numericalFile:
    numerical = line.split(',')
    if numerical[0] != 't':
      floatNumerical = map(float, numerical)
      if floatNumerical[1] < 2.5:
        Numerical.append(floatNumerical)
t, td, g, kl, kg, st, a, ad = zip(*Numerical)

e0 = st[1] + kl[1]

E0 = []
SumEg = []
Feed = []
Flow = []
Lost = []
Ae = []

sum_eg = 0.

for i in range(1, len(g)):      
  sum_eg += g[i]      
  E0.append(e0)
  SumEg.append(sum_eg)      
  Feed.append(e0 + sum_eg)
  Flow.append(kl[i] + kg[i] + st[i])      
  Lost.append(abs(Feed[i-1]-Flow[i-1]))
  Ae.append((kl[i] + kg[i]) / SIGMA )
   
   
plt.figure(1) 
plt.semilogy(td[1:], Lost, marker)
 
plt.figure(10) 
plt.plot(td[1:], Feed, marker) 
    
plt.figure(11)
plt.semilogy(td[1:], Flow, marker) 
   
   
plt.figure(2) 
plt.semilogy(td[1:], SumEg, marker)
   
plt.figure(3) 
plt.semilogy(td[1:], kl[1:], marker)
  
plt.figure(4) 
plt.semilogy(td[1:], kg[1:], marker)
  
plt.figure(5) 
plt.semilogy(td[1:], st[1:], marker)
   


# TODO CHECK LEG  
#legend = ('$l=9$', '$l=10$', '$l=11$', '$l=12$')
#legend = '$l=10$'


# TODO CHECK IMAGE FILE NAME
plt.figure(1) 
plt.ylabel(r'Lost    [ J ]', font_normal, rotation='vertical')
plt.xlabel(r'$t^*$', font_normal)
#plt.legend(legend, loc='lower right')
plt.tight_layout()
plt.savefig('c5-lostEnergy.png')
figure1.clear()
plt.close(figure1)

# TODO CHECK IMAGE FILE NAME
plt.figure(10) 
plt.ylabel(r'Feed    [ J ]', font_normal, rotation='vertical')
plt.ticklabel_format(axis='y', style='sci', scilimits=(1, 1))
plt.xlabel(r'$t^*$', font_normal)
#plt.legend(legend, loc='best')
plt.tight_layout()
plt.savefig('c5-feedEnergy.png')
figure10.clear()
plt.close(figure10)

# TODO CHECK IMAGE FILE NAME
plt.figure(11) 
plt.ylabel(r'Flow    [ J ]', font_normal, rotation='vertical')
plt.xlabel(r'$t^*$', font_normal)
#plt.legend(legend, loc='upper right')
plt.tight_layout()
plt.savefig('c5-flowEnergy.png')
figure11.clear()
plt.close(figure11)



# TODO CHECK IMAGE FILE NAME
plt.figure(2) 
plt.ylabel('$\sum{_t}$ ${{E_G}}$    [ J ]', font_normal)
plt.xlabel(r'$t^*$', font_normal)
#plt.legend(legend, loc='lower right')
plt.tight_layout()
plt.savefig('c5-sumEG.png')
figure2.clear()
plt.close(figure2)

# TODO CHECK IMAGE FILE NAME
plt.figure(3) 
plt.ylabel('$Ek_l$    [ J ]', font_normal)
plt.xlabel(r'$t^*$', font_normal)
#plt.legend(legend, loc='upper right')
plt.tight_layout()
plt.savefig('c5-Ekl.png')
figure3.clear()
plt.close(figure3)

# TODO CHECK IMAGE FILE NAME
plt.figure(4) 
plt.ylabel('$Ek_g$    [ J ]', font_normal)
plt.xlabel(r'$t^*$', font_normal)
#plt.legend(legend, loc='upper left')
plt.tight_layout()
plt.savefig('c5-Ekg.png')
figure4.clear()
plt.close(figure4)

# TODO CHECK IMAGE FILE NAME
plt.figure(5) 
plt.ylabel('$E\gamma$    [ J ]', font_normal)
plt.xlabel(r'$t^*$', font_normal)
#plt.legend(legend, loc='upper left')
plt.tight_layout()
plt.savefig('c5-Es.png')
figure5.clear()
plt.close(figure5)

