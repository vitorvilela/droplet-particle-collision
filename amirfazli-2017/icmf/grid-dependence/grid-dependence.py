import numpy as np
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

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

markers = ('ko-', 'ks-')

l =  np.array([0, 10, 11, 12]) # 10, 11, 12
area_for_t1_95 = np.array([0., 5.24, 5.41, 5.50])
area_for_t3_38 = np.array([0., 6.62, 6.98, 7.03])

figure = plt.figure(figsize=(6., 6.), dpi=300)

ax = figure.gca()
legend = ('$t^*=1.95$', '$t^*=2.38$')

opacity = 0.4
#bar_width = 0.35

plt.xlabel(r'$l$', font_normal)
plt.ylabel(r'$A^*$', font_italic, rotation='vertical')
#plt.xticks(range(len(l)),(l[0], l[1], l[2]), rotation=0)

#plt.axes([9,9,10,10,11,11])

#plt.bar(l, +area_for_t1_95, alpha=opacity+0.4, color='k', label='t^*=1.95')
#plt.bar(l, -area_for_t3_38, alpha=opacity, color='k', label='t^*=2.38')

#for x, y in zip(l, area_for_t1_95):
  #plt.text(x+0.4, y+0.05, '%.2f' % area_for_t1_95, ha='center', va= 'top')
    
#for x, y in zip(l, area_for_t3_38):
  #plt.text(x+0.4, -y-0.05, '%.2f' % y, ha='center', va= 'bottom')

#bar1 =  plt.bar(np.arange(len(area_for_t1_95)) + bar_width, area_for_t1_95, bar_width, align='center', alpha=opacity+0.4, color='k', label='t^*=1.95')
#bar2 = plt.bar(range(len(area_for_t3_38)), area_for_t3_38, bar_width, align='center', alpha=opacity, color='k', label='t^*=2.38')


plt.plot(l[1:], area_for_t1_95[1:], markers[0], markersize=12)
plt.plot(l[1:], area_for_t3_38[1:], markers[1], markersize=12)

##plt.xlim(-.5, n)
#plt.xticks([])
##plt.ylim(-1.25, +1.25) 
#plt.yticks([])


plt.ylim(ymin=5)
plt.ylim(ymax=7.25)

ax.xaxis.set_major_locator(MaxNLocator(integer=True))


plt.legend(legend, loc='best', fontsize=20)

plt.tight_layout()

plt.savefig('grid-dependence-2.png')
figure.clear()
plt.close(figure)
