from scipy import stats
import statsmodels.api as sm

import matplotlib.pyplot as plt
import numpy as np


font_title = { 'color'      : 'k',
                'fontweight' : 'bold',
                'fontsize'   : 16 }

font_normal = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontsize'   : 14 }

font_italic = { 'color'      : 'k',
                'fontweight' : 'normal',
                'fontstyle'  : 'italic',
                'fontsize'   : 12 }



# Simulation sample

#residue = np.array([0.006, 0.579, -0.574, -0.011, 0.206, 0.749, -0.344, 0.209, 0.406, 0.919, -0.094, 0.419])
residue = np.array([-0.2, 0.373, -0.780, -0.217, 0, 0.543, -0.55, 0.003, 0.2, 0.713, -0.3, 0.213])
x = residue

fig = plt.figure(figsize=(10., 7.), dpi=300)

ax1 = fig.add_subplot(121)
res = stats.probplot(x, dist=stats.norm, plot=ax1)
plt.xlabel('Normal', font_normal)
plt.ylabel('Residue', font_normal)
ax1.get_lines()[0].set_markerfacecolor('grey')
ax1.get_lines()[0].set_markeredgecolor('grey')
ax1.get_lines()[1].set_color('k')
ax1.set_title('')

ax2 = fig.add_subplot(122)
n, bins, patches = plt.hist(x, 7, density=False, facecolor='grey')
plt.xlabel('Residue', font_normal)
plt.ylabel('Count', font_normal)

plt.subplots_adjust(wspace=0.8)

plt.tight_layout()

plt.savefig('residue.png')
fig.clear()
plt.close(fig)




#nsample = 100
#np.random.seed(7)

## Sample with normal distribution

#x = np.random.normal(loc=0, scale=1, size=nsample)

#fig = plt.figure(figsize=(10., 7.), dpi=300)

#ax1 = fig.add_subplot(121)
#res = stats.probplot(x, dist=stats.norm, plot=ax1)
#plt.xlabel('Theoretical', font_normal)
#plt.ylabel('Sample', font_normal)
#ax1.get_lines()[0].set_markerfacecolor('grey')
#ax1.get_lines()[0].set_markeredgecolor('grey')
#ax1.get_lines()[1].set_color('k')
#ax1.set_title('')

#ax2 = fig.add_subplot(122)
#n, bins, patches = plt.hist(x, 10, density=False, facecolor='grey')
#plt.xlabel('Sample', font_normal)
#plt.ylabel('Count', font_normal)

#plt.subplots_adjust(wspace=0.4)
#plt.savefig('normal.png')
#fig.clear()
#plt.close(fig)


## Sample with uniform distribution

#x = np.random.uniform(low=0, high=1, size=nsample)

#fig = plt.figure(figsize=(10., 7.), dpi=300)

#ax1 = fig.add_subplot(121)
#res = stats.probplot(x, dist=stats.norm, plot=ax1)
#plt.xlabel('Theoretical', font_italic)
#plt.ylabel('Sample', font_italic)
#ax1.get_lines()[0].set_markerfacecolor('grey')
#ax1.get_lines()[0].set_markeredgecolor('grey')
#ax1.get_lines()[1].set_color('k')
#ax1.set_title('')

#ax2 = fig.add_subplot(122)
#n, bins, patches = plt.hist(x, 10, density=False, facecolor='grey')
#plt.xlabel('Sample', font_italic)
#plt.ylabel('Count', font_italic)


#plt.subplots_adjust(wspace=0.4)
#plt.savefig('uniform.png')
#fig.clear()
#plt.close(fig)


## Sample with beta distribution

#x = np.random.beta(0.5, 0.5, nsample)

#fig = plt.figure(figsize=(10., 7.), dpi=300)

#ax1 = fig.add_subplot(121)
#res = stats.probplot(x, dist=stats.norm, plot=ax1)
#plt.xlabel('Theoretical', font_italic)
#plt.ylabel('Sample', font_italic)
#ax1.get_lines()[0].set_markerfacecolor('grey')
#ax1.get_lines()[0].set_markeredgecolor('grey')
#ax1.get_lines()[1].set_color('k')
#ax1.set_title('')

#ax2 = fig.add_subplot(122)
#n, bins, patches = plt.hist(x, 10, density=False, facecolor='grey')
#plt.xlabel('Sample', font_italic)
#plt.ylabel('Count', font_italic)

#plt.subplots_adjust(wspace=0.4)
#plt.savefig('beta.png')
#fig.clear()
#plt.close(fig)




