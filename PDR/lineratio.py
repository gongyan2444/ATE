import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl


# initialize plot settings
mpl.rcParams['font.size'] = 16
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['patch.linewidth'] = 2
#mpl.rcParams['lines.markeredgewidth'] = 2
#mpl.rcParams['lines.markersize'] = 2
mpl.rcParams['axes.linewidth'] = 2

mpl.rcParams['xtick.direction'] = 'in' 
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.pad'] = 5
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.pad'] = 5

mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.width'] = 2
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.minor.size'] = 5

mpl.rcParams['figure.figsize'] = (18, 10)

#        RCW79A      RCW79B      RCW120A   RCW120B 
#[CII]   44.3+0.5    34.1+0.5   20.9+0.8   37.1+0.8
#        26.5+0.2    19.3+0.2   32.3+0.2   26.5+0.2
#        70.8+0.5    53.4+0.5   53.2+0.8   63.6+0.8
#[CI]    8.8+1.4     6.3+1.0    21.7+1.1   20.3+0.9
#
#CO(4-3) 49.9+3.4    21.8+1.2   47.9+1.8   58.6+2.6

def eratio(x, dx, y, dy):
    er = np.sqrt(dx**2/y**2+ x**2*dy**2/y**4)
    return er


dat1 = np.array([70.8, 53.4, 53.2, 63.6])
dat1e = np.array([0.5, 0.5, 0.8, 0.8])

dat2 = np.array([8.8, 6.3, 21.7, 20.9])
dat2e = np.array([1.4, 1.0, 1.1, 0.9])

dat3 = np.array([49.9, 6.3, 21.7, 20.9])
dat3e = np.array([3.4, 1.2, 1.8, 2.6])

r1 = dat1/dat2
er1 = eratio(dat1, dat1e, dat2, dat2e)
r2 = dat2/dat3
er2 = eratio(dat2, dat2e, dat3, dat3e)

g0  = np.array([172, 253, 29, 31])
g0e1 = np.array([64, 151, 11, 10])
g0e2 = np.array([48, 100, 9, 9])

# Define the labels for each point
labels = ['RCW79A', 'RCW79B', 'RCW120A', 'RCW120B']


fig, axes = plt.subplots(1, 3, figsize=(18, 5))

sc = axes[0].scatter(r1, r2, c=g0, s=100, cmap='inferno', edgecolor='k', zorder=2)
axes[0].errorbar(r1, r2, xerr=er1, yerr=er2, fmt='none', ecolor='gray', capsize=4, zorder=1)
axes[0].set_xlabel(r'$I_{\rm [CII]}/I_{\rm [CI]}$')
axes[0].set_ylabel(r'$I_{\rm [CI]}/I_{\rm CO}$')
# --- For axes[0]: r1 vs r2 ---
axes[0].annotate(labels[0], (r1[0], r2[0]), xytext=(5, 5), textcoords='offset points', fontsize=12)
axes[0].annotate(labels[1], (r1[1], r2[1]), xytext=(5, 5), textcoords='offset points', fontsize=12)
axes[0].annotate(labels[2], (r1[2], r2[2]), xytext=(-10, 25), textcoords='offset points', fontsize=12)
axes[0].annotate(labels[3], (r1[3], r2[3]), xytext=(10, -5), textcoords='offset points', fontsize=12)

cbax = plt.axes([0.126, 0.9, 0.225, 0.02])
cbar = plt.colorbar(sc, ax=axes[0],cax=cbax, label=r'$G_0$ [Habing]',  orientation="horizontal")
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.xaxis.set_label_position('top')

g0_err = [g0e2, g0e1]  # [lower, upper]

sc1 = axes[1].scatter(r1, g0, c=dat3, cmap='inferno', s=100, edgecolor='k', zorder=2)
axes[1].errorbar(r1, g0, xerr=er1, yerr=g0_err, fmt='none', ecolor='gray', capsize=4, zorder=1)
axes[1].set_xlabel(r'$I_{\rm [CII]}/I_{\rm [CI]}$')
axes[1].set_ylabel(r'$G_0$')
axes[1].annotate(labels[0], (r1[0], g0[0]), xytext=(5, -25), textcoords='offset points', fontsize=12)
axes[1].annotate(labels[1], (r1[1], g0[1]), xytext=(5, 5), textcoords='offset points', fontsize=12)
axes[1].annotate(labels[2], (r1[2], g0[2]), xytext=(-10, 15), textcoords='offset points', fontsize=12)
axes[1].annotate(labels[3], (r1[3], g0[3]), xytext=(10, -5), textcoords='offset points', fontsize=12)

cbax1 = plt.axes([0.4, 0.9, 0.225, 0.02])
cbar1 = plt.colorbar(sc1, ax=axes[1], cax=cbax1, label=r'$I_{\rm CO}$ (K.km/s)', orientation="horizontal")
cbar1.ax.xaxis.set_ticks_position('top')
cbar1.ax.xaxis.set_label_position('top')

# ---- r2 vs. g0, color by dat1 ----
sc2 = axes[2].scatter(r2, g0, c=dat1, cmap='inferno', s=100, edgecolor='k', zorder=2)
axes[2].errorbar(r2, g0, xerr=er2, yerr=g0_err, fmt='none', ecolor='gray', capsize=4, zorder=1)
axes[2].annotate(labels[0], (r2[0], g0[0]), xytext=(5, 5), textcoords='offset points', fontsize=12)
axes[2].annotate(labels[1], (r2[1], g0[1]), xytext=(5, 5), textcoords='offset points', fontsize=12)
axes[2].annotate(labels[2], (r2[2], g0[2]), xytext=(5, -12), textcoords='offset points', fontsize=12)
axes[2].annotate(labels[3], (r2[3], g0[3]), xytext=(5, 5), textcoords='offset points', fontsize=12)

axes[2].set_xlabel(r'$I_{\rm [CI]}/I_{\rm [CO]}$')
axes[2].set_ylabel(r'$G_0$')
cbax2 = plt.axes([0.673, 0.9, 0.225, 0.02])
cbar2 = plt.colorbar(sc2, ax=axes[2], cax=cbax2, label=r'$I_{\rm [CII]}$ (K.km/s)', orientation="horizontal")
cbar2.ax.xaxis.set_ticks_position('top')
cbar2.ax.xaxis.set_label_position('top')

plt.savefig('corr.pdf', pad_inches=0,bbox_inches='tight')

#plt.tight_layout()
#plt.show()
