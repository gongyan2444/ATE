from astropy.io import fits
import numpy as np
import astropy.wcs as wcs
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib as mpl
#from lmfit import Model
from FITS_tools.hcongrid import hcongrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
from matplotlib.colors import PowerNorm
from matplotlib.gridspec import GridSpec
import aplpy
from astropy.io import fits 
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
from reproject import reproject_interp


mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['patch.linewidth'] = 2
mpl.rcParams['lines.markeredgewidth'] = 2
mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['axes.linewidth'] = 2

mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.pad'] = 5
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.pad'] = 5

mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.minor.size'] = 5

mpl.rcParams['xtick.color'] = 'black'
mpl.rcParams['ytick.color'] = 'black'

fig = plt.figure(figsize=(12, 9))
#fig.subplots_adjust(wspace=0)
#gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
#gs.update(wspace=0)

ref_fits = 'rcw79_wise_4.fits'
hdul_ref = fits.open(ref_fits)
header = hdul_ref[0].header
#header['NAXIS1'] = 2000
#header['NAXIS2'] = 2000
#header['CRPIX1'] = 1000.5
#header['CRPIX2'] = 1000.5
wcs = header


data_r, _ = reproject_interp('G309.5+000IFx_Mosaic_Mom0.fits', header)

data_g, _ = reproject_interp('GLM_30900+0000_mosaic_I4.fits', header)
data_g[np.isnan(data_g)] = np.nanmax(data_g)

data_b, _ = reproject_interp("RCW79_MIPSGAL_024.fits", header)

data_b[np.isnan(data_b)] = np.nanmax(data_b)


fits.writeto('red_reprojected2.fits', data_r, header, overwrite=True)
fits.writeto('green_reprojected2.fits', data_g, header, overwrite=True)
fits.writeto('blue_reprojected2.fits', data_b, header, overwrite=True)


aplpy.make_rgb_image(['red_reprojected2.fits', 'green_reprojected2.fits', 'blue_reprojected2.fits'], 'wise.png',vmin_r=-0.000630765, vmax_r=0.00340979, vmin_g=19.7256, vmax_g=148.736, vmin_b=10.8519, vmax_b=233.272)

gc = aplpy.FITSFigure("rcw79_wise_4.fits", figure=fig, subplot=(1, 2, 1))
gc.show_rgb('wise.png')
gc.show_contour("/media/ygong/512G/2025DomeA/plots/rcw79-ci-sm6-Tmb.fits",levels=4.5+np.arange(6)*1.5, colors='tab:cyan')

# 13 39 54.0 -61 45 00
#gc.recenter(c.ra.deg, c.dec.deg,radius=10./60.)
gc.recenter(204.975, -61.75, radius=15./60.)
gc.tick_labels.set_xformat('hh:mm:ss')
gc.tick_labels.set_yformat('dd:mm')
#gc.tick_labels.set_xrotation(90)
#gc.tick_labels.set_yrotation(0)
gc.axis_labels.set_xtext("R.A. (J2000)")
gc.axis_labels.set_ytext("Dec. (J2000)")

#################




ref_fits = 'rcw120_wise_4.fits'
hdul_ref = fits.open(ref_fits)
header = hdul_ref[0].header
wcs = header

#data_b, header = fits.getdata("RCW120_MG3480p005_024.fits", header=True)

#data_b[np.isnan(data_b)] = np.nanmax(data_b)

data_r, _ = reproject_interp('G348.5+000I_Mosaic_Mom0.fits', header)

data_g, _ = reproject_interp('GLM_34800+0000_mosaic_I4_cutout_10788.fits', header)

data_b, _ = reproject_interp("RCW120_MG3480p005_024.fits", header)

data_b[np.isnan(data_b)] = np.nanmax(data_b)


fits.writeto('red_reprojected.fits', data_r, header, overwrite=True)
fits.writeto('green_reprojected.fits', data_g, header, overwrite=True)
fits.writeto('blue_reprojected.fits', data_b, header, overwrite=True)


aplpy.make_rgb_image(['red_reprojected.fits', 'green_reprojected.fits', 'blue_reprojected.fits'], 'wise.png',vmin_r=-0.000802011, vmax_r=0.0026693, vmin_g=31.6519, vmax_g=273.25, vmin_b=18.9215, vmax_b= 345.398)

gc = aplpy.FITSFigure("rcw120_wise_4.fits", figure=fig, subplot=(1, 2, 2))
gc.show_rgb('wise.png')
gc.show_contour("/media/ygong/512G/2025DomeA/plots/rcw120-ci-m0-6arcmin-Tmb.fits",levels=15+np.arange(6)*3, colors='tab:cyan')

# 17 12 24.0 -38 28 00
#gc.recenter(c.ra.deg, c.dec.deg,radius=10./60.)
gc.recenter(258.100, -38.466667, radius=10./60.)

gc.tick_labels.set_xformat('hh:mm:ss')
gc.tick_labels.set_yformat('dd:mm')
gc.axis_labels.hide_y()
gc.axis_labels.set_xtext("R.A. (J2000)")
#################

plt.savefig("rcw_3col.pdf", bbox_inches='tight',pad_inches=0)
