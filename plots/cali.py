from astropy.io import fits
import numpy as np

dat1, hd1 = fits.getdata("rcw120-co43-m0-6arcmin.fits", header=True)
dat2, hd2 = fits.getdata("rcw120-ci-m0-6arcmin.fits", header=True)

fits.writeto("rcw120-co43-m0-6arcmin-Tmb.fits", dat1/0.35, hd1, overwrite=True)
fits.writeto("rcw120-ci-m0-6arcmin-Tmb.fits", dat2/0.35, hd2, overwrite=True)

dat3, hd3 = fits.getdata("rcw79-co43-corr3-crop-sm6.fits", header=True)
dat4, hd4 = fits.getdata("rcw79-ci-sm6.fits", header=True)

fits.writeto("rcw79-co43-corr3-crop-sm6-Tmb.fits", dat3/0.35, hd3, overwrite=True)
fits.writeto("rcw79-ci-sm6-Tmb.fits", dat4/0.35, hd4, overwrite=True)
