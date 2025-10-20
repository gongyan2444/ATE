from pdrtpy.measurement import Measurement
from pdrtpy.modelset import ModelSet
import pdrtpy.pdrutils as utils
from pdrtpy.tool.lineratiofit import LineRatioFit
from pdrtpy.plot.lineratioplot import LineRatioPlot
from astropy.nddata import StdDevUncertainty
import astropy.units as u
import numpy as np
import corner
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.size'] = 10
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['patch.linewidth'] = 2
mpl.rcParams['lines.markeredgewidth'] = 2
mpl.rcParams['lines.markersize'] = 10
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

mpl.rcParams['figure.figsize'] = (20, 20)



m2 = Measurement(data=21.7,uncertainty = StdDevUncertainty([4.3]),
                 identifier="CI_609",restfreq="492.160651 GHz",unit="K km/s")
m3 = Measurement(data=47.9,uncertainty = StdDevUncertainty([10.2]),
                 identifier="CO_43",restfreq="461.04077 GHz", unit="K km/s")
m4 = Measurement(data=32.3,uncertainty = StdDevUncertainty([3.2]),
                 identifier="CII_158",restfreq="1900.5369 GHz",unit="K km/s")


m5 = Measurement(data=80.8,uncertainty = StdDevUncertainty([8.1]),
                 identifier="CO_10",restfreq="115.2712018 GHz",unit="K km/s")

m6 = Measurement(data=62.8,uncertainty = StdDevUncertainty([6.3]),
                 identifier="CO_32",restfreq="345.7959899 GHz",unit="K km/s")

m7 = Measurement(data=26.0,uncertainty = StdDevUncertainty([2.6]),
                 identifier="13CO_10",restfreq="110.2013541 GHz",unit="K km/s")

m8 = Measurement(data=14.8,uncertainty = StdDevUncertainty([1.5]),
                 identifier="13CO_32",restfreq="330.5879653 GHz",unit="K km/s")

a = [m2,m3,m4,m5,m6,m7,m8]


ms = ModelSet("kt2020",medium="clumpy",mass=1,z=1)
p = LineRatioFit(ms,measurements=a) 

p.run()
#p.fit_result[0].params['density'].min = 1e2
#p.fit_result[0].params['density'].max = 1e6
p.run(method='emcee',steps=2000,nwalkers=200)
#print(p.fit_result.fit_report())
#print(p.fit_result[0].params)
print(p.fit_result[0])
print(f' n = {p.density:.2e}\nI(FUV) = {p.radiation_field:.2e}')
print(f"{utils.toHabing(p.radiation_field):3.2f}")

res = p.fit_result[0]

scale = 1.6e-3 # 1 habing = 1.6E-3 erg s-1 cm-2
# copy the results table
rescopy = deepcopy(res.flatchain)

# scale the radiation_field column of the table.
rescopy['radiation_field'] /= scale
rescopy['density'] /= 1e3
# scale gas density rescopy['density'']/=1e4 

# now copy and scale the "best fit" values where the cross hairs are plotted.
truths=np.array(list(res.params.valuesdict().values()))
truths[1] /=scale
truths[0] /=1e3

#fig = corner.corner(rescopy, labels=["$n$","$G_0$"], truths=truths)

fig = corner.corner(rescopy, bins=20,
                    labels=[r"$n_{\rm gas}~{\rm [10^{3}~cm^{-3}]}$",r"$\chi~{\rm [Habing]}$"],
                    truths=truths, truth_color="#FFA833", show_titles=True,quantiles=(0.16, 0.5, 0.84), label_kwargs={"fontsize": 12}, labelpad=-.1)


fig.savefig("RCW120A-corner.pdf",bbox_inches='tight',pad_inches=0)



#plot = LineRatioPlot(p)

#plot.ratios_on_models(yaxis_unit="Habing",image=True,norm='simple',ncols=2,
#                      figsize=(20,10), meas_color=['red'],label=True,colorbar=True)
# Save the figure to a PNG
#plot.savefig("modelfits.png")


#plot.ratios_on_models(yaxis_unit="Habing", image=True, norm='simple', ncols=2,
#                      meas_color=['red'], label=True, colorbar=True)

#plot.overlay_all_ratios(yaxis_unit="Habing",figsize=(5,5))

#plot.overlay_all_ratios(yaxis_unit="Habing",figsize=(12,12),loc="upper left")

#plot.savefig("overlayratios.pdf")
#plt.tight_layout()  # Adjust layout to prevent overlap
#plt.savefig("combined_plot.png")
