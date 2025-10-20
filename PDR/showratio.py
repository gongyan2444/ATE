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


a = [m2,m3,m4]


#ms = ModelSet("wk2020",z=1)
ms = ModelSet("kt2020",medium="clumpy",mass=1,z=1)
p = LineRatioFit(ms,measurements=a) 

p.run()


plot = LineRatioPlot(p)

plot.ratios_on_models(yaxis_unit="Habing",image=True,norm='log',ncols=2,
                      figsize=(15,6), meas_color=['cyan'],label=True,colorbar=True)
# Save the figure to a PNG
plot.savefig("modelfits.pdf")


#plot.ratios_on_models(yaxis_unit="Habing", image=True, norm='simple', ncols=2,
#                      meas_color=['red'], label=True, colorbar=True)

#plot.overlay_all_ratios(yaxis_unit="Habing",figsize=(5,5))

#plot.overlay_all_ratios(yaxis_unit="Habing",figsize=(12,12),loc="upper left")

#plot.savefig("overlayratios.pdf")
#plt.tight_layout()  # Adjust layout to prevent overlap
#plt.savefig("combined_plot.png")
