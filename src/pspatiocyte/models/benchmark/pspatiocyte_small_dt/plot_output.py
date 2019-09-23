import os
import sys
import numpy
import csv
import math
from matplotlib import rc
from pylab import *
from collections import OrderedDict
#uncomment the following to create valid eps (scribus) and svg (inkscape):
#rc('svg', embed_char_paths=True)
#rc('mathtext', fontset=r'stixsans')

labelFontSize = 23
legendFontSize = 18
lineFontSize = 23

matplotlib.rcParams.update({'font.size': labelFontSize})
matplotlib.rcParams['figure.figsize'] = 9.1,7

path, file = os.path.split(os.path.abspath(__file__))
path = path+os.sep

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
    (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
    (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
    (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
    (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)  

colors = [1, 0, 4, 2, 3]

linestyles = OrderedDict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

#fileNames = ['spatiocyte/output.txt', 'pspatiocyte/output.txt', 'smoldyn/output.txt', 'readdy/output.txt', 'readdy_serial/output.txt']
#fileNames = ['spatiocyte/output.txt', 'pspatiocyte_fix/output.txt', 'smoldyn/output.txt', 'readdy/output.txt', 'readdy_serial/output.txt']
fileNames = ['output.txt']
legendTitles = ['Spatiocyte ($\Delta t=0.48\ \mathrm{ms}, T=166\ \mathrm{s}$)','pSpatiocyte (8 cores) ($\Delta t=0.22\ \mathrm{ms},\ T=22\ \mathrm{s}$)','Smoldyn ($\Delta t=1\ \mathrm{ms},\ T=455\ \mathrm{s}$)','Parallel ReaDDy (8 cores) ($\Delta t=1\ \mathrm{ms}, T=549\ \mathrm{s}$)','Serial ReaDDy ($\Delta t=1\ \mathrm{ms}, T=2237\ \mathrm{s}$)']
speciesList = ['E','S','ES','P']
lines = ['-','-','-','-','-','-','-']
opacity = [1, 1, 1, 1, 1]

volume = 909.
for f in range(len(fileNames)):
  if (os.path.isfile(path+fileNames[f])):
    deli = ','
    if (f == 2):
      deli = ' '
    data = genfromtxt(path+fileNames[f], delimiter=deli, skip_header=1).T
    colSize = len(data)-1
    for i in range(colSize):
      if (i == 0):
        plot(data[0], data[i+1]/volume, ls=lines[f], 
            color=tableau20[colors[f]], label=legendTitles[1],
            linewidth=3, alpha=opacity[f])
      else:
        plot(data[0], data[i+1]/volume, ls=lines[f],
            color=tableau20[colors[f]], linewidth=3, alpha=opacity[f])

n_particles_e = 9090
n_particles_s = 90910
duration = 100.0

from scipy.integrate import odeint

def f(x, t0, kf, kr, kcat):
    """
    x: state vector with concentrations of E, S, ES, P
    """
    return np.array([
        -kf * x[0] * x[1] + (kr+kcat)*x[2],
        -kf * x[0] * x[1] + kr * x[2],
        kf * x[0] * x[1]- (kr+kcat)*x[2],
        kcat*x[2]
    ])

init_state = np.array([n_particles_e, n_particles_s, 0., 0.]) / volume
ode_time = np.linspace(0.,duration,100000)
ode_result = odeint(f, y0=init_state, t=ode_time, args=(0.98e-2, 1., 1.))

plot(ode_time, ode_result[:,0], "--", color="k", alpha=.5, label="Mass Action")
plot(ode_time, ode_result[:,1], "--", color="k", alpha=.5)
plot(ode_time, ode_result[:,2], "--", color="k", alpha=.5)
plot(ode_time, ode_result[:,3], "--", color="k", alpha=.5)


annotate('ES', xy=(9, 0),  xycoords='data', xytext=(-29, -10), textcoords='offset points', color='k', size=lineFontSize)

annotate('E', xy=(9, 5),  xycoords='data', xytext=(-27, 15), textcoords='offset points', color='k', size=lineFontSize)

annotate('P', xy=(9, 23),  xycoords='data', xytext=(-27, 12), textcoords='offset points', color='k', size=lineFontSize)

annotate('S', xy=(9, 64),  xycoords='data', xytext=(-27, 12), textcoords='offset points', color='k', size=lineFontSize)


ax = gca()
leg = legend(loc=(0.07,0.35), labelspacing=0.12, handlelength=1.0, handletextpad=0.3, frameon=False)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   
xticks(size=labelFontSize)
yticks(size=labelFontSize)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(axis='both',which='both',direction='in',length=10,width=2)
ax.tick_params(axis='both',which='major',length=10,width=2)
for axis in ['top','bottom','left','right']:
     ax.spines[axis].set_linewidth(2)

xlabel('Time, $t$ (s)',size=labelFontSize)
ylabel("Concentration ($\mathrm{\mu m}^{-3}$)")
xlim(0,10)
tight_layout(pad=0)
show()

