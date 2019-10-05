import os
import sys
import numpy
import csv
import math
from matplotlib import rc
from pylab import *
from collections import OrderedDict

matplotlib.rcParams["mathtext.fontset"] = "stix"
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
    r'\usepackage[helvet]{sfmath}', r'\usepackage[utf8]{inputenc}',
    r'\usepackage{arev}', r'\usepackage{siunitx}',
    r'\sisetup{math-micro={\usefont{T1}{phv}{m}{n}\text{µ}}}']

labelFontSize = 15
legendFontSize = 13
lineFontSize = 15

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

colors = [1, 2, 4, 2, 3]

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

fileNames = ['output.txt']
legendTitles = ['pSpatiocyte, A', 'pSpatiocyte, B', 'pSpatiocyte, C']
lines = ['-','-','--','-','-','-','-']
opacity = [1, 1, 1, 1, 1]

volume = 6.21701e-16*1e+18
for f in range(len(fileNames)):
  if (os.path.isfile(path+fileNames[f])):
    data = genfromtxt(path+fileNames[f], delimiter=',', skip_header=1).T
    colSize = len(data)-1
    for i in range(colSize):
      if (i == 0):
        plot(data[0], data[i+1]/volume, ls=lines[f], 
            color=tableau20[colors[i]], label=legendTitles[i],
            linewidth=3, alpha=opacity[f])
      else:
        plot(data[0], data[i+1]/volume, ls=lines[i],
            color=tableau20[colors[i]], label=legendTitles[i],
            linewidth=3, alpha=opacity[f])

num_B = 64000.
num_C = 64000.
duration = 1.01
keff = 20e-19*1e+18
kr = 1.35

from scipy.integrate import odeint


# B + C -> A (keff)
# A -> B + C (kr)
#
# A = keff*B*C - kr*A
# B = -keff*B*C + kr*A
# C = -keff*B*C + kr*A
def f(x, t0, keff, kr):
  return np.array([
    keff * x[1] * x[2] - kr*x[0],
    -keff * x[1] * x[2] + kr*x[0],
    -keff * x[1] * x[2] + kr*x[0]
    ])

init_state = np.array([0, num_B, num_C]) / volume
ode_time = np.logspace(-6,0,100)
ode_result = odeint(f, y0=init_state, t=ode_time, args=(keff, kr))

plot(ode_time, ode_result[:,0], "--", color="k", label="ODE")
plot(ode_time, ode_result[:,1], "--", color="k")
plot(ode_time, ode_result[:,2], "--", color="k")


ax = gca()
handles, labels = ax.get_legend_handles_labels()
leg = legend(loc=(0.1,0.35), labelspacing=0.3, handlelength=2.5,
    handletextpad=0.8, frameon=False)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   

xticks(size=labelFontSize)
yticks(size=labelFontSize)

ax.tick_params(axis='both', which='major', direction='in', length=6, width=1,
    labelsize=lineFontSize)
ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1,
    labelsize=lineFontSize)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
xlabel('Time, $t$ (s)',size=labelFontSize)
ylabel("Concentration (\#\si{\micro}m$^{-3}$)",size=labelFontSize)
xlim(1e-6,5e-1)
plt.xscale('log')
tight_layout()
savefig('reversible_reaction.pdf', format='pdf', dpi=600)#, bbox_inches='tight')
show()

