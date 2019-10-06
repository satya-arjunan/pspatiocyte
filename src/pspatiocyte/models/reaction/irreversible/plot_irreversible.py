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

legendTitles = ['pSpatiocyte, A', 'pSpatiocyte, B', 'pSpatiocyte, C']
lines = ['-','-','--','-','-','-','-']
opacity = [1, 1, 1, 1, 1]
ks = [0.135, 1.35, 13.5]

volume = 6.21701e-16*1e+18

for f, k in enumerate(ks): 
  dirname = 'k_'+str(k)
  print(dirname)
  data = np.loadtxt(dirname+'/output.txt', delimiter=',', skiprows=1)
  plot(data[:,0], data[:,1]/volume, label="pSpatiocyte, $k$ = %s" %k,
      color=tableau20[colors[f]], linewidth=2.5)

num_A = 64000.
duration = 3.01

from scipy.integrate import odeint


# A -> B + C (kr)
#
# A = - kr*A
# B = kr*A
# C = kr*A
def f(x, t0, kr, keff):
  return np.array([
    -kr*x[0],
    kr*x[0],
    kr*x[0]
    ])

for i, k in enumerate(ks):
  init_state = np.array([num_A, 0, 0]) / volume
  ode_time = np.linspace(0,duration,500)
  ode_result = odeint(f, y0=init_state, t=ode_time, args=(k, 0))
  if (i == 0):
    plot(ode_time, ode_result[:,0], "--", color="k", label="ODE")
  else:
    plot(ode_time, ode_result[:,0], "--", color="k")

ax = gca()
handles, labels = ax.get_legend_handles_labels()
leg = legend(loc=(0.03,0.4), labelspacing=0.3, handlelength=2.5,
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
xlim(1e-2,3e+0)
plt.xscale('log')
tight_layout()
savefig('irreversible_reaction_A.pdf', format='pdf', dpi=600)
show()

