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

matplotlib.rcParams["mathtext.fontset"] = "stix"
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
    r'\usepackage[helvet]{sfmath}']

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

colors = [1, 0, 4, 2, 3, 5, 6, 7, 8]
cols = [1, 2, 3, 4, 5, 6, 7, 8]

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

#fileNames = ['D_0.06__ratio_1.3689/output.txt']
fileNames = ['D_0.06__ratio_2.5650/output.txt']
legendTitles = ['pSpatiocyte (8 cores) ($\Delta t=0.5\ \mathrm{ms},\ T=11\ \mathrm{s}$)','Spatiocyte ($\Delta t=0.5\ \mathrm{ms}, T=139\ \mathrm{s}$)','Smoldyn ($\Delta t=1\ \mathrm{ms},\ T=449\ \mathrm{s}$)','Parallel ReaDDy (8 cores) ($\Delta t=1\ \mathrm{ms}, T=549\ \mathrm{s}$)','Serial ReaDDy ($\Delta t=1\ \mathrm{ms}, T=2197\ \mathrm{s}$)']
speciesList = ['E','S','ES','P']
lines = ['-','-','-','-','-','-','-']
opacity = [1, 1, 1, 1, 1]

volume = 10.0144
for f in reversed(range(len(fileNames))):
  if (os.path.isfile(path+fileNames[f])):
    deli = ','
    data = genfromtxt(path+fileNames[f], delimiter=deli, skip_header=1).T
    colSize = len(data)-1
    for i in cols:
      if (i == 2):
        plot(data[0], (data[i]+data[9])/volume, ls=lines[i-1],
            color=tableau20[colors[i-1]], label=legendTitles[i-1], linewidth=3)
      else:
        plot(data[0], data[i]/volume, ls=lines[0], color=tableau20[colors[i-1]],
            linewidth=3)

molecule_radius = 0.0025
ka1, kd1, kcat1 = 0.04483455086786913, 1.35, 1.5
ka2, kd2, kcat2 = 0.09299017957780264, 1.73, 15.0
trel = 1e-6
k7 = math.log(2.)/trel
D = 0.06
duration = 300.
ratio = 1.3689

def kon(k):
  kD = 4*3.14*2*molecule_radius*2*D
  return k*kD/(k+kD)

def koff(kd,ka):
  return kon(ka)*kd/ka 

kon1 = kon(ka1)
koff1 = koff(kd1, ka1)

# (K + KK == K_KK | (kon1, koff1)
#     > Kp_KK | kcat1
#     > Kpp + KK | kcat2)
# (Kpp + PP == Kpp_PP | (kon1, koff1)
#     > Kp + PP | kcat1
#     > K + PP | kcat2)

# K = -kon1*KK*K + koff1*K_KK + kcat2*Kp*PP
# KK = -kon1*KK*K + koff1*K_KK + kcat2*Kp_KK
# K_KK = kon1*KK*K - koff1*K_KK - kcat1*K_KK
# Kp_KK = kcat1*K_KK - kcat2*Kp_KK 
# Kpp = kcat2*Kp_KK - kon1*Kpp*PP + koff1*Kpp_PP
# PP = -kon1*Kpp*PP + koff1*Kpp_PP + kcat1*Kpp_PP
# Kpp_PP = kon1*Kpp*PP - koff1*Kpp_PP - kcat1*Kpp_PP
# Kp = kcat1*Kpp_PP - kcat2*Kp*PP

#Kp_PP

#Kp_PP,KKa,PPa

# x[0] = -kon1*x[1]*x[0] + koff1*x[2] + kcat2*x[7]*x[5]
# x[1] = -kon1*x[1]*x[0] + koff1*x[2] + kcat2*x[3]
# x[2] = kon1*x[1]*x[0] - koff1*x[2] - kcat1*x[2]
# x[3] = kcat1*x[2] - kcat2*x[3] 
# x[4] = kcat2*x[3] - kon1*x[4]*x[5] + koff1*x[6]
# x[5] = -kon1*x[4]*x[5] + koff1*x[6] + kcat1*x[6]
# x[6] = kon1*x[4]*x[5] - koff1*x[6] - kcat1*x[6]
# x[7] = kcat1*x[6] - kcat2*x[7]*x[5]

from scipy.integrate import odeint

def f(x, t0):
  return np.array([ 
    -kon1*x[1]*x[0] + koff1*x[2] + kcat2*x[7]*x[5],
    -kon1*x[1]*x[0] + koff1*x[2] + kcat2*x[3],
    kon1*x[1]*x[0] - koff1*x[2] - kcat1*x[2],
    kcat1*x[2] - kcat2*x[3],
    kcat2*x[3] - kon1*x[4]*x[5] + koff1*x[6],
    -kon1*x[4]*x[5] + koff1*x[6] + kcat1*x[6],
    kon1*x[4]*x[5] - koff1*x[6] - kcat1*x[6],
    kcat1*x[6] - kcat2*x[7]*x[5]
    ])

NKT = int(120*volume) # total K
ratio = np.logspace(-1.5,1.5,12)[7]
print(ratio)
ode_time = np.linspace(0.,duration,100000)
result = []
NPP = np.rint(60*volume/(ratio+1)) # initial PP
NKK = int(60*volume-NPP) # initial KK 
init_state = np.array([NKT, NKK, 0., 0., 0., NPP, 0., 0.]) / volume
ode_result = odeint(f, y0=init_state, t=ode_time)
plot(ode_time, ode_result[:,0], "--", color="k", alpha=.5, label="Mass Action")
plot(ode_time, ode_result[:,1], "--", color="k", alpha=.5)
plot(ode_time, ode_result[:,2], "--", color="k", alpha=.5)
plot(ode_time, ode_result[:,3], "--", color="k", alpha=.5)
plot(ode_time, ode_result[:,4], "--", color="k", alpha=.5)
plot(ode_time, ode_result[:,5], "--", color="k", alpha=.5)
plot(ode_time, ode_result[:,6], "--", color="k", alpha=.5)
plot(ode_time, ode_result[:,7], "--", color="k", alpha=.5)

ax = gca()
handles, labels = ax.get_legend_handles_labels()
leg = legend(handles[::-1], labels[::-1], labelspacing=0.12, handlelength=1.0, handletextpad=0.3, frameon=False)
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

xlabel('Time, $t$',size=labelFontSize)
ylabel("Concentration")
xlim(0,300)
tight_layout(pad=0)
savefig('reaction.pdf', format='pdf', dpi=600)#, bbox_inches='tight')
show()


