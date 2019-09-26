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

matplotlib.rcParams.update({'font.size': labelFontSize})
#matplotlib.rcParams['figure.figsize'] = 9.1,7

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

colors = [1, 0, 4, 2, 16, 10, 6, 18, 8]
cols = [1, 2, 3, 4, 5, 6, 7, 8, 9]

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
#fileNames = ['D_0.06__ratio_2.5650/output.txt']
#fileNames = ['D_4.0__ratio_0.3899/output.txt']
#fileNames = ['D_4.0__ratio_0.2081/output.txt']
fileNames = ['D_4.0__ratio_1.0000/output.txt']
lines = ['-','-','-','-','-','-','-']
opacity = [1, 1, 1, 1, 1]
colormaps = {}

volume = 10.0144
for f in range(len(fileNames)):
  if (os.path.isfile(path+fileNames[f])):
    deli = ','
    af = open(fileNames[f], 'r')
    legendTitles = af.readline().strip().split(",")
    data = genfromtxt(path+fileNames[f], delimiter=deli, skip_header=1).T
    colSize = len(data)-1
    for i in cols:
      label = legendTitles[i].replace('_','\_')
      colormaps[legendTitles[i]] = colors[i-1]
      plot(data[0], data[i]/volume, ls=lines[0], color=tableau20[colors[i-1]],
          label=label, linewidth=2.5)

molecule_radius = 0.0025
ka1, kd1, kcat1 = 0.04483455086786913, 1.35, 1.5
ka2, kd2, kcat2 = 0.09299017957780264, 1.73, 15.0
trel = 1e-6
k7 = math.log(2.)/trel
D = 4.
duration = 300.

def kon(k):
  kD = 4*3.14*2*molecule_radius*2*D
  return k*kD/(k+kD)

def koff(kd,ka):
  return kon(ka)*kd/ka 

kon1 = kon(ka1)
koff1 = koff(kd1, ka1)
kon2 = kon(ka2)
koff2 = koff(kd2, ka2)

# (KK + K == K_KK | (kon1, koff1)
#     > Kp + KKa | kcat1)
# (KK + Kp == Kp_KK | (kon2, koff2)
#     > Kpp + KKa | kcat2)
# (KKa > KK | k7)
# (Kpp + PP == Kpp_PP | (kon1, koff1)
#     > Kp + PPa | kcat1)     
# (Kp + PP == Kp_PP | (kon2, koff2)
#     > K + PPa | kcat2)
# (PPa > PP | k7)

# KK = -kon1*KK*K + koff1*K_KK - kon2*KK*Kp + koff2*Kp_KK + k7*KKa
# K =  -kon1*KK*K + koff1*K_KK + kcat2*Kp_PP
# K_KK = kon1*KK*K - koff1*K_KK - kcat1*K_KK
# Kp = kcat1*K_KK - kon2*KK*Kp + koff2*Kp_KK + kcat1*Kpp_PP - kon2*Kp*PP + koff2*Kp_PP
# KKa = kcat1*K_KK + kcat2*Kp_KK - k7*KKa
# Kp_KK = kon2*KK*Kp - koff2*Kp_KK - kcat2*Kp_KK
# Kpp = kcat2*Kp_KK - kon1*Kpp*PP + koff1*Kpp_PP
# PP = -kon1*Kpp*PP + koff1*Kpp_PP - kon2*Kp*PP + koff2*Kp_PP + k7*PPa
# Kpp_PP = kon1*Kpp*PP - koff1*Kpp_PP - kcat1*Kpp_PP
# PPa = kcat1*Kpp_PP + kcat2*Kp_PP - k7*PPa
# Kp_PP = kon2*Kp*PP - koff2*Kp_PP - kcat2*Kp_PP

# x[0] = -kon1*x[0]*x[1] + koff1*x[2] - kon2*x[0]*x[3] + koff2*x[5] + k7*x[4]
# x[1] =  -kon1*x[0]*x[1] + koff1*x[2] + kcat2*x[10]
# x[2] = kon1*x[0]*x[1] - koff1*x[2] - kcat1*x[2]
# x[3] = kcat1*x[2] - kon2*x[0]*x[3] + koff2*x[5] + kcat1*x[8] - kon2*x[3]*x[7] + koff2*x[10]
# x[4] = kcat1*x[2] + kcat2*x[5] - k7*x[4]
# x[5] = kon2*x[0]*x[3] - koff2*x[5] - kcat2*x[5]
# x[6] = kcat2*x[5] - kon1*x[6]*x[7] + koff1*x[8]
# x[7] = -kon1*x[6]*x[7] + koff1*x[8] - kon2*x[3]*x[7] + koff2*x[10] + k7*x[9]
# x[8] = kon1*x[6]*x[7] - koff1*x[8] - kcat1*x[8]
# x[9] = kcat1*x[8] + kcat2*x[10] - k7*x[9]
# x[10] = kon2*x[3]*x[7] - koff2*x[10] - kcat2*x[10]

from scipy.integrate import odeint

def f(x, t0):
  return np.array([ 
    -kon1*x[0]*x[1] + koff1*x[2] - kon2*x[0]*x[3] + koff2*x[5] + k7*x[4],
    -kon1*x[0]*x[1] + koff1*x[2] + kcat2*x[10],
    kon1*x[0]*x[1] - koff1*x[2] - kcat1*x[2],
    kcat1*x[2] - kon2*x[0]*x[3] + koff2*x[5] + kcat1*x[8] - kon2*x[3]*x[7] + koff2*x[10],
    kcat1*x[2] + kcat2*x[5] - k7*x[4],
    kon2*x[0]*x[3] - koff2*x[5] - kcat2*x[5],
    kcat2*x[5] - kon1*x[6]*x[7] + koff1*x[8],
    -kon1*x[6]*x[7] + koff1*x[8] - kon2*x[3]*x[7] + koff2*x[10] + k7*x[9],
    kon1*x[6]*x[7] - koff1*x[8] - kcat1*x[8],
    kcat1*x[8] + kcat2*x[10] - k7*x[9],
    kon2*x[3]*x[7] - koff2*x[10] - kcat2*x[10],
    ])

NKT = int(120*volume) # total K
#ratio = np.logspace(-1.5,1.5,12)[3]
ratio = 1.
ode_time = np.linspace(0.,duration,100000)
result = []
NPP = np.rint(60*volume/(ratio+1)) # initial PP
NKK = int(60*volume-NPP) # initial KK 
init_state = np.array([NKK, NKT, 0., 0., 0., 0., 0., NPP, 0., 0., 0.])/volume
ode_result = odeint(f, y0=init_state, t=ode_time)
plot(ode_time, ode_result[:,0], "--", color='k', label="ODE (distributive)", linewidth=1)
plot(ode_time, ode_result[:,1], "--", color='k', linewidth=1)
plot(ode_time, ode_result[:,2], "--", color='k', linewidth=1)
plot(ode_time, ode_result[:,3], "--", color='k', linewidth=1)
#plot(ode_time, ode_result[:,4], "--", label="KKa")
plot(ode_time, ode_result[:,5], "--", color='k', linewidth=1)
plot(ode_time, ode_result[:,6], "--", color='k', linewidth=1)
plot(ode_time, ode_result[:,7], "--", color='k', linewidth=1)
plot(ode_time, ode_result[:,8], "--", color='k', linewidth=1)
#plot(ode_time, ode_result[:,9], "--", label="PPa")
plot(ode_time, ode_result[:,10], "--", color='k', linewidth=1)

ax = gca()
handles, labels = ax.get_legend_handles_labels()
leg = legend(labelspacing=0.12, handlelength=1.2, handletextpad=0.5, frameon=False, ncol=2)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   

xticks(size=labelFontSize)
yticks(size=labelFontSize)

ax.tick_params(axis='both', which='major', labelsize=lineFontSize)
ax.tick_params(axis='both', which='minor', labelsize=lineFontSize)
xlabel('Time, $t$ (s)',size=labelFontSize)
ylabel("Concentration (\#\si{\micro}m$^{-3}$)",size=labelFontSize)
xlim(0,100)
tight_layout(pad=0)
savefig('output.pdf', format='pdf', dpi=600)#, bbox_inches='tight')
show()


