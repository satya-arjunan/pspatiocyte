import subprocess
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

matplotlib.rcParams["mathtext.fontset"] = "stix"
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
    r'\usepackage[helvet]{sfmath}', r'\usepackage[utf8]{inputenc}',
    r'\usepackage{arev}', r'\usepackage{siunitx}',
    r'\sisetup{math-micro={\usefont{T1}{phv}{m}{n}\text{µ}}}']
ratios = np.append(np.logspace(-1.5,1.5,12), [1])

Ds = [4.0, 0.06]
#Ds = [0.06]
markers = ['o', 'o']
labels = ['pSpatiocyte ($D=4$ \si{\micro}m$^2$s$^{-1}$)', 'pSpatiocyte ($D=0.06$ \si{\micro}m$^2$s$^{-1}$)', 'ODE (distributive)', 'ODE (processive)']

fig,ax=plt.subplots()

#D=4
data = np.loadtxt('ode/ode_distributive.csv', delimiter=',')
ax.semilogx(data[:,0], data[:,1], '-', label=labels[2], lw=2)

#D=0.06
data = np.loadtxt('ode/ode_processive.csv', delimiter=',')
ax.semilogx(data[:,0], data[:,1], '-', label=labels[3], lw=2)


for i, D in enumerate(Ds):
  x = []
  y = []
  for j, ratio in enumerate(ratios):
    ratio_str = '{:.4f}'.format(ratio)
    dirname = 'D_'+str(D)+'__ratio_'+ratio_str
    subprocess.run(['python3', '../../scripts/gather_timecourse.py', dirname])
    data = np.loadtxt(dirname+'/output.txt', delimiter=',', skiprows=1)
    NKT = data[:1,4][0]
    NKK = data[:1,1][0]
    NPP = data[:1,3][0]
    #print("intial NKT:",NKT,"NKK:",NKK,"NPP:",NPP)
    KK_over_PP = (NKK/NPP) 
    x.append(KK_over_PP)
    rows = len(data)
    Kpp = np.mean(data[-int(rows*0.5):,2])
    y.append(Kpp/NKT)
  ax.semilogx(x,y, markers[i], label=labels[i], markersize=7)


labelFontSize = 16
legendFontSize = 15
lineFontSize = 15

handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles[::-1], labels[::-1], loc=(0.03,0.67), frameon=False)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   

ax.set_xlabel('KK$_0$/PP$_0$', size=labelFontSize)
ax.set_ylabel('Kpp/K$_0$', size=labelFontSize)
plt.xticks(size=labelFontSize)
plt.yticks(size=labelFontSize)

ax.tick_params(axis='both', which='major', direction='in', length=8, width=1,
    labelsize=lineFontSize)
ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1,
    labelsize=lineFontSize)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
fig.tight_layout()
plt.savefig('mapk_ratio.pdf', format='pdf', dpi=600)#, bbox_inches='tight')
plt.show()
