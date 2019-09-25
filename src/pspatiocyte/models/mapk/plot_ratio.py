import subprocess
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

matplotlib.rcParams["mathtext.fontset"] = "stix"
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
    r'\usepackage[helvet]{sfmath}']

ratios = np.logspace(-1.5,1.5,12)
Ds = [0.06, 4.0]
#Ds = [0.06]
markers = ['o', 'o']
labels = ['pSpatiocyte (D=0.06)', 'pSpatiocyte (D=4.0)', 'ODE (processive)', 'ODE (distributive)']

fig,ax=plt.subplots()

#D=0.06
data = np.loadtxt('ode/ode_processive.csv', delimiter=',')
ax.semilogx(data[:,0], data[:,1], '-', label=labels[2])

#D=4
data = np.loadtxt('ode/ode_distributive.csv', delimiter=',')
ax.semilogx(data[:,0], data[:,1], '-', label=labels[3])

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
    # KK/PP
    KK_over_PP = (NKK/NPP) 
    x.append(KK_over_PP)
    rows = len(data)
    Kpp = np.mean(data[-int(rows*0.3):,2])
    #Kp_PP = np.mean(data[-int(rows*0.3):,9]) # Kpp/NKT
    Kp_PP = 0
    y.append((Kpp+Kp_PP)/NKT)
  ax.semilogx(x,y, markers[i], label=labels[i], markersize=7)


labelFontSize = 15
legendFontSize = 13
lineFontSize = 15

handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(handles[::-1], labels[::-1], frameon=False)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   

ax.set_xlabel('KK$_0$/PP$_0$', size=labelFontSize)
ax.set_ylabel('Kpp/K$_0$', size=labelFontSize)
ax.tick_params(axis='both', which='major', labelsize=lineFontSize)
ax.tick_params(axis='both', which='minor', labelsize=lineFontSize)
fig.tight_layout()
plt.show()
