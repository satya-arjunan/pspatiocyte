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
Ds = [0.06, 4.0]
#Ds = [0.06]
markers = ['o', 'o']
labels = ['pSpatiocyte (D = 0.06 \si{\micro}m$^2$s$^{-1}$)', 'pSpatiocyte (D = 4 \si{\micro}m$^2$s$^{-1}$)', 'ODE (processive)', 'ODE (distributive)']

fig,ax=plt.subplots()

searches = [0]
#fractions = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
fractions = [0.5, 0.6, 0.7]
iterations = 10

volume = 909.
for search in searches:
  for fraction in fractions: 
    data = []
    for i in range(iterations):
      fraction_str = '{:.2f}'.format(fraction)
      dirname = 'forced_'+str(search)+'__fraction_'+fraction_str+'_'+str(i)
      subprocess.run(['python3', '../../scripts/gather_timecourse.py', dirname])
      if (i == 0):
        data = np.loadtxt(dirname+'/output.txt', delimiter=',', skiprows=1)
      else:
        data = data + np.loadtxt(dirname+'/output.txt', delimiter=',',
            skiprows=1)
    ax.plot(data[:,0]/iterations, data[:,1]/iterations/volume,
        label=dirname.replace('_','\_'))

n_particles_s = 50000
duration = 10.0

from scipy.integrate import odeint

def f(x, t0, kcat, kr):
    """
    x: state vector with concentrations of S, P
    """
    return np.array([
        -kcat*x[0],
        kcat*x[0]
    ])

init_state = np.array([n_particles_s, 0.]) / volume
ode_time = np.linspace(0.,duration,100000)
ode_result = odeint(f, y0=init_state, t=ode_time, args=(1.0, 0))

ax.plot(ode_time, ode_result[:,0], "--", color="k", alpha=.5, label="ODE")

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
ax.set_xlim(0.1,10.0)
plt.xscale('log')
plt.yscale('log')
fig.tight_layout()
plt.savefig('dissociation.pdf', format='pdf', dpi=600)#, bbox_inches='tight')
plt.show()

