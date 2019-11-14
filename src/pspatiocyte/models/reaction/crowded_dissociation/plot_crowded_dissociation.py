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


tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
    (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
    (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
    (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
    (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)  


fig,ax=plt.subplots()

searches = [0, 1]
fractions = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
lines = ['--', '-']
colors = [2, 8, 10, 0, 6, 4]
iterations = 100

volume = 909.
for s, search in enumerate(searches):
  for f, fraction in enumerate(fractions): 
    data = []
    for i in range(iterations):
      if (search == 0):
        label = "original, $\phi="
      else:
        label = "vacated voxel, $\phi="
      fraction_str = '{:.2f}'.format(fraction)
      dirname = 'vacated_'+str(search)+'__fraction_'+fraction_str+'_'+str(i)
      print(dirname)
      if (i == 0):
        data = np.loadtxt(dirname+'/output.txt', delimiter=',', skiprows=1)
      else:
        data = data + np.loadtxt(dirname+'/output.txt', delimiter=',',
            skiprows=1)
    ax.plot(data[:,0]/iterations, data[:,1]/iterations/volume,
        label=label+str(fraction)+"$", color=tableau20[colors[f]], ls=lines[s],
        linewidth=2.5)

n_particles_s = 50000
duration = 12.0

from scipy.integrate import odeint

def f(x, t0, kcat, kr):
    return np.array([
        -kcat*x[0],
        kcat*x[0]
    ])

init_state = np.array([n_particles_s, 0.]) / volume
ode_time = np.logspace(-3,2,1000)
ode_result = odeint(f, y0=init_state, t=ode_time, args=(1.0, 0))

ax.plot(ode_time, ode_result[:,0], "--", color="k", label="ODE", linewidth=2.5)

labelFontSize = 16
legendFontSize = 15
lineFontSize = 15

handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(frameon=False, loc=(0.07,0.01), handlelength=2.5, handletextpad=0.8)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   

plt.xticks(size=labelFontSize)
plt.yticks(size=labelFontSize)

ax.tick_params(axis='both', which='major', direction='in', length=8, width=1,
    labelsize=lineFontSize)
ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1,
    labelsize=lineFontSize)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
plt.xlabel('Time, $t$ (s)',size=labelFontSize)
plt.ylabel("Concentration (\#\si{\micro}m$^{-3}$)",size=labelFontSize)
plt.xscale('log')
plt.yscale('log')
ax.set_ylim(1e-3,1e+2)
ax.set_xlim(0.01,11)
fig.tight_layout()
plt.savefig('dissociation_reaction.pdf', format='pdf', dpi=600)#, bbox_inches='tight')
plt.show()

