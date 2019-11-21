import subprocess
import numpy as np
import matplotlib
import os
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd

matplotlib.rcParams["mathtext.fontset"] = "stix"
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
    r'\usepackage[helvet]{sfmath}', r'\usepackage[utf8]{inputenc}',
    r'\usepackage{arev}', r'\usepackage{siunitx}',
    r'\sisetup{math-micro={\usefont{T1}{phv}{m}{n}\text{µ}}}']

labelFontSize = 16
legendFontSize = 15
lineFontSize = 15
matplotlib.rcParams.update({'font.size': labelFontSize})

fig,ax=plt.subplots()

fractions = [0.0, 0.3, 0.5]
D = 10e-12
times = []

for fraction in fractions:
  filename = 'fraction_%.1f.csv' %fraction
  df = pd.read_csv(filename)
  times = df.iloc[:,0]
  y = df.iloc[:,1]
  ax.plot(times/1e-6, y/1e-12, 
      label=r'pSpatiocyte ($\phi=%0.1f$, $D_0=%d$ \si{\micro}m$^{2}$s$^{-1}$)'
      %(fraction, int(D/1e-12)), lw=2)
  
ax.plot(times/1e-6, 6.0*D*times/1e-12, color='darkblue', ls='--', label='6$Dt$, $D=10$ \si{\micro}m$^{2}$s$^{-1}$')
ax.plot(times/1e-6, 6.0*6.3e-12*times/1e-12, color='sienna', ls='--', label='6$Dt$, $D=6.3$ \si{\micro}m$^{2}$s$^{-1}$')
ax.plot(times/1e-6, 6.0*3.65e-12*times/1e-12, color='darkgreen', ls='--', label='6$Dt$, $D=3.6$ \si{\micro}m$^{2}$s$^{-1}$')

handles, labels = ax.get_legend_handles_labels()
leg = ax.legend(frameon=False, loc=2, handletextpad=0.3, labelspacing=0.2, handlelength=1.5)
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
ax.set_xlim(0,160)
ax.set_ylim(0,0.012)
ax.set_xlabel('Time, $t$ (\si{\micro}s)', size=labelFontSize)
ax.set_ylabel('Mean-squared displacement (\si{\micro}m$^{2}$)', size=labelFontSize)
fig.tight_layout()
plt.savefig('crowded_diffusion_csv.pdf', format='pdf', dpi=600)
plt.show()





