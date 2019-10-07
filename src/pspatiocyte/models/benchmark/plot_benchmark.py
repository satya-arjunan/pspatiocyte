import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc


tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
    (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
    (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
    (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
    (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)  

colors = [0, 1, 8, 4, 2, 3]
skips = [1, 1, 0, 0, 1, 0]
sims = ['pspatiocyte', 'pspatiocyte_small_dt', 'spatiocyte', 'smoldyn', 'readdy', 'readdy_serial']
labels = ['pSpatiocyte ($\Delta t=0.5$ ms)', 'pSpatiocyte ($\Delta t=0.2$ ms)', 'Spatiocyte ($\Delta t=0.5$ ms)', 'Smoldyn ($\Delta t=1$ ms)', 'Parallel ReaDDy ($\Delta t=1$ ms)', 'Serial ReaDDy ($\Delta t=1$ ms)']
data = []
height = 1.0  # the height of the bars

r = [[6,10,14,18],
     [5,9,13,17],
     [4],
     [3],
     [2,8,12,16],
     [1]]

labelFontSize = 15
legendFontSize = 15
lineFontSize = 12

def autolabel(sim, rects, xpos='center'):
  ha = {'center': 'center', 'right': 'left', 'left': 'right'}
  offset = {'center': 0, 'right': 1, 'left': -1} 
  for i, rect in enumerate(rects):
    width = rect.get_width()
    if (width >= 1000):
      xytext=(-19,-3)
    elif(width >= 100):
      xytext=(-15,-3)
    if (sim == 'pspatiocyte' or sim == 'pspatiocyte_small_dt'):
      if(width >= 100):
        xytext=(14, -3)
      else:
        xytext=(12, -3)
    elif (sim == 'spatiocyte'): 
      xytext=(14, -3)
    ax.annotate('{}'.format(width),
        xy=(width, rect.get_y() + rect.get_height()/ 4),
        xytext=xytext,  # use 3 points offset
        textcoords="offset points",  # in both directions
        ha=ha[xpos], va='bottom', size=lineFontSize)

fig, ax = plt.subplots()

for i in range(len(sims)): 
  if (skips[i]):
    data = np.around(np.genfromtxt(sims[i]+'/elapsed_time.txt',
      delimiter=',', skip_header=1).T[1]).astype(int)
  else:
    data = np.around(np.genfromtxt(sims[i]+'/elapsed_time.txt').T).astype(int) 
  rects = ax.barh(r[i], data, height, label=labels[i], color=tableau20[colors[i]])
  autolabel(sims[i], rects)

leg = ax.legend(frameon=False)
for t in leg.get_texts():
  t.set_fontsize(legendFontSize)   


plt.xticks(size=labelFontSize)
plt.yticks(size=labelFontSize)

ax.tick_params(axis='x', which='major', direction='in', length=8, width=1,
    labelsize=labelFontSize)
ax.tick_params(axis='y', which='major', length=8, width=1, labelsize=labelFontSize)

ax.xaxis.set_ticks_position('both')
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_xlabel('Run time, $T$ (s)', size=16)
ax.set_ylabel('CPU cores (threads)', size=16)
ax.tick_params(axis='both', which='major', labelsize=labelFontSize)
ax.tick_params(axis='both', which='minor', labelsize=labelFontSize)
ax.set_yticks([3, 9, 13, 17], )
ax.set_yticklabels(('1', '2', '4', '8'))
ax.set_xlim(0,2200)
fig.tight_layout()
plt.savefig('benchmark.pdf', format='pdf', dpi=600)#, bbox_inches='tight')
plt.show()
