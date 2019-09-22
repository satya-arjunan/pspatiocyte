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

colors = [1, 0, 4, 2, 3]
skips = [1, 0, 0, 1, 0]
sims = ['pspatiocyte', 'spatiocyte', 'smoldyn', 'readdy', 'readdy_serial']
labels = ['pSpatiocyte', 'Spatiocyte', 'Smoldyn', 'Parallel ReaDDy', 'Serial ReaDDy']
data = []
height = 1.0  # the height of the bars

r = [[5,8,11,14],
     [4],
     [3],
     [2,7,10,13],
     [1]]

labelFontSize = 14
legendFontSize = 14
lineFontSize = 12

def autolabel(sim, rects, xpos='center'):
  ha = {'center': 'center', 'right': 'left', 'left': 'right'}
  offset = {'center': 0, 'right': 1, 'left': -1} 
  for i, rect in enumerate(rects):
    width = rect.get_width()
    xytext=(-16,-1)
    if (sim == 'pspatiocyte'):
      xytext=(10, -1)
    elif (sim == 'spatiocyte'): 
      xytext=(14, -1)
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

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_xlabel('Run time, $T$ (s)', size=labelFontSize)
ax.set_ylabel('CPU cores (threads)', size=labelFontSize)
ax.tick_params(axis='both', which='major', labelsize=lineFontSize)
ax.tick_params(axis='both', which='minor', labelsize=lineFontSize)
ax.set_yticks([3, 7.5, 10.5, 13.5], )
ax.set_yticklabels(('1', '2', '4', '8'))
ax.set_xlim(0,2200)
fig.tight_layout()
plt.savefig('benchmark.pdf', format='pdf', dpi=600)#, bbox_inches='tight')
plt.show()
