import matplotlib
import matplotlib.pyplot as plt
import numpy as np
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

colors = [0, 2, 8, 4, 1, 3]
sims = ['elapsed_time_scalability.txt']
labels = ['pSpatiocyte ($D=0.06$ \si{\micro}m$^{2}$s$^{-1}$)', 'pSpatiocyte ($D=4$ \si{\micro}m$^{2}$s$^{-1}$)']
data = []
height = 1.0  # the height of the bars

r = [[2,5,8,11],
     [1,4,7,10]]

labelFontSize = 15
legendFontSize = 15
lineFontSize = 15

def autolabel(rects, xpos='center'):
  ha = {'center': 'center', 'right': 'left', 'left': 'right'}
  offset = {'center': 0, 'right': 1, 'left': -1} 
  for i, rect in enumerate(rects):
    width = rect.get_width()
    if (width >= 1000):
      xytext=(-19,-3)
    elif(width >= 100):
      xytext=(-15,-3)
    ax.annotate('{}'.format(width),
        xy=(width, rect.get_y() + rect.get_height()/ 4),
        xytext=xytext,  # use 3 points offset
        textcoords="offset points",  # in both directions
        ha=ha[xpos], va='bottom', size=lineFontSize)

fig, ax = plt.subplots()

for i in range(len(sims)): 
  data = np.around(np.genfromtxt(sims[i], delimiter=',',
    skip_header=1).T[2]).astype(int)
  first = data[:-4]
  second = data[4:]
  rects = ax.barh(r[0], first, height, label=labels[0],
      color=tableau20[colors[0]])
  autolabel(rects)
  rects = ax.barh(r[1], second, height, label=labels[1],
      color=tableau20[colors[1]])
  autolabel(rects)

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
ax.set_xlabel('Runtime, $T$ (s)', size=16)
ax.set_ylabel('CPU cores (processes)', size=16)
ax.tick_params(axis='both', which='major', labelsize=labelFontSize)
ax.tick_params(axis='both', which='minor', labelsize=labelFontSize)
ax.set_yticks([1.5, 4.5, 7.5, 10.5], )
ax.set_yticklabels(('1', '2', '4', '8'))
#ax.set_xlim(0,2200)
fig.tight_layout()
plt.savefig('mapk_scalability.pdf', format='pdf', dpi=600)#, bbox_inches='tight')
plt.show()
