import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib
from matplotlib import rc

matplotlib.rcParams["mathtext.fontset"] = "stix"
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}',
    r'\usepackage[helvet]{sfmath}', r'\usepackage[utf8]{inputenc}',
    r'\usepackage{arev}', r'\usepackage{siunitx}',
    r'\sisetup{math-micro={\usefont{T1}{phv}{m}{n}\text{µ}}}']

labelFontSize = 10
legendFontSize = 13
lineFontSize = 15

matplotlib.rcParams.update({'font.size': labelFontSize})

fig = plt.figure()
ax = plt.axes(projection='3d')

inputf = 'coordinates.spa'
f = open(inputf, 'r')
header = f.readline().strip().split(',')
voxel_radius = float(header[-1].split('=')[1])
data = np.loadtxt(inputf, delimiter=',', skiprows=1, dtype=np.str)
tracks = {}
for row in data:
  for species_id, col in enumerate(row[1:]): #skip the first time col
    if col:
      coords = col.split(' ')
      for i in range(0, len(coords), 4):
        x = float(coords[i])
        y = float(coords[i+1])
        z = float(coords[i+2])
        mol_id = int(coords[i+3])
        if mol_id not in tracks:
          tracks[mol_id] = []
        tracks[mol_id].append([x, y, z])

for mol_id in tracks:
  coords = tracks[mol_id]
  track = np.asarray(coords)
  ax.plot3D(track[:,0], track[:,1], track[:,2])

ax.set_xlim(0, 2.5)
ax.set_ylim(0, 2.5)
ax.set_zlim(0, 2.0)
ax.set_xlabel('x (\si{\micro}m)', labelpad=10)
ax.set_ylabel('y (\si{\micro}m)', labelpad=8)
ax.set_zlabel('z (\si{\micro}m)', labelpad=5)
plt.gca().invert_yaxis()
plt.savefig('trajectory.pdf', format='pdf', dpi=600)#, bbox_inches='tight')
plt.show()


