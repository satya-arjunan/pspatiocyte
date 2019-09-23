import subprocess
import numpy as np
import matplotlib.pyplot as plt

ratios = np.logspace(-1.5,1.5,6)
Ds = [0.06, 4.0]
markers = ['o', 'x']
fig,ax=plt.subplots(1,1,figsize=(12,8))

for i, D in enumerate(Ds):
  x = []
  y = []
  for ratio in ratios:
    ratio_str = '{:.4f}'.format(ratio)
    dirname = 'D_'+str(D)+'__ratio_'+ratio_str
    subprocess.run(['python3', '../../scripts/gather_timecourse.py', dirname])
    data = np.loadtxt(dirname+'/output.txt', delimiter=',', skiprows=1)
    NKT = data[:1,4][0]
    # KK/PP
    KK_over_PP = (data[:1,1]/data[:1,3])[0] 
    x.append(KK_over_PP)
    rows = len(data)
    print(rows)
    Kpp = np.mean(data[-int(rows*0.3):,2])
    # Kpp/NKT
    y.append(Kpp/NKT)
  x = np.sort(x)
  y = np.sort(y)
  ax.semilogx(x,y, markers[i], markersize=10)

plt.show()
