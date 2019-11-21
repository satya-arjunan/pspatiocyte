#!/usr/bin/env python3

import subprocess
import csv
import numpy as np

ratio = 1.0
Ds = [0.06, 4.0]
threads = [1,2,4,8]

data = []

for D in Ds:
  for i in threads:
    time = 0.
    print("running pSpatiocyte MAPK model with", i,
        "threads, diffusion coefficient, D:", D, "ratio:", ratio)
    ratio_str = '{:.4f}'.format(ratio)
    dirname = 'D_'+str(D)+'__ratio_'+ratio_str
    result = subprocess.run(['mpirun', '-np', str(i), 'mapk', dirname, str(D),
      str(ratio)], stdout=subprocess.PIPE)
    time = float(result.stdout.decode('utf-8').split('\n')[-2])
    data.append([i, D, time, ratio])
    print('\telapsed time:',time)

with open("elapsed_time_scalability.txt", "w+") as f:
  writer = csv.writer(f)
  writer.writerow(['threads', 'D', 'elapsed_time', 'ratio'])
  for row in data:
    writer.writerow(row)


