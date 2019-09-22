#!/usr/bin/env python3

import subprocess

iterations = 3
time = 0.
for i in range(1, iterations+1):
  print("running Smoldyn, iteration", i, "of", iterations)
  result = subprocess.run(['smoldyn', 'model.txt', '-w'],
      stdout=subprocess.PIPE)
  time = time + float(result.stdout.decode('utf-8').split('\n')[-2].
      split(' ')[-2])
time = time/iterations
print("elapsed time (s):", time)

with open("elapsed_time.txt", "w+") as f:
  f.write("%f" % time)
