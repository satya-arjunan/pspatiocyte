#!/usr/bin/env python3

import subprocess
import csv

#threads = [1,2,4,8]
threads = [8]
iterations = 1
times = []

for i in threads:
  subprocess.run(['rm', '-rf', 'output'])
  subprocess.run(['mkdir', 'output'])
  time = 0.
  for j in range(1, iterations+1):
    print("running pSpatiocyte with", i, "thread(s), iteration", j, "of",
        iterations)
    result = subprocess.run(['mpirun', '-np', str(i), 'benchmark'],
        stdout=subprocess.PIPE)
    time = time + float(result.stdout.decode('utf-8').split('\n')[-2])
    print("elapsed time:", time/j)
  times.append(time/iterations)

subprocess.run(['python3', '../../scripts/gather_timecourse.py', '8'])

rows = zip(threads, times)
with open("elapsed_time.txt", "w+") as f:
  writer = csv.writer(f)
  writer.writerow(["threads", "elapsed_time"])
  for row in rows:
    writer.writerow(row)

