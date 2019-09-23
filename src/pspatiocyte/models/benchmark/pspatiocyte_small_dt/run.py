#!/usr/bin/env python3

import subprocess
import csv

threads = [1,2,4,8]
iterations = 3
times = []

for i in threads:
  time = 0.
  for j in range(1, iterations+1):
    print("running pSpatiocyte with", i, "thread(s), iteration", j, "of",
        iterations)
    result = subprocess.run(['mpirun', '-np', str(i), 'pspatiocyte_small_dt'],
        stdout=subprocess.PIPE)
    time = time + float(result.stdout.decode('utf-8').split('\n')[-2])
    print("elapsed time:", time/j)
  times.append(time/iterations)

subprocess.run(['python3', '../../../scripts/gather_timecourse.py'])

rows = zip(threads, times)
with open("elapsed_time.txt", "w+") as f:
  writer = csv.writer(f)
  writer.writerow(["threads", "elapsed_time"])
  for row in rows:
    writer.writerow(row)

