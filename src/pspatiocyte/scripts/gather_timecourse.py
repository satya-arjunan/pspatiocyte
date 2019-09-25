import os
import sys
import numpy as np

def main(): 
  dirname = 'output'
  outputdir = ''
  if (len(sys.argv) >= 2):
    dirname = sys.argv[1]
    outputdir = dirname + os.sep
  print("Gathering timecourse data from directory:",dirname)
  inputf = dirname + os.sep + "timecourses_%06d"%(0) 
  if (os.path.exists(inputf) == False):
    return 0
  f = open(inputf, 'r')
  header = f.readline().strip()
  nproc = int(header.split(',')[0].split(' ')[1])
  result = np.loadtxt(inputf, delimiter=',', skiprows=1)
  ofname = outputdir + "output.txt"
  # repeat
  for i in range(1, nproc):
    inputf  = dirname + os.sep + "timecourses_%06d"%(i)
    next_result = np.loadtxt(inputf, delimiter=',', skiprows=1)
    result[:,1:] = result[:,1:] + next_result[:,1:] 
  # Total timecourse
  np.savetxt(ofname, result, fmt=''.join(['%.10f']+[',%i']*(result.shape[1]-1)), header=header, delimiter=',')

if __name__ == "__main__":
  main()
	

