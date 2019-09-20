import os,sys
import numpy

def main(): 
  if len(sys.argv) < 2:
    print("Nproc, No is need")
    quit() 
  # Number of procs
  nproc = int(sys.argv[1])
  inputf = "output/timecourses_%06d"%(0) 
  f = open(inputf, 'r')
  header = f.readline()
  result = numpy.loadtxt(inputf, delimiter=',', skiprows=1)
  ofname = "output.txt"
  # repeat
  for i in range(1, nproc):
    input  = "./output/timecourses_%06d"%(i)
    temp   = numpy.loadtxt(input, delimiter=',', skiprows=1)
    result[:,1:] = result[:,1:] + temp[:,1:] 
  # Total timecourse
  numpy.savetxt(ofname, result, fmt=''.join(['%.10f']+[',%i']*(result.shape[1]-1)), header=header, delimiter=',')

if __name__ == "__main__":
  main()
	

