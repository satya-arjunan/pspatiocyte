import os
import sys
import math
import pandas as pd
import numpy as np

# constants
sqr31 = math.sqrt(3.0)
sqr13 = math.sqrt(1.0/3.0)
sqr83 = math.sqrt(8.0/3.0)

class RealCoord:
  def __init__(self, nproc, radius, size, decomp, psize, nspecies, nlines):
    self.nproc  = nproc
    self.radius = radius
    self.size   = size
    self.decomp = decomp
    self.psize  = psize
    self.nspecies = nspecies
    self.nlines = nlines

  def get_subvolume_origin(self, rank):
    zmod = int(rank%self.decomp[0])
    ymod = int(rank/self.decomp[0])%self.decomp[1]
    xmod = int(rank/(self.decomp[0]*self.decomp[1]))
    return [xmod*self.psize[0], ymod*self.psize[1], zmod*self.psize[2]] 

  def get_real_coordinate(self, indexlist): 
    xind = indexlist[0]
    yind = indexlist[1]
    zind = indexlist[2] 
    rcx = self.radius * ( xind*2 - ((yind)%2) )
    rcy = self.radius * ( yind*sqr31 + 2*(zind%2)*sqr13 )
    rcz = self.radius * ( zind*sqr83 ) 
    rcoord = (rcx,rcy,rcz)
    return rcoord

  def index_to_coordinate(self, index, origin):
    pxg = self.psize[0]+2
    pyg = self.psize[1]+2
    pzg = self.psize[2]+2
    zind  = int(index%pzg)
    yind  = int((index%(pyg*pzg))/pzg)
    xind  = int(index/(pyg*pzg))
    xind += origin[0]
    yind += origin[1]
    zind += origin[2]
    return self.get_real_coordinate((xind,yind,zind))
	
  def write(self, outputdir):
    outfile = open(outputdir+"coordinates.spa", 'w')
    files = [open(outputdir+"coordinates_%06d"%(x), 'r') 
        for x in range(self.nproc)]
    for afile in files:
      header = afile.readline().split(',')[:-1]
    header.append("voxel_radius=%e"%(self.radius))
    outfile.write(','.join(header)+os.linesep) #write header
    origins = [self.get_subvolume_origin(x) for x in range(self.nproc)]
    for i in range(self.nlines):
      outline = [[]]*(1+self.nspecies)
      for rank in range(self.nproc):
        line = files[rank].readline().strip().split(',')
        outline[0] = line[0]
        for species in range(len(line[1:])):
          voxel_index = line[1:][species]
          if (voxel_index):
            voxel_index = map(int,voxel_index.split(' '))
            coordinates = []
            for index in voxel_index:
              coordinates.append(self.index_to_coordinate(index, origins[rank]))
              coordinates.append(next(voxel_index))
            if (outline[species+1]):
              outline[species+1].extend(coordinates)
            else:
              outline[species+1] = coordinates
      outlinestr = map(str,outline) 
      outline = [
        x.strip('[').strip(']').replace(',', '').replace('(','').replace(')','')
        for x in outlinestr]
      outfile.write(','.join(outline)+os.linesep)

def main():
  dirname = 'output' + os.sep
  outputdir = dirname
  if (len(sys.argv) >= 2):
    dirname = sys.argv[1]
    outputdir = dirname + os.sep
  inputf = outputdir+"coordinates_%06d"%(0) 
  f = open(inputf, 'r')
  header = f.readline().strip().split(',')
  nlines = 0
  for line in f:
    nlines += 1
  radius, nproc, nx, ny, nz = map(float, header[-1].split(' '))
  nproc = int(nproc)
  nspecies = len(header)-2
  size = map(int, (nx, ny, nz))
  #decomposition
  nd = math.log(nproc,2)
  dx = 2**(int(nd/3))
  dy = 2**(int((nd-(int(nd/3)))/2))
  dz = int(nproc/(dy*dx))
  decomp = (dx,dy,dz)
  # process size
  psize = (int(nx/dx), int(ny/dy), int(nz/dz))
  obj = RealCoord(nproc, radius, size, decomp, psize, nspecies, nlines)
  obj.write(outputdir)

if __name__ == "__main__":
  main()




