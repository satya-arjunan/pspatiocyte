#include "Species.hpp"
#include "Lattice.hpp"
#include "Compartment.hpp"

inline int max3( int a, int b, int c )
{
    const int t = ( a>b ? a : b );
    return ( t>c ? t : c );
}

// parallel constructor with ghost cells
Lattice::Lattice(string name, double r, ParallelEnvironment &pe,
                 const int invalid_id, const int vacant_id, const int ghost_id):
  invalid_id_(invalid_id),
  vacant_id_(vacant_id),
  ghost_id_(ghost_id),
  name_(name),
  nx_(pe.getispan()), 
  ny_(pe.getjspan()),
  nz_(pe.getkspan()),
  Nx_(nx_ + 2*GHOST_SIZE),
  Ny_(ny_ + 2*GHOST_SIZE),
  Nz_(nz_ + 2*GHOST_SIZE),
  radius_(r),
  SQR31_(sqrt(3.0)),
  SQR13_(sqrt(1.0 / 3.0)),
  SQR83_(sqrt(8.0 / 3.0)) {
    voxelVector_.resize(Nx_ * Ny_ * Nz_);
    // buffers for quater ghost
    const int Nxh = Nx_/2;
    const int Nyh = Ny_/2;
    const int Nzh = Nz_/2;
    // quarter ghost, see Lattice::loadGhost()
    //bufsize = 4*max3( (Nxh+2)*(Nyh+2), (Nyh+2)*(Nzh+2), (Nzh+2)*(Nxh+2) );
    // full ghost, see Lattice::loadConfiguration()
    bufsize = max3( Nx_*Ny_, Ny_*Nz_, Nz_*Nx_ );
    inbound = new int [bufsize]; 
    outboundx = new int [bufsize];
    outboundy = new int [bufsize];
    outboundz = new int [bufsize];

    // preallocate enough spaces to suppress system calls
    const int molesize = bufsize/4;
    jumpincoords.reserve( molesize );
    spillcoordsX.reserve( molesize );
    spillcoordsY.reserve( molesize );
    spillcoordsZ.reserve( molesize );

    for (int i=0; i<Nx_; ++i) {
      for (int j=0; j<Ny_; ++j) {
        for(int k=0; k<Nz_; ++k) {
          const int lc = linearCoord(i,j,k);
          if ( i<GHOST_SIZE || Nx_-GHOST_SIZE<=i ||
               j<GHOST_SIZE || Ny_-GHOST_SIZE<=j ||
               k<GHOST_SIZE || Nz_-GHOST_SIZE<=k ) {
            voxelVector_[lc].species_id = ghost_id_;
          }
          else {
            voxelVector_[lc].species_id = vacant_id_;
          }
        }
      }
    }
  }

/* 
 *  twelve adjacent voxels on twisted cartesian coordinate system
 *  
 *           top-plane            mid-plane              bottom-plane
 *                                                     
 *  if(k%2==0 && j%2==0) 
 *   j+1                             NW(3)  NE(4)
 *   j      TN(0)             W(8)    x      E(5)       BN(9)
 *   j-1    TW(2)  TE(1)             SW(7)  SE(6)       BW(11)   BE(10)
 *           i     i+1        i-1     i     i+1          i      i+1
 *                                                      
 *  if(k%2==0 && j%2==1) 
 *   j+1                      NW(3)  NE(4)
 *   j             TN(0)       W(8)   x     E(5)                BN(9)
 *   j-1    TW(2)  TE(1)      SW(7)  SE(6)              BW(11)  BE(10)
 *          i-1     i         i-1     i     i+1         i-1      i
 *
 *  if(k%2==1 && j%2==0) 
 *   j+1    TW(0)  TE(1)             NW(3)  NE(4)       BW(9)   BE(10)
 *   j      TN(2)             W(8)    x      E(5)       BN(11)
 *   j-1                             SW(7)  SE(6)
 *           i     i+1        i-1     i     i+1          i      i+1
 *  
 *  if(k%2==1 && j%2==1) 
 *   j+1    TW(0)  TE(1)      NW(3)  NE(4)              BW(9)   BE(10)
 *   j             TN(2)       W(8)   x     E(5)                BN(11)
 *   j-1                      SW(7)  SE(6)
 *          i-1     i         i-1     i     i+1         i-1      i
 */

int Lattice::get_adjacent_coord(const int src_coord, const unsigned index) {
  int i, j, k;
  coord_to_ijk(src_coord, i, j, k);
  switch (index) {
  case 0:
    if (k%2) {
      return linearCoordFast((j%2==0)?i:i-1, j+1, k+1); // TW 
    }
    return linearCoordFast(i, j, k+1); // TN
  case 1:
    if (k%2) {
      return linearCoordFast((j%2==0)?i+1:i, j+1, k+1); // TE
    }
    return linearCoordFast((j%2==0)?i+1:i, j-1, k+1); // TE
  case 2:
    if (k%2) {
      return linearCoordFast(i, j, k+1); // TS
    }
    return linearCoordFast((j%2==0)?i:i-1, j-1, k+1); // TW
  case 3:
    return linearCoordFast((j%2==0)?i:i-1, j+1, k);  // NW
  case 4:
    return linearCoordFast((j%2==0)?i+1:i, j+1, k);  // NE
  case 5:
    return linearCoordFast(i+1, j, k);  // E 
  case 6:
    return linearCoordFast((j%2==0)?i+1:i, j-1, k);  // SE
  case 7:
    return linearCoordFast((j%2==0)?i:i-1, j-1, k);  // SW
  case 8:
    return linearCoordFast(i-1, j, k);  // W
  case 9:
    if (k%2) {
      return linearCoordFast((j%2==0)?i:i-1, j+1, k-1);// BW
    }
    return linearCoordFast(i, j, k-1);// BN
  case 10:
    if (k%2) {
      return linearCoordFast((j%2==0)?i+1:i, j+1, k-1);// BE
    }
    return linearCoordFast((j%2==0)?i+1:i, j-1, k-1);// BE
  }
  //case 11:
  if (k%2) {
    return linearCoordFast(i, j, k-1);// BS
  }
  return linearCoordFast((j%2==0)?i:i-1, j-1, k-1);// BW
}


/*
 *   three dimensional lock-step communication:
 *
 *   X-phase
 *     if(xbit==1) get (   1, ybegin:yend, zbegin:zend) from x+1 node
 *                  as (Nx+1, ybegin:yend, zbegin:zend) of myself
 *     if(xbit==0) get (  Nx, ybegin:yend, zbegin:zend) from x-1 node
 *                  as (   0, ybegin:yend, zbegin:zend) of myself
 *
 *   Y-phase
 *     if(ybit==1) get (xbegin:xend,    1, zbegin:zend) from y+1 node
 *                  as (xbegin:xend, Ny+1, zbegin:zend) of myself
 *     if(ybit==0) get (xbegin:xend,   Ny, zbegin:zend) from y-1 node
 *                  as (xbegin:xend,    0, zbegin:zend) of myself
 *
 *
 *   copy from (*bit==1 ? plus : minus);
 *   xghost  1:Ny,   1:Nz
 *   yghost  0:Nx+1, 1:Nz
 *   zghost  0:Nx+1, 0:Ny+1
 *   
 *   if(xbit==1) get (   1,1:Ny,1:Nz) from right node
 *                as (Nx+1,1:Ny,1:Nz) of myself
 *   if(xbit==0) get (  Nx,1:Ny,1:Nz) from left node
 *                as (   0,1:Ny,1:Nz) of myself
 *  
 *   if(ybit==1) get (0:Nx+1,   1,1:Nz) from right node
 *                as (0:Nx+1,Ny+1,1:Nz) of myself
 *   if(ybit==0) get (0:Nx+1,  Ny,1:Nz) from left node
 *                as (0:Nx+1,   0,1:Nz) of myself
 *  
 *   if(zbit==1) get (0:Nx+1,0:Ny+1,   1) from right node
 *                as (0:Nx+1,0:Ny+1,Nz+1) of myself
 *   if(zbit==0) get (0:Nx+1,0:Ny+1,  Nz) from left node
 *                as (0:Nx+1,0:Ny+1,   0) of myself
 *
 *   CAUTION: following function depend on the condition GHOST_SIZE=1 !!!
 *            setting comp_ID=-1 in clear-ghost causes semantic error !!! 
 */

void Lattice::clear_ghost(std::vector<SpillMolecule>& spill_coords) {
  for (unsigned i(0); i < spill_coords.size(); ++i) {
    int& species_id(voxelVector_[spill_coords[i].coord].species_id);
    if (species_id > -out_id_) {
      species_id = ghost_id_;
    }
    else { //out ghost voxel
      species_id = (species_id/out_id_)*out_id_;
    }
  }

  for (unsigned i(0); i < occupied_ghosts_.size(); ++i) {
    int& species_id(voxelVector_[occupied_ghosts_[i]].species_id);
    if (species_id > -out_id_) {
      species_id = ghost_id_;
    }
    else { //out ghost voxel
      species_id = (species_id/out_id_)*out_id_;
    }
  }
}

void Lattice::set_out_properties(Vector<int> begin[], Vector<int> end[],
                                 Vector<int> bit[]) {
  for (unsigned i(0); i < 8; ++i) {
    Boundary& b(boundary_[i]);
    b.begin = begin[i];
    b.end = end[i];
    b.bit = bit[i];
    b.xout = (bit[i].x == 1 ? 1:nx_);
    b.xin = (bit[i].x == 1 ? nx_+1:0);
    b.yout = (bit[i].y ==1 ? 1:ny_);
    b.yin = (bit[i].y ==1 ? ny_+1:0);
    b.zout = (bit[i].z == 1 ? 1:nz_);
    b.zin = (bit[i].z == 1 ? nz_+1:0);
    b.yx_first = (begin[i].y-2+ny_)%ny_+1;
    b.yx_last = (end[i].y+ny_)%ny_+1;
    b.zx_first = (begin[i].z-2+nz_)%nz_+1;
    b.zx_last = (end[i].z+nz_)%nz_+1;
  }
}

/*
void Lattice::set_out_voxels(const int &xbegin, const int &xend,
                             const int &ybegin, const int &yend,
                             const int &zbegin, const int &zend,
                             const int &xbit, const int &ybit,
                             const int &zbit, ParallelEnvironment &pe,
                             const int sub) {
  int Nx = pe.getispan(); 
  int Ny = pe.getjspan(); 
  int Nz = pe.getkspan(); 

  //x direction:
  const int i0 = ( xbit==1 ?    1 : Nx );
  const int i1 = ( xbit==1 ? Nx+1 :  0 );
  for (int jj=ybegin-1; jj<=yend+1; ++jj) {
    for (int kk=zbegin-1; kk<=zend+1; ++kk) {
      const int j = (jj-1+Ny)%Ny+1;
      const int k = (kk-1+Nz)%Nz+1;
      const int lc0 = linearCoordFast(i0,j,k);
      int& species_id(voxelVector_[lc0].species_id);
        if (species_id != invalid_id_) {
          if (species_id == ghost_id_) {  //out ghost voxel:
            species_id = -out_id_*(1+sub);
          }
          else if (species_id == vacant_id_) { //out normal voxel
            species_id = out_id_*(1+sub);
          }
          //shared out ghost voxel
          else if (species_id < 0 && species_id != -out_id_*(1+sub)) {
            species_id = -9*out_id_;
          }
          //shared out normal voxel
          else if (species_id > 0 && species_id != out_id_*(1+sub)) {
            species_id = 9*out_id_;
          }
          else {
            std::cout << "species_id:" << species_id << std::endl;
          }
        }
    }
  }

  //y direction:
  const int j0 = ( ybit==1 ?    1 : Ny );
  const int j1 = ( ybit==1 ? Ny+1 :  0 );
  for (int ii=xbegin-1; ii<=xend+1; ++ii) {
    for (int kk=zbegin-1; kk<=zend+1; ++kk) {
      const int i = ii;
      const int k = (kk-1+Nz)%Nz+1;
      const int lc0 = linearCoordFast(i,j0,k);
      int& species_id(voxelVector_[lc0].species_id);
        if (species_id != invalid_id_) {
          if (species_id == ghost_id_) {  //out ghost voxel:
            species_id = -out_id_*(1+sub);
          }
          else if (species_id == vacant_id_) { //out normal voxel
            species_id = out_id_*(1+sub);
          }
          //shared out ghost voxel
          else if (species_id < 0 && species_id != -out_id_*(1+sub)) {
            species_id = -9*out_id_;
          }
          //shared out normal voxel
          else if (species_id > 0 && species_id != out_id_*(1+sub)) {
            species_id = 9*out_id_;
          }
          else {
            std::cout << "species_id:" << species_id << std::endl;
          }
        }
    }
  }

  //z direction:
  const int k0 = ( zbit==1 ?    1 : Nz );  // donar
  const int k1 = ( zbit==1 ? Nz+1 :  0 );  // ghost
  for (int ii=xbegin-1; ii<=xend+1; ++ii) {
    for (int jj=ybegin-1; jj<=yend+1; ++jj) {
      const int i = ii;
      const int j = jj;
      const int lc0 = linearCoordFast(i,j,k0);
      int& species_id(voxelVector_[lc0].species_id);
        if (species_id != invalid_id_) {
          if (species_id == ghost_id_) {  //out ghost voxel:
            species_id = -out_id_*(1+sub);
          }
          else if (species_id == vacant_id_) { //out normal voxel
            species_id = out_id_*(1+sub);
          }
          //shared out ghost voxel
          else if (species_id < 0 && species_id != -out_id_*(1+sub)) {
            species_id = -9*out_id_;
          }
          //shared out normal voxel
          else if (species_id > 0 && species_id != out_id_*(1+sub)) {
            species_id = 9*out_id_;
          }
          else {
            std::cout << "species_id:" << species_id << std::endl;
          }
        }
    }
  }
}
*/
int Lattice::check_outmolecules(const std::vector<OutMolecule> outmolecules[],
                                const int sub) {
  Boundary& b(boundary_[sub]);
  for (unsigned n(1); n < 9; ++n) { 
    for (unsigned x(0); x < outmolecules[n].size(); ++x) { 
      int i, j, k;
      coord_to_ijk(outmolecules[n][x].coord, i, j, k);
      if (i == b.xout &&
          ((j >= b.begin.y && j <= b.end.y) || 
           j == b.yx_first || j == b.yx_last) &&
          ((k >= b.begin.z && k <= b.end.z) || 
           k == b.zx_first || k == b.zx_last)) {
        if (n != sub+1 && n != 9) {
          std::cout << "error in outmolecules xn:" << n << " sub:" << sub+1 <<
            std::endl;
          return 1;
        }
      }
      if (j == b.yout && i >= b.begin.x-1 && i <= b.end.x+1 &&
          ((k >= b.begin.z && k <= b.end.z) || 
           k == b.zx_first || k == b.zx_last)) {
        if (n != sub+1 && n != 9) {
          std::cout << "error in outmolecules xn:" << n << " sub:" << sub+1 <<
            std::endl;
          return 1;
        }
      }
      if (k == b.zout && i >= b.begin.x-1 && i <= b.end.x+1 &&
          j >= b.begin.y-1 && j <= b.end.y+1) { 
        if (n != sub+1 && n != 9) {
          std::cout << "error in outmolecules xn:" << n << " sub:" << sub+1 <<
            std::endl;
          return 1;
        }
      }
    }
  }
  return 0;
}

//set voxels that will send out its occupancy state to neighbour processes:
void Lattice::set_out_voxels() {
  for (int sub(0); sub < 8; ++sub) {
    const Boundary& b(boundary_[sub]);

    //x direction:
    for (int jj(b.begin.y-1); jj <= b.end.y+1; ++jj) {
      for (int kk(b.begin.z-1); kk <= b.end.z+1; ++kk) {
        const int j((jj-1+ny_)%ny_+1);
        const int k((kk-1+nz_)%nz_+1);
        const int lc0(linearCoordFast(b.xout, j, k));
        int& species_id(voxelVector_[lc0].species_id);
        if (species_id != invalid_id_) {
          if (species_id == ghost_id_) {  //out ghost voxel:
            species_id = -out_id_*(1+sub);
          }
          else if (species_id == vacant_id_) { //out normal voxel
            species_id = out_id_*(1+sub);
          }
          //shared out ghost voxel
          else if (species_id < 0 && species_id != -out_id_*(1+sub)) {
            species_id = -9*out_id_;
          }
          //shared out normal voxel
          else if (species_id > 0 && species_id != out_id_*(1+sub)) {
            species_id = 9*out_id_;
          }
        }
      }
    }

    //y direction:
    for (int i(b.begin.x-1); i <= b.end.x+1; ++i) {
      for (int kk(b.begin.z-1); kk <= b.end.z+1; ++kk) {
        const int k((kk-1+nz_)%nz_+1);
        const int lc0(linearCoordFast(i, b.yout, k));
        int& species_id(voxelVector_[lc0].species_id);
        if (species_id != invalid_id_) {
          if (species_id == ghost_id_) {  //out ghost voxel:
            species_id = -out_id_*(1+sub);
          }
          else if (species_id == vacant_id_) { //out normal voxel
            species_id = out_id_*(1+sub);
          }
          //shared out ghost voxel
          else if (species_id < 0 && species_id != -out_id_*(1+sub)) {
            species_id = -9*out_id_;
          }
          //shared out normal voxel
          else if (species_id > 0 && species_id != out_id_*(1+sub)) {
            species_id = 9*out_id_;
          }
        }
      }
    }

    //z direction:
    for (int i(b.begin.x-1); i <= b.end.x+1; ++i) {
      for (int j(b.begin.y-1); j <= b.end.y+1; ++j) {
        const int lc0 = linearCoordFast(i, j, b.zout);
        int& species_id(voxelVector_[lc0].species_id);
        if (species_id != invalid_id_) {
          if (species_id == ghost_id_) {  //out ghost voxel:
            species_id = -out_id_*(1+sub);
          }
          else if (species_id == vacant_id_) { //out normal voxel
            species_id = out_id_*(1+sub);
          }
          //shared out ghost voxel
          else if (species_id < 0 && species_id != -out_id_*(1+sub)) {
            species_id = -9*out_id_;
          }
          //shared out normal voxel
          else if (species_id > 0 && species_id != out_id_*(1+sub)) {
            species_id = 9*out_id_;
          }
        }
      }
    }
  }

  for (int sub(0); sub < 8; ++sub) {
    const Boundary& b(boundary_[sub]);

    //x direction:
    for (int jj(b.begin.y-1); jj <= b.end.y+1; ++jj) {
      for (int kk(b.begin.z-1); kk <= b.end.z+1; ++kk) {
        const int j((jj-1+ny_)%ny_+1);
        const int k((kk-1+nz_)%nz_+1);
        const int lc0(linearCoordFast(b.xout, j, k));
        int species_id(voxelVector_[lc0].species_id);
        if (species_id == 9*out_id_) {
          is_shared_[sub] = 1;
          break;
        }
      }
    }

    //y direction:
    for (int i(b.begin.x-1); i <= b.end.x+1; ++i) {
      for (int kk(b.begin.z-1); kk <= b.end.z+1; ++kk) {
        const int k((kk-1+nz_)%nz_+1);
        const int lc0(linearCoordFast(i, b.yout, k));
        int species_id(voxelVector_[lc0].species_id);
        if (species_id == 9*out_id_) {
          is_shared_[sub] = 1;
          break;
        }
      }
    }

    //z direction:
    for (int i(b.begin.x-1); i <= b.end.x+1; ++i) {
      for (int j(b.begin.y-1); j <= b.end.y+1; ++j) {
        const int lc0 = linearCoordFast(i, j, b.zout);
        int species_id(voxelVector_[lc0].species_id);
        if (species_id == 9*out_id_) {
          is_shared_[sub] = 1;
          break;
        }
      }
    }
  }
}

//all inbounds are confirmed not invalid_id_ and
//are all ghost voxels (species_id < 0)
void Lattice::load_ghost(ParallelEnvironment &pe,
                         const std::vector<OutMolecule>& outmolecules,
                         const std::vector<OutMolecule>& outmolecules_shared,
                         const int sub) { 
  occupied_ghosts_.clear();
  out_ghost_molecules_.clear();
  Boundary& b(boundary_[sub]);
  MPI::Cartcomm cart = pe.getcart();
  MPI::Request reqs[2];
  MPI::Status stat[2];

  int nsentx(0);
  int nsenty(0);
  int nsentz(0);

  for (unsigned x(0); x < outmolecules.size(); ++x) { 
    int i, j, k;
    coord_to_ijk(outmolecules[x].coord, i, j, k);
    if (i == b.xout) {
      outboundx[nsentx  ] = j;
      outboundx[nsentx+1] = k;
      outboundx[nsentx+2] = outmolecules[x].species_id;
      nsentx += 3; 
    }
    else if (j == b.yout) {
      outboundy[nsenty  ] = i;
      outboundy[nsenty+1] = k;
      outboundy[nsenty+2] = outmolecules[x].species_id;
      nsenty += 3;
    }
    else {
      outboundz[nsentz  ] = i;
      outboundz[nsentz+1] = j;
      outboundz[nsentz+2] = outmolecules[x].species_id;
      nsentz += 3;
    }
  }

  //if (is_shared_[sub]) {
    for (unsigned x(0); x < outmolecules_shared.size(); ++x) { 
      const int coord(outmolecules_shared[x].coord);
      if (getXind(coord) == b.xout) {
        outboundx[nsentx  ] = getYind(coord);
        outboundx[nsentx+1] = getZind(coord);
        outboundx[nsentx+2] = outmolecules_shared[x].species_id;
        nsentx += 3; 
      }
      else if (getYind(coord) == b.yout) {
        outboundy[nsenty  ] = getXind(coord);
        outboundy[nsenty+1] = getZind(coord);
        outboundy[nsenty+2] = outmolecules_shared[x].species_id;
        nsenty += 3;
      }
      else if (getZind(coord) == b.zout) {
        outboundz[nsentz  ] = getXind(coord);
        outboundz[nsentz+1] = getYind(coord);
        outboundz[nsentz+2] = outmolecules_shared[x].species_id;
        nsentz += 3;
      }
    }
  //}


  {
    //x direction:
    const int left  = pe.getneighbor_i_minus();
    const int right = pe.getneighbor_i_plus();
    const int to   = ( b.bit.x==1 ?  left : right );
    const int from = ( b.bit.x==1 ? right :  left );
    cart.Sendrecv(outboundx, nsentx, MPI::INT, to, 0, inbound, bufsize,
                  MPI::INT, from, 0, stat[1]); 

    // unpack received data
    const int nrecv = stat[1].Get_count( MPI::INT );
    for (int irecv=0; irecv<nrecv; irecv+=3) {
      const int j = inbound[irecv  ];
      const int k = inbound[irecv+1];
      const int lc1 = linearCoordFast(b.xin, j, k);
      occupied_ghosts_.push_back(lc1);
      int& species_id(voxelVector_[lc1].species_id);
      //if normal ghost:
      if (species_id == ghost_id_) {
        species_id = -inbound[irecv+2];
      }
      else if (species_id <= -out_id_) { //out ghost:
        species_id = (species_id/out_id_)*out_id_-inbound[irecv+2];
        out_ghost_molecules_.push_back(SpeciesCoord(inbound[irecv+2], lc1));
      }
    }
  }

  {
    //y direction:
    for (unsigned x(0); x < out_ghost_molecules_.size(); ++x) { 
      int i, j, k;
      coord_to_ijk(out_ghost_molecules_[x].coord, i, j, k);
      if (j == b.yout) {/* && i >= b.begin.x-1 && i <= b.end.x+1 &&
          ((k >= b.begin.z && k <= b.end.z) || 
           k == b.zx_first || k == b.zx_last)) {*/
        outboundy[nsenty  ] = i;
        outboundy[nsenty+1] = k;
        outboundy[nsenty+2] = out_ghost_molecules_[x].species_id;
        nsenty += 3;
      }
    }

    const int left  = pe.getneighbor_j_minus();
    const int right = pe.getneighbor_j_plus();
    const int to   = ( b.bit.y==1 ?  left : right );
    const int from = ( b.bit.y==1 ? right :  left );
    cart.Sendrecv( outboundy,   nsenty, MPI::INT, to,   0,
                    inbound, bufsize, MPI::INT, from, 0, stat[1] ); 
   
    // unpack received data
    const int nrecv = stat[1].Get_count( MPI::INT );
    for(int irecv=0; irecv<nrecv; irecv+=3) {
      const int i = inbound[irecv  ];
      const int k = inbound[irecv+1];
      const int lc1 = linearCoordFast(i, b.yin, k);
      occupied_ghosts_.push_back(lc1);
      int& species_id(voxelVector_[lc1].species_id);
      //if normal ghost:
      if (species_id == ghost_id_) {
        species_id = -inbound[irecv+2];
      }
      else if (species_id <= -out_id_) { //out ghost:
        species_id = (species_id/out_id_)*out_id_-inbound[irecv+2];
        out_ghost_molecules_.push_back(SpeciesCoord(inbound[irecv+2], lc1));
      }
    }
  }


  {
    //z direction:
    for (unsigned x(0); x < out_ghost_molecules_.size(); ++x) {
      int i, j, k;
      coord_to_ijk(out_ghost_molecules_[x].coord, i, j, k);
      if (k == b.zout) {/* && i >= b.begin.x-1 && i <= b.end.x+1 &&
          j >= b.begin.y-1 && j <= b.end.y+1) {*/
        outboundz[nsentz  ] = i;
        outboundz[nsentz+1] = j;
        outboundz[nsentz+2] = out_ghost_molecules_[x].species_id;
        nsentz += 3;
      }
    }

    // send and receive
    const int left  = pe.getneighbor_k_minus();
    const int right = pe.getneighbor_k_plus();
    const int to   = ( b.bit.z==1 ?  left : right );
    const int from = ( b.bit.z==1 ? right :  left );
    cart.Sendrecv( outboundz, nsentz, MPI::INT, to,   0,
                    inbound, bufsize, MPI::INT, from, 0, stat[1] ); 
   
    // unpack received data
    const int nrecv = stat[1].Get_count( MPI::INT ); 
    for(int irecv=0; irecv<nrecv; irecv+=3) {
      const int i = inbound[irecv  ];
      const int j = inbound[irecv+1];
      const int lc1 = linearCoordFast(i, j, b.zin);
      occupied_ghosts_.push_back(lc1);
      int& species_id(voxelVector_[lc1].species_id);
      //if normal ghost:
      if (species_id == ghost_id_) {
        species_id = -inbound[irecv+2];
      }
      else if (species_id <= -out_id_) { //out ghost:
        species_id = (species_id/out_id_)*out_id_-inbound[irecv+2];
      }
    }
  }
}

void Lattice::feedGhost(const int &xbegin, const int &xend,
                        const int &ybegin, const int &yend,
                        const int &zbegin, const int &zend,
                        const int &xbit,
                        const int &ybit,
                        const int &zbit, ParallelEnvironment &pe) {
  // dimension without boundary
  const int Nx = pe.getispan(); 
  const int Ny = pe.getjspan(); 
  const int Nz = pe.getkspan(); 
  MPI::Cartcomm cart = pe.getcart();
  MPI::Request  reqs[2];
  MPI::Status   stat[2];

  // identify ghost region
  const int xghost = (xbit==0 ? 0 : Nx+1);
  const int yghost = (ybit==0 ? 0 : Ny+1);
  const int zghost = (zbit==0 ? 0 : Nz+1);

  //z direction:
  {
    int nsent = 0;
    const int nspill = spillcoordsZ.size();
    for(int m=0; m<nspill; ++m) {
      const int lc = spillcoordsZ[m].coord;
      outboundx[nsent  ] = getXind(lc);
      outboundx[nsent+1] = getYind(lc);
      outboundx[nsent+2] = spillcoordsZ[m].species_id;
      outboundx[nsent+3] = spillcoordsZ[m].mol_id;
      nsent += 4;
    }

    // send and receive
    const int left  = pe.getneighbor_k_minus();
    const int right = pe.getneighbor_k_plus();
    const int to   = ( zbit==1 ? right :  left );   // upside-down
    const int from = ( zbit==1 ?  left : right );   // upside-down
    cart.Sendrecv(outboundx,   nsent, MPI::INT, to,   0,
                  inbound, bufsize, MPI::INT, from, 0, stat[1] ); 

    // unpack received data
    const int k1 = ( zbit==1 ? 1 : Nz );  // target
    const int nrecv = stat[1].Get_count( MPI::INT ); 
    for(int irecv=0; irecv<nrecv; irecv+=4) {
      const int i = inbound[irecv  ];
      const int j = inbound[irecv+1];
      const int species_id = inbound[irecv+2];
      const unsigned mol_id = inbound[irecv+3];
      const int lc = linearCoordFast(i,j,k1);

      if (j == yghost) {
        spillcoordsY.push_back(SpillMolecule(species_id, lc, mol_id));
      }
      else if (i == xghost) {
        spillcoordsX.push_back(SpillMolecule(species_id, lc, mol_id));
      }
      else {
        jumpincoords.push_back(SpillMolecule(species_id, lc, mol_id));
      }
    }
  }

  {
    //y direction:
    int nsent = 0;
    const int nspill = spillcoordsY.size();
    for(int m=0; m<nspill; ++m) {
      const int lc = spillcoordsY[m].coord;
      outboundx[nsent  ] = getXind(lc);
      outboundx[nsent+1] = getZind(lc);
      outboundx[nsent+2] = spillcoordsY[m].species_id;
      outboundx[nsent+3] = spillcoordsY[m].mol_id;
      nsent += 4;
    }

    // send and receive
    const int left  = pe.getneighbor_j_minus();
    const int right = pe.getneighbor_j_plus();
    const int to   = ( ybit==1 ? right :  left );   // leftside-right
    const int from = ( ybit==1 ?  left : right );   // leftside-right
    cart.Sendrecv( outboundx,   nsent, MPI::INT, to,   0,
                    inbound, bufsize, MPI::INT, from, 0, stat[1] ); 
    // unpack received data
    const int j1 = ( ybit==1 ? 1 : Ny );  // target
    const int nrecv = stat[1].Get_count( MPI::INT );
    for(int irecv=0; irecv<nrecv; irecv+=4) {
      const int i = inbound[irecv  ];
      const int k = inbound[irecv+1];
      const int species_id = inbound[irecv+2];
      const unsigned mol_id = inbound[irecv+3];
      const int lc = linearCoordFast(i,j1,k);

      if (i == xghost) {
        spillcoordsX.push_back(SpillMolecule(species_id, lc, mol_id));
      }
      else {
        jumpincoords.push_back(SpillMolecule(species_id, lc, mol_id));
      }
    }
  }

  {
    //x direction:
    int nsent = 0;
    const int nspill = spillcoordsX.size();
    for(int m=0; m<nspill; ++m) {
      const int lc = spillcoordsX[m].coord;
      outboundx[nsent  ] = getYind(lc);
      outboundx[nsent+1] = getZind(lc);
      outboundx[nsent+2] = spillcoordsX[m].species_id;
      outboundx[nsent+3] = spillcoordsX[m].mol_id;
      nsent += 4;
    }

    // send and receive
    const int left  = pe.getneighbor_i_minus();
    const int right = pe.getneighbor_i_plus();
    const int to   = ( xbit==1 ? right :  left );   // backside-fore
    const int from = ( xbit==1 ?  left : right );   // backside-fore
    cart.Sendrecv( outboundx,   nsent, MPI::INT, to,   0,
                    inbound, bufsize, MPI::INT, from, 0, stat[1] ); 

    // unpack received data
    const int i1 = ( xbit==1 ? 1 : Nx );  // target
    const int nrecv = stat[1].Get_count( MPI::INT );
    for(int irecv=0; irecv<nrecv; irecv+=4) {
      const int j = inbound[irecv  ];
      const int k = inbound[irecv+1];
      const int species_id = inbound[irecv+2];
      const unsigned mol_id = inbound[irecv+3];
      const int lc = linearCoordFast(i1,j,k);
      jumpincoords.push_back(SpillMolecule(species_id, lc, mol_id));
    }
  }
}

