#ifndef __LATTICE_HPP
#define __LATTICE_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include "Vector.hpp"
#include "ParallelEnvironment.hpp"
#include "Molecule.hpp"
#include "Species.hpp"
#include "Common.hpp"
using namespace std;

struct Boundary {
  Vector<int> begin;
  Vector<int> end;
  Vector<int> bit;
  int xout;
  int xin;
  int yout;
  int yin;
  int zout;
  int zin;
  int yx_first;
  int yx_last;
  int zx_first;
  int zx_last;
};

class Lattice
{
public:
    Lattice(string s, double r, ParallelEnvironment &pe,
            const int invalid_id, const int vacant_id, const int ghost_id);
  
    ~Lattice()
    {
       delete []  inbound;
       delete [] outboundx;
       delete [] outboundy;
       delete [] outboundz;
    }

    string getName()
    { 
        return name_;
    }
   
    int getXdim()
    {
        return LNx_; 
    }

    int getYdim()
    {
        return LNy_;
    }

    int getZdim()
    { 
        return LNz_;
    }

    double getradius()
    { 
        return radius_;
    }
  
    Voxel& get_voxel(int index)
    {
        return voxels_[index]; 
    }
  
    int linearCoord(int i, int j, int k) 
    {
       return ( (i<0 || i>=LNx_ ||
                 j<0 || j>=LNy_ ||
                 k<0 || k>=LNz_) ? -1 : k + ( j + i * LNy_ ) * LNz_ );
    }
    // convert cartesian to linear coordinates without range check
    int linearCoordFast(int i, int j, int k) 
    {
        return k + ( j + i * LNy_ ) * LNz_;
    }
  
    // convert cartesian to linear coordinates (overloaded)
    int linearCoord(Coordinate position) 
    {
        return position.ck + ( position.cj + position.ci * LNy_ ) * LNz_;
    }

    int get_adjacent_coord(const int src_coord, const unsigned index);
  
    int getXind(int index) 
    {
        return index / (LNz_ * LNy_); 
    }

    int getYind(int index)
    {
        return (index % (LNz_ * LNy_)) / LNz_; 
    }

    int getZind(int index)
    { 
        return index % LNz_; 
    }

    void coord_to_ijk(const int coord, int& i, int& j, int& k) {
        i = coord / (LNz_ * LNy_); 
        j = (coord % (LNz_ * LNy_)) / LNz_; 
        k = coord % LNz_; 
    }
    void clear_ghost(std::vector<SpillMolecule>& spill_coords);

    // real-x = rv*i*2 - rv*(j%2) 
    double getXpos(int i, int j, int k) 
    {
        return radius_ * ( i * 2 - (j % 2) );
    }
  
    // real-y = rv*j*sqrt(3) + 2*rv*sqrt(1/3)*(k%2)  
    double getYpos(int i, int j, int k)
    {
        return radius_ * ( j * SQR31_ + (k % 2) * 2 * SQR13_ );
    }
  
    // real-z = rv*k*sqrt(8/3);
    double getZpos(int i, int j, int k) 
    {
        return radius_ * ( k * SQR83_ );
    }
  
    int getID()
    {
        return latticeID_; 
    } 
  
    void setID(int ID)
    {
        latticeID_ = ID;
    }   

    void set_out_id(const int out_id) {
      out_id_ = out_id;
    }

    void set_out_properties(Vector<int> begin[], Vector<int> end[],
                            Vector<int> bit[]);
    void set_out_voxels();
    void set_out_voxels(const int &xbegin, const int &xend,
                        const int &ybegin, const int &yend,
                        const int &zbegin, const int &zend,
                        const int &xbit, const int &ybit,
                        const int &zbit, ParallelEnvironment &pe,
                        const int sub);

    void load_ghost(ParallelEnvironment &pe,
                    const std::vector<OutMolecule>& outmolecules,
                    const std::vector<OutMolecule>& outmolecules_shared,
                    const int sub);

    void feedGhost(const int &xbegin, const int &xend, 
                   const int &ybegin, const int &yend, 
                   const int &zbegin, const int &zend, 
                   const int &xbit,   const int &ybit,   const int &zbit,
                   ParallelEnvironment &pe);
    int check_outmolecules(const std::vector<OutMolecule> outmolecules[],
                            const int sub);

    vector<SpillMolecule> jumpincoords;

    // followings are used to update Lattice::voxels_      
    vector<SpillMolecule> spillcoordsX;
    vector<SpillMolecule> spillcoordsY;
    vector<SpillMolecule> spillcoordsZ;

private:
  const int invalid_id_;
  const int vacant_id_;
  const int ghost_id_;
  const string name_;
  const int nx_;          // comp dimensions
  const int ny_;
  const int nz_;
  const int LNx_;          // system dimensions
  const int LNy_;
  const int LNz_;
  const double radius_;
  const double SQR31_;
  const double SQR13_;
  const double SQR83_;
  int out_id_;
  int latticeID_;
  std::vector<Voxel> voxels_;
  std::vector<int> occupied_ghosts_;
  std::vector<SpeciesCoord> out_ghost_molecules_;
  int *outboundx;
  int *outboundy;
  int *outboundz;
  int  *inbound;
  int bufsize;  // size of in/outbound
  Boundary boundary_[8];
  std::vector<OutMolecule> om[10];
  unsigned is_shared_[8] = {0, 0, 0, 0, 0, 0, 0, 0};
};



#endif /* __LATTICE_HPP */
