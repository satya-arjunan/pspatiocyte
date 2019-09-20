#include <limits>
#include "Molecule.hpp"
#include "Compartment.hpp"
#include "World.hpp"


//Before calling this partition method, Lattice::Lattice has already set
//all edge voxels of this process' domain to ghost. Depending on the location
//of each domain, the ghost voxels can be invalid. We make them invalid here
//in this partition method below. So after calling this method, some
//ghost voxels will become invalid.
void Compartment::initialize(Lattice &g, ParallelEnvironment &pe,
                             std::vector<Species*>& species_list) {
  species_list_ = species_list;
  species_molecules_.resize(species_list_.size());
  Nx_ = pe.getispan();    // 0 and Nx+1 serve as ghost cells
  Ny_ = pe.getjspan();    // 0 and Ny+1 serve as ghost cells
  Nz_ = pe.getkspan();    // 0 and Nz+1 serve as ghost cells 
  const int NE = MPI::PROC_NULL;  // no neighbor
  const double r = g.getradius();
  const double SQR2 = 1.414213562;

  //set all invalid ids:
  if (type_==VOLUME) {
    // if right-boundary node
    if(pe.getneighbor_i_plus()==NE)
    for(int j=0; j<=Ny_+1; ++j)
    for(int k=0; k<=Nz_+1; ++k)
    {
        const int lc0 = g.linearCoord(Nx_  , j, k);
        const int lc1 = g.linearCoord(Nx_+1, j, k);
        g.getVoxel(lc0).species_id =
        g.getVoxel(lc1).species_id = invalid_id_;
    }

    // if left-boundary node
    if(pe.getneighbor_i_minus()==NE)
    for(int j=0; j<=Ny_+1; ++j)
    for(int k=0; k<=Nz_+1; ++k)
    {
        const int lc0 = g.linearCoord(0, j, k);
        const int lc1 = g.linearCoord(1, j, k);
        g.getVoxel(lc0).species_id =
        g.getVoxel(lc1).species_id = invalid_id_;
    }

    // if rear-boundary node
    if(pe.getneighbor_j_plus()==NE)
    for(int i=0; i<=Nx_+1; ++i)
    for(int k=0; k<=Nz_+1; ++k)
    {
        const int lc0 = g.linearCoord(i, Ny_  , k);
        const int lc1 = g.linearCoord(i, Ny_+1, k);
        g.getVoxel(lc0).species_id =
        g.getVoxel(lc1).species_id = invalid_id_;
    }

    // if front-boundary node
    if(pe.getneighbor_j_minus()==NE)
    for(int i=0; i<=Nx_+1; ++i)
    for(int k=0; k<=Nz_+1; ++k)
    {
        const int lc0 = g.linearCoord(i, 0, k);
        const int lc1 = g.linearCoord(i, 1, k);
        g.getVoxel(lc0).species_id =
        g.getVoxel(lc1).species_id = invalid_id_;
    }

    // if top boundary node
    if(pe.getneighbor_k_plus()==NE)
    for(int i=0; i<=Nx_+1; ++i)
    for(int j=0; j<=Ny_+1; ++j)
    {
        const int lc0 = g.linearCoord(i, j, Nz_  );
        const int lc1 = g.linearCoord(i, j, Nz_+1);
        g.getVoxel(lc0).species_id =
        g.getVoxel(lc1).species_id = invalid_id_;
    }

    // if bottom boundary node
    if(pe.getneighbor_k_minus()==NE)
    for(int i=0; i<=Nx_+1; ++i)
    for(int j=0; j<=Ny_+1; ++j)
    {
        const int lc0 = g.linearCoord(i, j, 0);
        const int lc1 = g.linearCoord(i, j, 1);
        g.getVoxel(lc0).species_id =
        g.getVoxel(lc1).species_id = invalid_id_;
    }

    // append inner voxels
    for(int i=1; i<=Nx_; ++i)
    for(int j=1; j<=Ny_; ++j)
    for(int k=1; k<=Nz_; ++k)
    {
        const int lc = g.linearCoord(i, j, k);
        if(g.getVoxel(lc).species_id == vacant_id_)
            voxelVector_.push_back(lc);
    }
    numberOfVoxel_ = voxelVector_.size();
    local_volume_ = numberOfVoxel_*4.0*SQR2*r*r*r;
    volume_ = local_volume_;
    const MPI::Cartcomm cart = pe.getcart();
    cart.Allreduce(MPI::IN_PLACE, &volume_, 1, MPI::DOUBLE, MPI::SUM);
  }

  mid_span_.x = Nx_/2 + Nx_%2;
  mid_span_.y = Ny_/2 + Ny_%2;
  mid_span_.z = Nz_/2 + Nz_%2; 
  for (unsigned N(0); N < 8; ++N) {
    bit_[N].x = ((N/4)%2);
    bit_[N].y = ((N/2)%2);
    bit_[N].z = ((N/1)%2);
    begin_[N].x = (bit_[N].x*mid_span_.x + GHOST_SIZE);
    begin_[N].y = (bit_[N].y*mid_span_.y + GHOST_SIZE);
    begin_[N].z = (bit_[N].z*mid_span_.z + GHOST_SIZE);
    end_[N].x = (bit_[N].x*(Nx_-mid_span_.x)+mid_span_.x + GHOST_SIZE - 1);
    end_[N].y = (bit_[N].y*(Ny_-mid_span_.y)+mid_span_.y + GHOST_SIZE - 1);
    end_[N].z = (bit_[N].z*(Nz_-mid_span_.z)+mid_span_.z + GHOST_SIZE - 1);
    ghosts_[N].x = (bit_[N].x==0 ? begin_[N].x-1 : end_[N].x+1);
    ghosts_[N].y = (bit_[N].y==0 ? begin_[N].y-1 : end_[N].y+1);
    ghosts_[N].z = (bit_[N].z==0 ? begin_[N].z-1 : end_[N].z+1);
  }
  g.set_out_properties(begin_, end_, bit_);
  g.set_out_voxels();

  for (unsigned i(0); i < 10; ++i) {
    out_cnts_[i] = 0;
    out_ghost_cnts_[i] = 0;
  }
  for(int i=0; i<=Nx_+1; ++i) {
    for(int j=0; j<=Ny_+1; ++j) {
      for(int k=0; k<=Nz_+1; ++k) {
        const unsigned coord(g.linearCoord(i, j, k));
        Voxel& voxel(g.getVoxel(coord));
        int species_id(voxel.species_id);
        if (species_id >= out_id_) {
          ++n_out_voxels_;
          ++out_cnts_[species_id/out_id_];
        }
        if (species_id <= -out_id_) {
          ++n_out_ghost_voxels_;
          ++out_ghost_cnts_[abs(species_id/out_id_)];
        }
        if (species_id == ghost_id_) {
          ++n_ghost_voxels_;
        }
      }
    }
  }
  /*
  std::cout << "n_out:" << n_out_voxels_ << " n_out_ghost:" <<
    n_out_ghost_voxels_ << " n_ghost:" << n_ghost_voxels_ << std::endl;
    */
}

void Compartment::check_voxels(Lattice& g, const double id) {
  unsigned n_out_voxels(0);
  unsigned n_out_ghost_voxels(0);
  unsigned n_ghost_voxels(0);
  unsigned out_cnts[10];
  unsigned out_ghost_cnts[10];
  for (unsigned i(0); i < 10; ++i) {
    out_cnts[i] = 0;
    out_ghost_cnts[i] = 0;
  }
  for(int i=0; i<=Nx_+1; ++i) {
    for(int j=0; j<=Ny_+1; ++j) {
      for(int k=0; k<=Nz_+1; ++k) {
        const unsigned coord(g.linearCoord(i, j, k));
        Voxel& voxel(g.getVoxel(coord));
        int species_id(voxel.species_id);
        if (species_id >= out_id_) {
          ++n_out_voxels;
          ++out_cnts[species_id/out_id_];
        }
        if (species_id <= -out_id_) {
          ++n_out_ghost_voxels;
          ++out_ghost_cnts[abs(species_id/out_id_)];
        }
        if (species_id < 0 && species_id > -out_id_) {
          ++n_ghost_voxels;
        }
        if (species_id > out_id_ && species_id%out_id_ != 0) {
          const int sid(species_id);
          const int out_id(species_id/out_id_);
          species_id = species_id%out_id_;
          std::vector<Molecule>& mols(species_molecules_[species_id]);
          if (voxel.mol_index >= outmolecules_[out_id].size()) {
            std::cout << id << " out id:" << out_id << " species_id:" << sid <<
              " voxel.mol_index:" << voxel.mol_index << " size:" <<
              outmolecules_[out_id].size() <<
              std::endl;
            abort();
          }
          if (outmolecules_[out_id][voxel.mol_index].mol_index >= mols.size()) {
            std::cout << id << " a voxel.mol_index:" << voxel.mol_index << 
              " size:" << outmolecules_[out_id].size() << std::endl;
            std::cout << id << " outmolecules_[voxel.mol_index].mol_index:" << 
              outmolecules_[out_id][voxel.mol_index].mol_index << " size:" <<
              mols.size() << " out_id:" << out_id << " sid:" << sid <<
              " species_id:" << species_id << std::endl;
            abort();
          }
          if (mols[outmolecules_[out_id][voxel.mol_index].mol_index].coord != 
              coord) {
            std::cout << id << " error not consistent species id on out_voxel"
              << " coord:" << coord << " another coord:" << 
              mols[outmolecules_[out_id][voxel.mol_index].mol_index].coord <<
              std::endl;
              std::cout << "\tout_id:" << out_id << std::endl;
            abort();
          }

        }
        else if (species_id > 0 && species_id != invalid_id_ &&
                 species_id != vacant_id_ && species_id < out_id_) {
          std::vector<Molecule>& mols(species_molecules_[species_id]);
          if (voxel.mol_index >= mols.size()) {
            std::cout << id << " nonoutvoxel voxel.mol_index:" <<
              voxel.mol_index << " size:" << mols.size() << " species:" 
              << species_id << " name:" <<
              species_list_[species_id]->getName() << std::endl;
            abort();
          }
          if (mols[voxel.mol_index].coord != coord) {
            std::cout << id << 
              " error not consistent species id on normal voxel" <<
              " coord:" << coord << " other coord:" << 
              mols[voxel.mol_index].coord << " species id:" << species_id <<
              std::endl;
            abort();
          }
        }
      }
    }
  }
  if (n_out_voxels != n_out_voxels_) {
    std::cout << id << " error: counted n_out_voxels:" << n_out_voxels <<
      " not equals to init number:" << n_out_voxels_ <<
      "\n\tcounted out ghost:" << n_out_ghost_voxels <<
      " init out ghost:" << n_out_ghost_voxels_ << " counted ghost:" <<
      n_ghost_voxels << " init ghost:" << n_ghost_voxels_ << std::endl;
    abort();
  }
  if (n_out_ghost_voxels != n_out_ghost_voxels_) {
    std::cout << id << " error: counted n_out_ghost_voxels:" << 
      n_out_ghost_voxels << " not equals to init number:" << 
      n_out_ghost_voxels_ << "\n\tcounted out normal:" << n_out_voxels <<
      " init out normal:" << n_out_voxels_ << " counted ghost:" <<
      n_ghost_voxels << " init ghost:" << n_ghost_voxels_ << std::endl;
    abort();
  }
  if (n_ghost_voxels != n_ghost_voxels_) {
    std::cout << id << " error: counted n_ghost_voxels:" << 
      n_ghost_voxels << " not equals to init number:" << 
      n_ghost_voxels_ << "\n\tcounted out normal:" << n_out_voxels <<
      " init out normal:" << n_out_voxels_ << " counted out ghost:" <<
      n_out_ghost_voxels << " init out ghost:" << n_out_ghost_voxels_ <<
      std::endl;
    abort();
  }
  for (unsigned i(0); i < 10; ++i) {
    /*
    std::cout << "i:" << i << " " << out_cnts[i] << " " << out_cnts_[i] <<
      std::endl;
    std::cout << "i:" << i << " " << out_ghost_cnts[i] << " " <<
      out_ghost_cnts_[i] << std::endl;
      */
    if (out_cnts[i] != out_cnts_[i]) {
      std::cout << id << " error: counted out_cnts in " << i << ", " <<
        out_cnts[i] << " not equals to init number:" << out_cnts_[i] <<
        std::endl;
      abort();
    }
    if (out_ghost_cnts[i] != out_ghost_cnts_[i]) {
      std::cout << id << " error: counted out_ghost_cnts in " << i << ", " <<
        out_ghost_cnts[i] << " not equals to init number:" << 
        out_ghost_cnts_[i] << std::endl;
      abort();
    }
  }

  for (unsigned i(1); i < 10; ++i) {
    for (unsigned j(0); j < outmolecules_[i].size(); ++j) {
      if (outmolecules_[i][j].species_id%out_id_ == 0) {
        std::cout << "error in out molecules out id:" << i << std::endl;
        abort();
      }
    }
  }
  std::vector<OutMolecule> om[10];
  for(int i=0; i<=Nx_+1; ++i) {
    for(int j=0; j<=Ny_+1; ++j) {
      for(int k=0; k<=Nz_+1; ++k) {
        const unsigned coord(g.linearCoord(i, j, k));
        Voxel& voxel(g.getVoxel(coord));
        int species_id(voxel.species_id);
        if (species_id >= out_id_) {
          om[species_id/out_id_].push_back(OutMolecule(coord,0,0));
        }
      }
    }
  }
  for (unsigned i(0); i < 8; ++i) {
    if (g.check_outmolecules(om, i)) {
      std::cout << id << " error in coord of outmolecules" << std::endl;
      abort();
    }
  }
}

void Compartment::throwinMolecules(Species& s, const unsigned size,
                                   Lattice &g) {
  const int species_id(s.getID());
  if (size > numberOfVoxel_) {
    fout << "ERROR: too many molecules" << endl;
    abort();
  }
  std::vector<Molecule>& molecules(species_molecules_[species_id]);
  double rv = g.getradius();
  if (type_==VOLUME) {
    s.setVolumeDt(rv);
  }
  for (unsigned i(0); i < size; ++i) {
    unsigned trial = 0;
    unsigned index;
    unsigned coord;
    int sid;
    do {
      if(trial++ > numberOfVoxel_) {
        fout << "ABORT: hardly able to find vacant voxel! (" << s.getName() <<
          ")" << endl;
        abort();
      }
      index = (int)(numberOfVoxel_*(*randdbl_)());
      coord = voxelVector_[index];
      sid = g.getVoxel(coord).species_id;
    } while((sid < out_id_ && sid != vacant_id_) ||
            (sid > out_id_ && sid%out_id_ != 0));
    add_molecule(species_id, coord, get_new_mol_id(), g.getVoxel(coord)); 
  }
}

void Compartment::attachIndependentReaction(Reaction& r) {
  independent_reactions_.push_back(&r);
}

bool Compartment::processingReaction(Reaction &f, Lattice &g,
                                     ParallelEnvironment &pe) {
  const int r0 = f.getR0();       // moles of reactant 0
  const int r1 = f.getR1();       // moles of reactant 1
  const int r2 = f.getP2();       // moles of product  2
  const int r3 = f.getP3();       // moles of product  3
  Species *s0 = f.getS0();        // species of reactant 0
  Species *s1 = f.getS1();        // species of reactant 1
  Species *s2 = f.getS2();        // species of product  2
  Species *s3 = f.getS3();        // species of product  3

  if(r0 && r1) {
    fout << "ERROR: bimolecular independent reaction is not allowed"
      << ", reactant1 is omitted!" << endl;
    abort();
  }

  int center=-1;    // coordinate of voxel of which molecule is removed
  int neighbor=-1;  // coordinate of neighbor voxel for binary products
  bool result = false;  // whether reaction is fired
  const bool binary = (r2 && r3) ? true : false; // two products

  if (r0) { 
    result = eliminateReactant(*s0, g, binary, center, neighbor);
  }

  if (r1) {
    result = eliminateReactant(*s1, g, binary, center, neighbor);
  }

  if (result) {
    if (binary) { 
      const double location = (*randdbl_)();
      if(location>0.5) { 
        generateProduct(*s2, g, center);
        generateProduct(*s3, g, neighbor);
      }
      else {
        generateProduct(*s2, g, neighbor);
        generateProduct(*s3, g, center);
      }
    }
    else {
      if (r2) {
        generateProduct(*s2, g, center);
      }
      else
      if (r3) {
        generateProduct(*s3, g, center);
      }
    }
  }
  return result;
}

bool Compartment::eliminateReactant(Species &s, Lattice &g, bool binary,
                                    int &center, int &neighbor) {
  std::vector<Molecule>& molecules(species_molecules_[s.getID()]); 
  const int size(molecules.size());
  const int target((size==1) ? 0 : (int)(size*(*randdbl_)()));
  Molecule& molecule(molecules[target]);
  const int coord(molecule.coord);
  if(binary) {
    bool noneighbor(searchAllNeighbors(g, coord, center, neighbor));
    if (noneighbor) {
      return false;
    } 
  }
  else {
    center = coord;
  }
  Voxel& v0(g.getVoxel(coord));
  remove_molecule(g, s.getID(), v0);
  return true;
}


bool Compartment::searchAllNeighbors(Lattice &g, const int &coord, int &center,
                                     int &neighbor) {
  int order[12];
  bool intact[12] = {true, true, true, true, true, true, true, true, true,
    true, true, true};
  int slot = 0;
  while( slot<12 ) {
    const int pos = (int)(12.0*(*randdbl_)());
    if( intact[pos] ) {
      intact[pos] = false;
      order[slot++] = pos;
    }
  }

  bool nongb = true; // all neighbor voxels are not selected
  // search for all neighbor voxels
  for(int i=0; i<12; i++) {
    const int ngb = order[i];
    const int adj = g.get_adjacent_coord(coord, ngb);
    // if No coordinate (-1)
    if (adj < 0) {
      continue;
    }

    const int species_id(g.getVoxel(adj).species_id);
    if(species_id == vacant_id_ || (species_id >= out_id_ &&
                                    species_id%out_id_ == 0)) {
      center   = coord;
      neighbor = adj;
      nongb = false;
      break;
    }
  }
  return nongb;
}


void Compartment::generateProduct(Species &s, Lattice &g, int coord) {
  Voxel& voxel(g.getVoxel(coord));
  add_molecule(s.getID(), coord, get_new_mol_id(), voxel); 
}

void Compartment::calculatePropensityNoParallel(Reaction &f) {
  const int r0 = f.getR0();        // moles of reactant 0
  const int r1 = f.getR1();        // moles of reactant 1
  Species  *s0 = f.getS0();        // species of reactant 0
  Species  *s1 = f.getS1();        // species of reactant 1

  double propensity = 0.0;
  double concentration0 = 0.0;  // number of reactant S0 / volume
  double concentration1 = 0.0;  // number of reactant S1 / volume

  if(r0) {
    std::vector<Molecule>& molecules(species_molecules_[(*s0).getID()]); 
    concentration0 = molecules.size()/volume_;
  }

  if(r1) {
    std::vector<Molecule>& molecules(species_molecules_[(*s1).getID()]); 
    concentration1 = molecules.size()/volume_;
  }
  propensity = f.getK()*pow(concentration0, r0)*pow(concentration1, r1)*volume_;
  f.setOldPropensity(propensity);
}

void Compartment::calculateLocalPropensity(Reaction &f, double current_time) {
  const int r0 = f.getR0();        // moles of reactant 0
  const int r1 = f.getR1();        // moles of reactant 1
  Species  *s0 = f.getS0();        // species of reactant 0
  Species  *s1 = f.getS1();        // species of reactant 1

  double propensity = 0.0;
  double concentration0 = 0.0;  // number of reactant S0 / local_volume
  double concentration1 = 0.0;  // number of reactant S1 / local_volume

  if(r0) { 
    std::vector<Molecule>& molecules(species_molecules_[(*s0).getID()]); 
    concentration0 = molecules.size()/local_volume_;
  }

  if(r1) {
    std::vector<Molecule>& molecules(species_molecules_[(*s1).getID()]); 
    concentration1 = molecules.size()/local_volume_;
  }

  propensity = f.getK()*pow( concentration0, r0)*
    pow(concentration1, r1 )*local_volume_;
  f.setPropensity(propensity);
}


double Compartment::getNewPropensity(double current_time) {
  vector<Reaction*>::iterator p(independent_reactions_.begin());
  vector<Reaction*>::iterator e(independent_reactions_.end());
  totalPropensity_ = 0.0;
  while(p!=e) {
      calculateLocalPropensity(*(*p), current_time);
      totalPropensity_ += (*p)->getPropensity();
      p++;
  }
  vector<Reaction*>::iterator p3(independent_reactions_.begin());
  vector<Reaction*>::iterator e3(independent_reactions_.end());
  while(p3!=e3) {
      calculatePropensityNoParallel(*(*p3));
      p3++;
  }
  return totalPropensity_;
}

double Compartment::getOldTotalPropensity() {
  vector<Reaction*>::iterator p(independent_reactions_.begin());
  vector<Reaction*>::iterator e(independent_reactions_.end());
  double totalPropensity = 0.0;
  while(p!=e) {
    totalPropensity += (*p)->getOldPropensity();
    p++;
  }
  return totalPropensity;
}

double Compartment::getNextTime(ParallelEnvironment &pe, double current_time) {
  if(next_react_time_ == std::numeric_limits<double>::infinity() || 
     next_react_time_ == current_time) {
    return getNewNextTime(pe, current_time);
  }
  double dt(std::numeric_limits<double>::infinity());
  getNewPropensity(current_time);
  const double oldPropensity(globalTotalPropensity_);
  globalTotalPropensity_ = totalPropensity_;
  pe.getcart().Allreduce(MPI::IN_PLACE, &globalTotalPropensity_, 1,
                         MPI::DOUBLE, MPI::SUM);
  dt = oldPropensity/globalTotalPropensity_*(next_react_time_-current_time);
  next_react_time_ = current_time + dt;
  return next_react_time_;
}

double Compartment::getNewNextTime(ParallelEnvironment &pe,
                                   double current_time) {
  double dt(std::numeric_limits<double>::infinity());
  getNewPropensity(current_time);
  globalTotalPropensity_ = totalPropensity_;
  pe.getcart().Allreduce(MPI::IN_PLACE, &globalTotalPropensity_, 1,
                         MPI::DOUBLE, MPI::SUM);
  if( pe.getrank()==0 ) {
    dt = -log((*randdbl_)())/globalTotalPropensity_;
  }
  pe.getcart().Bcast(&dt, 1 , MPI_DOUBLE, 0);
  next_react_time_ = current_time + dt;
  return next_react_time_;
}


double Compartment::executeDirectMethodReaction(Lattice &g,
                                                ParallelEnvironment &pe,
                                                double current_time) {
  double propensities[pe.getsize()];
  for (unsigned i(0); i < pe.getsize(); ++i) {
    propensities[i] = 0;
  }
  propensities[pe.getrank()] = totalPropensity_;
  pe.getcart().Allreduce(MPI::IN_PLACE, &propensities, pe.getsize(),
                         MPI::DOUBLE, MPI::MAX);
  double random(0);
  if( pe.getrank()==0 ) {
    random = (*randdbl_)();
  }
  pe.getcart().Bcast(&random, 1 , MPI_DOUBLE, 0);
  double global(0);
  for (unsigned i(0); i < pe.getsize(); ++i) {
    global += propensities[i];
  }
  random *= global;
  double local(0);
  for (unsigned i(0); i < pe.getsize(); ++i) {
    local += propensities[i];
    if(local >= random) {
      if(i == pe.getrank()) {
        double randPropensity = 0.0;
        double accuPropensity = 0.0;
        randPropensity = totalPropensity_*(*randdbl_)();
        vector<Reaction*>::iterator p(
                                               independent_reactions_.begin());
        vector<Reaction*>::iterator e(independent_reactions_.end());
        while(p!=e) {
          accuPropensity += (*p)->getPropensity();  // accumulating propensity
          if(accuPropensity >= randPropensity) break;
          p++;
        }
        processingReaction( *(*p), g, pe );
      }
      break;
    }
  }
  next_react_time_ = current_time;
  return 0;
}

void Compartment::calculateProbability(Reaction &f, Lattice &g) {
  const int r0 = f.getR0();        // moles of reactant 0
  const int r1 = f.getR1();        // moles of reactant 1
  Species  *s0 = f.getS0();        // species of reactant 0
  Species  *s1 = f.getS1();        // species of reactant 1
  const double kAB = f.getK();
  const double DA = s0->getD();
  const double DB = s1->getD();
  const double rv = g.getradius();
  const double SQR2 = 1.414213562;
  if(r0==1 && r1==1) {
    const double probability = kAB /(6.0*SQR2*rv*(DA+DB));
    f.setProbability(probability);
  }
}

void Compartment::findMaxProbability(Species &s) {
  int m = influenced_reactions_.size(); 
  double maxPi = 0.0;
  for(int j=0; j<m; ++j) {
    Reaction& Rj(*influenced_reactions_[j]);
    Species* s0 = Rj.getS0();
    Species* s1 = Rj.getS1(); 
    if( s==*s0 || s==*s1 ) {
      const double Pi = Rj.getProbability();
      maxPi = ( maxPi<Pi ? Pi : maxPi );
    }
  }
  s.set_rho(maxPi);
}


//Type of vacant ids:
//vacant_id_
//ghost_id_
//+/- out_id_

void Compartment::walk(std::vector<Species*>& species_list,
                       Lattice &g, ParallelEnvironment &pe) {
  //check_voxels(g, 0);
  std::vector<unsigned> ghost_species_ids;
  std::vector<unsigned> ghost_src_coords;
  std::vector<unsigned> ghost_tar_coords;
  for (unsigned i(0); i < species_list.size(); ++i) {
    Species& s(*species_list[i]);
    std::vector<Molecule>& molecules(species_molecules_[s.getID()]);
    const unsigned begin_mol_size(molecules.size());
    unsigned index(0);
    std::vector<unsigned>& is_reactive(is_reactive_[s.getID()]);
    std::vector<double>& reaction_probabilities(reaction_probabilities_[
                                                s.getID()]);
    while (index < begin_mol_size && index < molecules.size()) { 
      Molecule& src_molecule(molecules[index]);
      const unsigned src_coord(src_molecule.coord);
      Voxel& src_voxel(g.getVoxel(src_coord)); 
      const unsigned adj_index((*adj_rand_)());
      const unsigned tar_coord(g.get_adjacent_coord(src_coord, adj_index));
      const int src_voxel_species_id(src_voxel.species_id);
      Voxel& tar_voxel(g.getVoxel(tar_coord));
      int tar_voxel_species_id(tar_voxel.species_id);
      if(tar_voxel_species_id == vacant_id_) { //vacant and not ghost voxel
        if (s.get_walk_probability() == 1 ||
            (*randdbl_)() < s.get_walk_probability()) {
          if (src_voxel_species_id < out_id_) { //src in normal voxel
            tar_voxel = src_voxel;
            src_voxel.species_id = vacant_id_;
          }
          else {
            //src leaving out_voxel to occupy normal vacant voxel:
            //populate a normal molecule in the tar voxel:
            const unsigned oindex(src_voxel.mol_index);
            const unsigned src_out_id(src_voxel_species_id/out_id_);
            tar_voxel.species_id = src_voxel_species_id%out_id_;
            tar_voxel.mol_index = outmolecules_[src_out_id][oindex].mol_index;
            vacate_outmolecule_in_voxel(g, src_out_id, src_voxel);
          }
          src_molecule.coord = tar_coord;
        }
      } //tar is a normal vacant out voxel:
      else if (tar_voxel_species_id >= out_id_ && 
               tar_voxel_species_id%out_id_ == 0) {
        if (s.get_walk_probability() == 1 ||
            (*randdbl_)() < s.get_walk_probability()) {
          //src is occupying a normal voxel:
          if (src_voxel_species_id < out_id_) { 
            const unsigned tar_out_id(tar_voxel_species_id/out_id_);
            populate_outmolecule_in_voxel(tar_out_id,
                                          src_voxel_species_id,
                                          src_voxel.mol_index, tar_coord,
                                          tar_voxel);
            src_voxel.species_id = vacant_id_;
          }
          else { //src is occupying an out_voxel and tar is a vacant out voxel:
            const unsigned src_out_id(src_voxel_species_id/out_id_);
            const unsigned tar_out_id(tar_voxel_species_id/out_id_);
            if (src_out_id == tar_out_id) {
              tar_voxel = src_voxel;
              //make src voxel vacant:
              src_voxel.species_id = tar_voxel_species_id;
              outmolecules_[src_out_id][src_voxel.mol_index].coord = tar_coord;
            }
            else {
              const int src_species_id(src_voxel_species_id%out_id_);
              const unsigned oindex(src_voxel.mol_index);
              populate_outmolecule_in_voxel(tar_out_id,
                                            src_species_id,
                      outmolecules_[src_out_id][oindex].mol_index, tar_coord,
                                            tar_voxel);
              vacate_outmolecule_in_voxel(g, src_out_id, src_voxel);
            }
          }
          src_molecule.coord = tar_coord;
        }
      }
      else if (is_parallel_ && tar_voxel_species_id < 0) { //ghost voxel
        ghost_species_ids.push_back(s.getID());
        ghost_src_coords.push_back(src_coord);
        ghost_tar_coords.push_back(tar_coord);
      }
      else {
        if (tar_voxel_species_id > out_id_) {
          tar_voxel_species_id = tar_voxel_species_id%out_id_;
        }
        if (is_reactive[tar_voxel_species_id] &&
            (reaction_probabilities[tar_voxel_species_id] == 1 ||
             (*randdbl_)() < reaction_probabilities[tar_voxel_species_id])) {
          const unsigned old_size(molecules.size());
          react(g, s, tar_voxel_species_id, src_voxel, tar_voxel, src_coord,
                tar_coord);
          if (molecules.size() < old_size) {
            --index;
          }
        }
      }
      ++index;
    }
  }
  if (is_parallel_) {
    walk_on_ghost(ghost_species_ids, ghost_src_coords, ghost_tar_coords, g, pe);
  }
}

//combine reaction_on_out() with react(()
void Compartment::react(Lattice& g, Species& s, const unsigned tarID,
                        Voxel& src_voxel, Voxel& tar_voxel,
                        const unsigned src_coord, const unsigned tar_coord) {
  const unsigned reaction_id(influenced_reaction_ids_[s.getID()][tarID]);
  Reaction& reaction(*influenced_reactions_[reaction_id]);
  Species& s0(*reaction.getS0());
  if (&s0 == &s) {
    do_reaction(g, reaction, src_voxel, tar_voxel, src_coord, tar_coord);
  }
  else {
    do_reaction(g, reaction, tar_voxel, src_voxel, tar_coord, src_coord);
  }
}

void Compartment::do_reaction(Lattice& g, Reaction& r, Voxel& v0, Voxel& v1,
                              const unsigned c0, const unsigned c1) {
  Species& s0(*r.getS0());
  Species& s1(*r.getS1());
  if (&s0 != &s1) {
    if (r.getP2()) {
      Species& p0(*r.getS2());
      if (&p0 != &s0 && &p0 != &s1) {
        //must remove before adding since in the removed species, the replacing
        //molecule might be the one removed and it will be occupied
        //by the new molecule:
        remove_molecule(g, s0.getID(), v0);
        add_molecule(p0.getID(), c0, get_new_mol_id(), v0);
      }
    }
    else {
      remove_molecule(g, s0.getID(), v0);
    }
    if (r.getP3()) {
      Species& p1(*r.getS3());
      if (&p1 != &s1 && &p1 != &s0) {
        //must remove before adding since in the removed species, the replacing
        //molecule might be the one removed and it will be occupied
        //by the new molecule:
        remove_molecule(g, s1.getID(), v1);
        add_molecule(p1.getID(), c1, get_new_mol_id(), v1);
      }
    }
    else {
      remove_molecule(g, s1.getID(), v1);
    }
  }
  else {
  }
}

void Compartment::add_molecule(const int species_id, const unsigned coord,
                               const unsigned mol_id, Voxel& voxel) {
  std::vector<Molecule>& molecules(species_molecules_[species_id]);
  const unsigned mol_index(molecules.size());
  molecules.push_back(Molecule(coord, mol_id));
  if (voxel.species_id >= out_id_) {//out voxel
    const unsigned out_id(voxel.species_id/out_id_);
    populate_outmolecule_in_voxel(out_id, species_id, mol_index, coord, voxel);
  }
  else {
    voxel.species_id = species_id;
    voxel.mol_index = mol_index;
  }
}

unsigned Compartment::remove_molecule(Lattice& g, const int species_id,
                                      Voxel& voxel) {
  int mol_index(voxel.mol_index);
  if (voxel.species_id > out_id_) { //out voxel
    const unsigned out_id(voxel.species_id/out_id_);
    mol_index = outmolecules_[out_id][mol_index].mol_index;
    vacate_outmolecule_in_voxel(g, out_id, voxel); 
  }
  else {
    voxel.species_id = vacant_id_;
  }
  return remove_molecule_from_molecules(g, species_id, mol_index);
}

void Compartment::populate_outmolecule_in_voxel(const int out_id,
                                                const int species_id,
                                                const unsigned mol_index,
                                                const unsigned coord,
                                                Voxel& voxel) {
  std::vector<OutMolecule>& outmolecules(outmolecules_[out_id]);
  voxel.species_id = out_id*out_id_+species_id;
  voxel.mol_index = outmolecules.size();
  outmolecules.push_back(OutMolecule(coord, mol_index, species_id));
}

void Compartment::vacate_outmolecule_in_voxel(Lattice& g, const int out_id,
                                              Voxel& voxel) {
  std::vector<OutMolecule>& outmolecules(outmolecules_[out_id]);
  const unsigned oindex(voxel.mol_index);
  OutMolecule& back_outmolecule(outmolecules.back());
  Voxel& back_voxel(g.getVoxel(back_outmolecule.coord));
  back_voxel.mol_index = oindex;
  outmolecules[oindex] = back_outmolecule;
  outmolecules.pop_back();
  voxel.species_id = out_id*out_id_;
}

unsigned Compartment::remove_molecule_from_molecules(Lattice& g,
                                                     const int species_id,
                                                     const unsigned mol_index) {
  //remove molecule from molecule list:
  //Since we will be updating the molecule list (molecules), if the back 
  //molecule in the list is an outmolecule, we need to update its mol_index
  //to point to the new index in molecules. The back molecule's 
  //voxel.mol_index is unchanged because outmolecules_ is unchanged.
  std::vector<Molecule>& molecules(species_molecules_[species_id]);
  const Molecule& back_molecule(molecules.back());
  const unsigned back_coord(back_molecule.coord);
  Voxel& back_voxel(g.getVoxel(back_coord));
  //we could have removed this back outmolecule before calling this function
  //so need ensure that the voxel is still not vacant before removing the
  //outmolecule:
  if (back_voxel.species_id > out_id_ && 
      back_voxel.species_id%out_id_ != 0) {
    const unsigned out_id(back_voxel.species_id/out_id_);
    outmolecules_[out_id][back_voxel.mol_index].mol_index = mol_index; 
  }
  else {
    back_voxel.mol_index = mol_index;
  }
  const unsigned mol_id = molecules[mol_index].mol_id;
  molecules[mol_index] = back_molecule;
  molecules.pop_back();
  return mol_id;
}


void Compartment::calculateCollisionTime(Species& s) {
  s.calcCollisionTime();
  const unsigned id0(s.getID());
  for (unsigned i(0); i < is_reactive_[id0].size(); ++i) {
    if (is_reactive_[id0][i]) {
      const unsigned reaction_id(influenced_reaction_ids_[id0][i]);
      Reaction& reaction(*influenced_reactions_[reaction_id]);
      double probability(reaction.getProbability()*s.get_walk_probability());
      reaction_probabilities_[id0][i] = probability;
      std::cout << "reaction:" << reaction.getName() << std::endl;
      std::cout << "\tspecies:" << s.getName() << std::endl;
      std::cout << "\tp:" << probability << std::endl;
      std::cout << "\twalk probability:" << s.get_walk_probability() <<
        std::endl;
      std::cout << "\twalk interval:" << s.getDt() <<
        std::endl;
    }
  }
}

void Compartment::attachInfluencedReaction(Reaction& r) {
  Species& s0(*r.getS0());
  Species& s1(*r.getS1());
  const unsigned reaction_id(influenced_reactions_.size());
  influenced_reactions_.push_back(&r);
  const unsigned id0(s0.getID());
  const unsigned id1(s1.getID());
  if (!influenced_reaction_ids_.size()) {
    const unsigned species_size(species_list_.size());
    influenced_reaction_ids_.resize(species_size);
    is_reactive_.resize(species_size);
    reaction_probabilities_.resize(species_size);
    for (unsigned i(0); i < species_size; ++i) {
      influenced_reaction_ids_[i].resize(species_size);
      is_reactive_[i].resize(species_size, 0);
      reaction_probabilities_[i].resize(species_size, 0);
    }
  }
  influenced_reaction_ids_[id0][id1] = reaction_id;
  influenced_reaction_ids_[id1][id0] = reaction_id;
  is_reactive_[id0][id1] = 1;
  is_reactive_[id1][id0] = 1;
}

unsigned Compartment::check_species_size(Lattice& g, ParallelEnvironment& pe,
                                         const unsigned sid,
                                         const double id) {
  Species& s(*species_list_[sid]);
  const std::vector<Molecule>& molecules(species_molecules_[s.getID()]);
  unsigned size(molecules.size());
  pe.getcart().Allreduce(MPI::IN_PLACE, &size, 1, MPI::INT, MPI::SUM);
  fout << id << " " << s.getName() << " size:" << size << std::endl;
  return size;
}


void Compartment::walk_on_ghost(std::vector<unsigned>& species_ids,
                                std::vector<unsigned>& src_coords,
                                std::vector<unsigned>& tar_coords, Lattice& g,
                                ParallelEnvironment& pe) { 
  std::vector<int> order = {0, 1, 2, 3, 4, 5, 6, 7};
  if (!pe.getrank()) {
    std::shuffle(order.begin(), order.end(), rng_);
  }
  pe.getcart().Bcast(order.data(), order.size(), MPI::INT, 0); 
  
  //prepare sub-vector of molecules for each sub-volume 
  std::vector<unsigned> sub_indices[8];
  for (unsigned i(0); i < 8; ++i) {
    sub_indices[i].reserve(src_coords.size()/4);
  }

  for (unsigned i(0); i < src_coords.size(); ++i) {
    // get linear coordinate of voxel which molecule resides
     const int coord(src_coords[i]);
     // identify sub-volume, GHOST_SIZE = 1
     const int xbit(g.getXind(coord) > (mid_span_.x + GHOST_SIZE - 1) ? 1 : 0);
     const int ybit(g.getYind(coord) > (mid_span_.y + GHOST_SIZE - 1) ? 1 : 0);
     const int zbit(g.getZind(coord) > (mid_span_.z + GHOST_SIZE - 1) ? 1 : 0);
     const int index(xbit*4 + ybit*2 + zbit);
     // append molecule index to sub-vector
     sub_indices[index].push_back(i);
  } 
  
  for(unsigned sv(0); sv < 8; ++sv) {
    g.jumpincoords.clear(); 
    std::vector<SpillMolecule> spill_coords;
    const unsigned N(order[sv]);
    // load ghost cells from adjacent process using MPI 
    // this is the costliest operation, need to further optimize
    // by keeping a list of molecules that are in the in-compartment-ghost:
    g.load_ghost(pe, outmolecules_[N+1], outmolecules_[9], N);
    // loop over sub-moleculeVector
    for (unsigned i(0); i < sub_indices[N].size(); ++i) {
      const unsigned index(sub_indices[N][i]);
      const unsigned tar_coord(tar_coords[index]);
      const unsigned src_coord(src_coords[index]);
      Voxel& tar_voxel(g.getVoxel(tar_coord));
      int tar_voxel_species_id(tar_voxel.species_id);
      Voxel& src_voxel(g.getVoxel(src_coord));
      const int src_voxel_species_id(src_voxel.species_id);
      const int src_species_id(species_ids[index]);
      //some molecules that are supposed to move into ghost may have reacted
      //so we need to ensure the id is consistent:
      Species& src(*species_list_[src_species_id]);
      if (src_species_id == src_voxel_species_id ||
          src_species_id == src_voxel_species_id%out_id_) {
        //vacant voxel of ghost or out ghost voxel, tar_voxel_species_id < 0:
        if (tar_voxel_species_id == ghost_id_ || 
            tar_voxel_species_id%out_id_ == 0) {
          if (src.get_walk_probability() == 1 ||
              (*randdbl_)() < src.get_walk_probability()) {
            if (tar_voxel_species_id == ghost_id_) {
              tar_voxel.species_id = -src_species_id;
            }
            else { //tar is a vacant out_ghost_voxel
              tar_voxel.species_id = tar_voxel_species_id-src_species_id;
            }
            const unsigned mol_id(remove_molecule(g, src_species_id,
                                                  src_voxel));
            spill_coords.push_back(SpillMolecule(src_species_id, tar_coord,
                                                mol_id));
          }
        }
        else {
          if (tar_voxel_species_id < -out_id_) {
            tar_voxel_species_id = tar_voxel_species_id%out_id_;
          }
          tar_voxel_species_id = -tar_voxel_species_id;
          if (is_reactive_[src_species_id][tar_voxel_species_id] && 
              (reaction_probabilities_[src_species_id][tar_voxel_species_id] ==
               1 || (*randdbl_)() <
               reaction_probabilities_[src_species_id][tar_voxel_species_id])) {
            react_on_ghost(g, src, tar_voxel_species_id, src_voxel, tar_voxel,
                           src_coord, tar_coord, spill_coords);
          }
        }
      }
    }
    //changed state of voxels in ghost (molecule added or removed,
    //or changed species):
    //There could be duplicates in the spill_coords but it is ok.
    //It is only an indicator that the voxel in the local compartment
    //should be updated according to the species_id when reading jumpin 
    //molecules:
    spill_molecules(g, spill_coords, ghosts_[N]);
    g.clear_ghost(spill_coords);
    g.feedGhost(begin_[N].x, end_[N].x, begin_[N].y, end_[N].y,
                begin_[N].z, end_[N].z, bit_[N].x, bit_[N].y, bit_[N].z,
                pe);
    jumpin_molecules(g, pe);
  }
}

void Compartment::react_on_ghost(Lattice& g, Species& s, const unsigned tarID,
                                 Voxel& src_voxel, Voxel& tar_voxel,
                                 const unsigned src_coord,
                                 const unsigned tar_coord, 
                                 std::vector<SpillMolecule>& spill_coords) {
  const unsigned reaction_id(influenced_reaction_ids_[s.getID()][tarID]);
  Reaction& reaction(*influenced_reactions_[reaction_id]);
  Species& s0(*reaction.getS0());
  if (&s0 == &s) {
    do_reaction_on_ghost(g, reaction, src_voxel, tar_voxel, src_coord,
                         tar_coord, spill_coords);
  }
  else {
    do_reaction_on_ghost(g, reaction, tar_voxel, src_voxel, tar_coord,
                         src_coord, spill_coords);
  }
}


void Compartment::do_reaction_on_ghost(Lattice& g, Reaction& r, Voxel& v0,
                                       Voxel& v1, const unsigned c0,
                                       const unsigned c1,
                                     std::vector<SpillMolecule>& spill_coords) {
  Species& s0(*r.getS0());
  Species& s1(*r.getS1());
  if (&s0 != &s1) {
    if (r.getP2()) {
      Species& p0(*r.getS2());
      if (&p0 != &s0 && &p0 != &s1) {
        //must remove before adding since in the removed species, the replacing
        //molecule might be the one removed and it will be occupied
        //by the new molecule:
        remove_molecule_on_ghost(g, s0.getID(), c0, v0, spill_coords);
        add_molecule_on_ghost(p0.getID(), c0, v0, spill_coords);
      }
    }
    else {
      remove_molecule_on_ghost(g, s0.getID(), c0, v0, spill_coords);
    }
    if (r.getP3()) {
      Species& p1(*r.getS3());
      if (&p1 != &s1 && &p1 != &s0) {
        //must remove before adding since in the removed species, the replacing
        //molecule might be the one removed and it will be occupied
        //by the new molecule:
        remove_molecule_on_ghost(g, s1.getID(), c1, v1, spill_coords);
        add_molecule_on_ghost(p1.getID(), c1, v1, spill_coords);
      }
    }
    else {
      remove_molecule_on_ghost(g, s1.getID(), c1, v1, spill_coords);
    }
  }
  else {
  }
}

void Compartment::add_molecule_on_ghost(const int species_id,
                                        const unsigned coord, Voxel& v,
                                    std::vector<SpillMolecule>& spill_coords) {
  if (v.species_id < 0) { //ghost voxel
    //if out_ghost_voxel,convert to out ghost species:
    if (v.species_id <= -out_id_) {
      v.species_id = (v.species_id/out_id_)*out_id_-species_id;
    }
    else { //convert to ghost species:
      v.species_id = -species_id;
    }
    spill_coords.push_back(SpillMolecule(species_id, coord, get_new_mol_id()));
  }
  else {
    add_molecule(species_id, coord, get_new_mol_id(), v);
  }
}

void Compartment::remove_molecule_on_ghost(Lattice& g, const int species_id,
                                           const unsigned coord, Voxel& v,
                                     std::vector<SpillMolecule>& spill_coords) {
  if (v.species_id < 0) { //ghost voxel
    if (v.species_id < -out_id_) { // out_ghost_voxel
      v.species_id = (v.species_id/out_id_)*out_id_;
    }
    else {
      v.species_id = ghost_id_;
    }
    spill_coords.push_back(SpillMolecule(vacant_id_, coord, 0));
  }
  else {
    remove_molecule(g, species_id, v);
  }
}

void Compartment::spill_molecules(Lattice& g,
                                const std::vector<SpillMolecule>& spill_coords,
                                const Vector<unsigned>& ghost) {
  // clear spillout molecules for every sub-volume
  g.spillcoordsX.clear();
  g.spillcoordsY.clear();
  g.spillcoordsZ.clear();
  for (unsigned i(0); i < spill_coords.size(); ++i) {
    const SpillMolecule& species_coord(spill_coords[i]);
    const unsigned coord(species_coord.coord);
    const int xadj(g.getXind(coord));
    const int yadj(g.getYind(coord));
    const int zadj(g.getZind(coord)); 
    if (zadj == ghost.z) {
      g.spillcoordsZ.push_back(species_coord);
    }
    else if (yadj == ghost.y) {
      g.spillcoordsY.push_back(species_coord);
    }
    else if (xadj == ghost.x) {
      g.spillcoordsX.push_back(species_coord);
    }
  }
}


//not only jumpin, but also changed species, or deleted
//jump only into out_voxels, never into out_ghost_voxels or into normal voxels:
void Compartment::jumpin_molecules(Lattice& g, ParallelEnvironment& pe) {
  int cnt(0);
  const unsigned size(g.jumpincoords.size());
  for (unsigned i(0); i < size; ++i) {
    const SpillMolecule& sc(g.jumpincoords[i]);
    int species_id(sc.species_id);
    const unsigned coord(sc.coord);
    Voxel& v(g.getVoxel(coord));
    const int old_species_id(v.species_id);
    //if voxel state changed:
    if (old_species_id%out_id_ != species_id) {
      if (old_species_id%out_id_ == 0) { //if previously vacant out voxel:
        add_molecule(species_id, coord, sc.mol_id, v);
      }
      else { //if species has changed:
        remove_molecule(g, old_species_id%out_id_, v); 
        if (species_id != vacant_id_) { //replace with new species:
          add_molecule(species_id, coord, sc.mol_id, v);
        }
      }
    }
  }
}

void Compartment::throwinOneMolecule(Species &s, Lattice &g) {
  const int i0 = g.getXdim();
  const int j0 = g.getYdim();
  const int k0 = g.getZdim();
  const int index = g.linearCoord(i0/2,j0/2,k0/2);    // central point
  const int specID(s.getID());
  const double x = g.getXpos(i0/2, j0/2, k0/2);
  const double y = g.getYpos(i0/2, j0/2, k0/2);
  const double z = g.getZpos(i0/2, j0/2, k0/2);
  species_molecules_[specID].push_back(Molecule(index, get_new_mol_id()));
  g.getVoxel(index).species_id = specID;
  g.getVoxel(index).mol_index = species_molecules_[specID].size()-1;
}


void Compartment::addCoordinatesSpecies(Species& species) {
  output_coord_species_.push_back(species);
}

void Compartment::addNumberSpecies(Species& species) {
  output_number_species_.push_back(species);
}

void Compartment::outputCoordinatesHeader(Lattice& lattice,
                                          ParallelEnvironment &pe) {
  fout2 << "Time";
  for (unsigned i(0); i < output_coord_species_.size(); ++i) {
    Species& s(output_coord_species_[i]);
    fout2 << "," << s.getName();
  }
  fout2 << "," << lattice.getradius() << " " << pe.getsize() << " " << 
    pe.getNx() << " " << pe.getNy() << " " << pe.getNz() << endl;
}

void Compartment::outputCoordinates(const double current_time) {
  fout2 << setprecision(15) << current_time;
  for (unsigned i(0); i < output_coord_species_.size(); ++i) {
    Species& s(output_coord_species_[i]);
    std::vector<Molecule>& molecules(species_molecules_[s.getID()]);
    vector<Molecule>::iterator p2 = molecules.begin();
    vector<Molecule>::iterator e2 = molecules.end();
    fout2 << "," ;
    while (p2!=e2) {
      // get linear coordinate of voxel
      if (p2 == molecules.begin()) {
        fout2 << p2->coord << " " << p2->mol_id;
      }
      else {
        fout2 << " " << p2->coord << " " << p2->mol_id;
      }
      p2++;
    }
  }
  fout2 << endl;
}


void Compartment::outputNumbersHeader() {
  fout3 << setprecision(15) << "#time(sec)";
  for (unsigned i(0); i < output_number_species_.size(); ++i) {
    Species& s(output_number_species_[i]);
    fout3 << "," << s.getName();
  }
  fout3 << endl;
}


void Compartment::outputNumbers(const double current_time) {
  fout3 << setprecision(15) << current_time;
  for (unsigned i(0); i < output_number_species_.size(); ++i) {
    Species& s(output_number_species_[i]);
    fout3 << "," << species_molecules_[s.getID()].size();
  }
  fout3 << endl;
}
