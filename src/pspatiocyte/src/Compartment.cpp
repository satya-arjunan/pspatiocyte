//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of pSpatiocyte
//
//        Copyright (C) 2019 Satya N.V. Arjunan
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// Motocyte is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// Motocyte is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with Motocyte -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Satya Arjunan <satya.arjunan@gmail.com>
// and Atsushi Miyauchi
//


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
  nx_ = pe.getispan(); //x size allocated to local process from global lattice 
  ny_ = pe.getjspan(); //y size allocated to local process from global lattice 
  nz_ = pe.getkspan(); //z size allocated to local process from global lattice 
  //LNx = nx_+2 (2 for ghost voxels at the sides)
  //LNy = ny_+2 (2 for ghost voxels at the sides)
  //LNz = nz_+2 (2 for ghost voxels at the sides)
  //local lattice size is (nx_+2)*(ny_+2)*(nz_*2) or LNx*LNy*LNz
  //molecules can only occupy coordinates:
  //(x=1 to nx_),(y=1 to ny_),(z=1 to nz)
  fout << "max:" << nx_ << " " << ny_ << " " << nz_ << std::endl;
  // x=0 and x=nx_+1 are occupied by ghost or invalid voxels
  // y=0 and y=ny_+1 are occupied by ghost or invalid voxels
  // z=0 and z=nz_+1 are occupied by ghost or invalid voxels
  const int NE = MPI::PROC_NULL;  // no neighbor
  const double rv = g.getradius();
  const double SQR2 = 1.414213562;

  //set all invalid ids:
  if (type_==VOLUME) {
    // if right-boundary node
    if(pe.getneighbor_i_plus()==NE)
    for(int j=0; j<=ny_+1; ++j)
    for(int k=0; k<=nz_+1; ++k)
    {
        const int lc0 = g.linearCoord(nx_  , j, k);
        const int lc1 = g.linearCoord(nx_+1, j, k);
        g.get_voxel(lc0).species_id =
        g.get_voxel(lc1).species_id = invalid_id_;
    }

    // if left-boundary node
    if(pe.getneighbor_i_minus()==NE)
    for(int j=0; j<=ny_+1; ++j)
    for(int k=0; k<=nz_+1; ++k)
    {
        const int lc0 = g.linearCoord(0, j, k);
        const int lc1 = g.linearCoord(1, j, k);
        g.get_voxel(lc0).species_id =
        g.get_voxel(lc1).species_id = invalid_id_;
    }

    // if rear-boundary node
    if(pe.getneighbor_j_plus()==NE)
    for(int i=0; i<=nx_+1; ++i)
    for(int k=0; k<=nz_+1; ++k)
    {
        const int lc0 = g.linearCoord(i, ny_  , k);
        const int lc1 = g.linearCoord(i, ny_+1, k);
        g.get_voxel(lc0).species_id =
        g.get_voxel(lc1).species_id = invalid_id_;
    }

    // if front-boundary node
    if(pe.getneighbor_j_minus()==NE)
    for(int i=0; i<=nx_+1; ++i)
    for(int k=0; k<=nz_+1; ++k)
    {
        const int lc0 = g.linearCoord(i, 0, k);
        const int lc1 = g.linearCoord(i, 1, k);
        g.get_voxel(lc0).species_id =
        g.get_voxel(lc1).species_id = invalid_id_;
    }

    // if top boundary node
    if(pe.getneighbor_k_plus()==NE)
    for(int i=0; i<=nx_+1; ++i)
    for(int j=0; j<=ny_+1; ++j)
    {
        const int lc0 = g.linearCoord(i, j, nz_  );
        const int lc1 = g.linearCoord(i, j, nz_+1);
        g.get_voxel(lc0).species_id =
        g.get_voxel(lc1).species_id = invalid_id_;
    }

    // if bottom boundary node
    if(pe.getneighbor_k_minus()==NE)
    for(int i=0; i<=nx_+1; ++i)
    for(int j=0; j<=ny_+1; ++j)
    {
        const int lc0 = g.linearCoord(i, j, 0);
        const int lc1 = g.linearCoord(i, j, 1);
        g.get_voxel(lc0).species_id =
        g.get_voxel(lc1).species_id = invalid_id_;
    }

    // append inner voxels
    for(int i=1; i<=nx_; ++i)
    for(int j=1; j<=ny_; ++j)
    for(int k=1; k<=nz_; ++k)
    {
        const int lc = g.linearCoord(i, j, k);
        if(g.get_voxel(lc).species_id == vacant_id_)
            voxel_coords_.push_back(lc);
    }
    n_voxels_ = voxel_coords_.size(); // = (nx_-2)*(ny_-2)*(nz_-2)
    unsigned total_voxels(n_voxels_);
    local_volume_ = n_voxels_*4.0*SQR2*rv*rv*rv;
    global_volume_ = local_volume_;
    const MPI::Cartcomm cart = pe.getcart();
    cart.Allreduce(MPI::IN_PLACE, &global_volume_, 1, MPI::DOUBLE, MPI::SUM);
    cart.Allreduce(MPI::IN_PLACE, &total_voxels, 1, MPI::INT, MPI::SUM);
    if (!pe.getrank()) {
      std::cout << "voxel radius:" << rv << std::endl;
      std::cout << "total voxels:" << total_voxels << std::endl;
      std::cout << "actual total volume:" << global_volume_ << std::endl;
    }
  }

  mid_span_.x = nx_/2 + nx_%2;
  mid_span_.y = ny_/2 + ny_%2;
  mid_span_.z = nz_/2 + nz_%2; 
  for (unsigned N(0); N < 8; ++N) {
    bit_[N].x = ((N/4)%2);
    bit_[N].y = ((N/2)%2);
    bit_[N].z = ((N/1)%2);
    begin_[N].x = (bit_[N].x*mid_span_.x + GHOST_SIZE);
    begin_[N].y = (bit_[N].y*mid_span_.y + GHOST_SIZE);
    begin_[N].z = (bit_[N].z*mid_span_.z + GHOST_SIZE);
    end_[N].x = (bit_[N].x*(nx_-mid_span_.x)+mid_span_.x + GHOST_SIZE - 1);
    end_[N].y = (bit_[N].y*(ny_-mid_span_.y)+mid_span_.y + GHOST_SIZE - 1);
    end_[N].z = (bit_[N].z*(nz_-mid_span_.z)+mid_span_.z + GHOST_SIZE - 1);
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
  for(int i=0; i<=nx_+1; ++i) {
    for(int j=0; j<=ny_+1; ++j) {
      for(int k=0; k<=nz_+1; ++k) {
        const unsigned coord(g.linearCoord(i, j, k));
        Voxel& voxel(g.get_voxel(coord));
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
  for(int i=0; i<=nx_+1; ++i) {
    for(int j=0; j<=ny_+1; ++j) {
      for(int k=0; k<=nz_+1; ++k) {
        const unsigned coord(g.linearCoord(i, j, k));
        Voxel& voxel(g.get_voxel(coord));
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
              species_list_[species_id]->get_name() << std::endl;
            abort();
          }
          if (mols[voxel.mol_index].coord != coord) {
            std::cout << id << 
              " error not consistent species id on normal voxel" <<
              " coord:" << coord << " other coord:" << 
              mols[voxel.mol_index].coord << " species id:" << species_id <<
              " name:" << species_list_[species_id]->get_name() <<
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
  for(int i=0; i<=nx_+1; ++i) {
    for(int j=0; j<=ny_+1; ++j) {
      for(int k=0; k<=nz_+1; ++k) {
        const unsigned coord(g.linearCoord(i, j, k));
        Voxel& voxel(g.get_voxel(coord));
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

unsigned Compartment::coord_to_subvolume(Lattice& g, const unsigned coord) {
  int i, j, k;
  g.coord_to_ijk(coord, i, j, k);
  const int xbit(i > (mid_span_.x + GHOST_SIZE - 1) ? 1 : 0);
  const int ybit(j > (mid_span_.y + GHOST_SIZE - 1) ? 1 : 0);
  const int zbit(k > (mid_span_.z + GHOST_SIZE - 1) ? 1 : 0);
  return xbit*4 + ybit*2 + zbit;
}

void Compartment::populate_molecules(Species& s, unsigned size,
                                     Lattice &g, ParallelEnvironment &pe) {
  double rv = g.getradius();
  if (type_ == VOLUME) {
    s.setVolumeDt(rv);
  }

  const Vector<float>& origin(s.get_populate_origin());
  const Vector<float>& range(s.get_populate_range());
  if (!size && origin.x == 0.5 && origin.y == 0.5 && origin.z == 0.5 &&
      range.x == 1 && range.y == 1 && range.z == 1) { 
    return;
  }

  if (size > n_voxels_) {
    std::cout << "Too many molecules to be populated, size:" << size << " "
      << " available:" << n_voxels_ << " " << s.get_name() << std::endl;
    abort();
  }

  if (origin.x == 0.5 && origin.y == 0.5 && origin.z == 0.5 &&
      range.x == 1 && range.y == 1 && range.z == 1) { 
    populate_molecules(s, size, g);
  }
  //log ghost, vacant and out voxels of each subvolume in the local subdomain:
  else if (range.x == -1) {
    std::string rank(std::to_string(pe.getrank()));;
    for(int i=0; i<=nx_+1; ++i) {
      for(int j=0; j<=ny_+1; ++j) {
        for(int k=0; k<=nz_+1; ++k) {
          const unsigned coord(g.linearCoord(i, j, k));
          Voxel& voxel(g.get_voxel(coord));
          int species_id(voxel.species_id); 
          std::string sv(std::to_string(coord_to_subvolume(g, coord)));
          if (species_id >= out_id_ &&
              s.get_name() == std::string("Out")+rank+sv) {
            add_molecule(s.get_id(), coord, get_new_mol_id(), voxel); 
          }
          else if (species_id <= -out_id_ &&
                   s.get_name() == std::string("OutGhost")+rank+sv) { 
            add_molecule(s.get_id(), coord, get_new_mol_id(), voxel); 
          }
          else if (species_id == ghost_id_ &&
                   s.get_name() == std::string("Ghost")+rank+sv) { 
            add_molecule(s.get_id(), coord, get_new_mol_id(), voxel); 
          }
          else if (species_id == vacant_id_ &&
                   s.get_name() == std::string("Vacant")+rank+sv) { 
            add_molecule(s.get_id(), coord, get_new_mol_id(), voxel); 
          }
        }
      }
    }
  }
  else { 
    pe.getcart().Allreduce(MPI::IN_PLACE, &size, 1, MPI::INT, MPI::SUM);
    if (pe.getrank() == 0) {
      const int species_id(s.get_id());
      //int coord(g.linearCoordFast(nx_-20-(*adj_rand_)(),
      //ny_-20-(*adj_rand_)(), nz_-20-(*adj_rand_)()));
      int coord(g.linearCoordFast(nx_, ny_, nz_));
      int sid(g.get_voxel(coord).species_id);
      if (!((sid < out_id_ && sid != vacant_id_) ||
            (sid > out_id_ && sid%out_id_ != 0))) {
        add_molecule(species_id, coord, get_new_mol_id(), g.get_voxel(coord)); 
      }
      else {
        std::cout << "Unable to populate ranged" << s.get_name() << std::endl;
        abort();
      }
    }
  }
}


void Compartment::populate_molecules(Species& s, const unsigned size,
                                     Lattice &g) {
  const int species_id(s.get_id());
  for (unsigned i(0); i < size; ++i) {
    unsigned trial = 0;
    unsigned index;
    unsigned coord;
    int sid;
    do {
      if(trial++ > n_voxels_) {
        fout << "ABORT: unable to find a vacant voxel to populate! ("
          << s.get_name() << ")" << endl;
        abort();
      }
      index = (int)(n_voxels_*(*randdbl_)());
      coord = voxel_coords_[index];
      sid = g.get_voxel(coord).species_id;
    } while((sid < out_id_ && sid != vacant_id_) ||
            (sid > out_id_ && sid%out_id_ != 0));
    add_molecule(species_id, coord, get_new_mol_id(), g.get_voxel(coord)); 
  }
}

void Compartment::add_direct_method_reaction(Reaction& r) {
  direct_method_reactions_.push_back(&r);
}

bool Compartment::do_direct_method_reaction(Reaction &f, Lattice &g,
                                                 ParallelEnvironment &pe) {
  const int r0 = f.getR0();       // moles of reactant 0
  const int r1 = f.getR1();       // moles of reactant 1
  const int r2 = f.getP2();       // moles of product  2
  const int r3 = f.getP3();       // moles of product  3
  Species& s0 = f.getS0();        // species of reactant 0
  Species& s1 = f.getS1();        // species of reactant 1
  Species& s2 = f.getS2();        // species of product  2
  Species& s3 = f.getS3();        // species of product  3

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
    result = remove_reactant(s0, g, binary, center, neighbor);
  }

  if (r1) {
    result = remove_reactant(s1, g, binary, center, neighbor);
  }

  if (result) {
    if (binary) { 
      const double location = (*randdbl_)();
      if(location>0.5) { 
        add_product(s2, g, center);
        add_product(s3, g, neighbor);
      }
      else {
        add_product(s2, g, neighbor);
        add_product(s3, g, center);
      }
    }
    else {
      if (r2) {
        add_product(s2, g, center);
      }
      else
      if (r3) {
        add_product(s3, g, center);
      }
    }
  }
  return result;
}

bool Compartment::remove_reactant(Species &s, Lattice &g, bool binary,
                                    int &center, int &neighbor) {
  std::vector<Molecule>& molecules(species_molecules_[s.get_id()]); 
  const int size(molecules.size());
  const int target((size==1) ? 0 : (int)(size*(*randdbl_)()));
  Molecule& molecule(molecules[target]);
  const int coord(molecule.coord);
  if(binary) { //if two products, find a vacant neighbor voxel:
    if (!get_vacant_neighbor(g, coord, center, neighbor)) {
      return false;
    }
  }
  else {
    center = coord;
  }
  Voxel& v0(g.get_voxel(coord));
  remove_molecule(g, s.get_id(), v0);
  return true;
}


bool Compartment::get_vacant_neighbor(Lattice& g, const int& coord, int& center,
                                      int& neighbor) {
  unsigned adj_index((*adj_rand_)());
  for (unsigned i(0); i < 12; ++i) {
    const int tar_coord(g.get_adjacent_coord(coord, adj_index));
    const int species_id(g.get_voxel(tar_coord).species_id);
    if (species_id == vacant_id_ || (species_id >= out_id_ &&
                                     species_id%out_id_ == 0)) {
      center   = coord;
      neighbor = tar_coord;
      return true;
    }
    ++adj_index;
    if (adj_index > 11) {
      adj_index = 0;
    }
  }
  if (is_force_search_vacant_) {
    for (unsigned i(0); i < 12; ++i) {
      const int tar_coord(g.get_adjacent_coord(coord, adj_index));
      if (walk_molecule(g, tar_coord, coord)) {
        center   = coord;
        neighbor = tar_coord;
        return true;
      }
      ++adj_index;
      if (adj_index > 11) {
        adj_index = 0;
      }
    }
  }
  return false;
}

bool Compartment::walk_molecule(Lattice& g, const int src_coord,
                                const int curr_coord) {
  Voxel& src_voxel(g.get_voxel(src_coord)); 
  const int src_voxel_species_id(src_voxel.species_id);
  if (src_voxel_species_id < 0 || src_voxel_species_id == invalid_id_) {
    return false;
  }
  const int species_id(src_voxel_species_id%out_id_);
  std::vector<Molecule>& molecules(species_molecules_[species_id]);

  unsigned adj_index((*adj_rand_)());
  for (unsigned i(0); i < 12; ++i) {
    const unsigned tar_coord(g.get_adjacent_coord(src_coord, adj_index));
    if (tar_coord != curr_coord) {
      Voxel& tar_voxel(g.get_voxel(tar_coord));
      int tar_voxel_species_id(tar_voxel.species_id);
      if(tar_voxel_species_id == vacant_id_) { //vacant and not ghost voxel
        if (src_voxel_species_id < out_id_) { //src in normal voxel
          Molecule& src_molecule(molecules[src_voxel.mol_index]);
          tar_voxel = src_voxel;
          src_voxel.species_id = vacant_id_;
          src_molecule.coord = tar_coord;
        }
        else {
          //src leaving out_voxel to occupy normal vacant voxel:
          //populate a normal molecule in the tar voxel:
          const unsigned oindex(src_voxel.mol_index);
          const unsigned src_out_id(src_voxel_species_id/out_id_);
          tar_voxel.species_id = src_voxel_species_id%out_id_;
          const unsigned mol_index(outmolecules_[src_out_id][oindex].mol_index);
          tar_voxel.mol_index = mol_index;
          vacate_outmolecule_in_voxel(g, src_out_id, src_voxel);
          Molecule& src_molecule(molecules[mol_index]);
          src_molecule.coord = tar_coord;
          check_voxels(g, 1);
        }
        return true;
      } //tar is a normal vacant out voxel:
      else if (tar_voxel_species_id >= out_id_ && 
               tar_voxel_species_id%out_id_ == 0) {
          //src is occupying a normal voxel:
        if (src_voxel_species_id < out_id_) { 
          Molecule& src_molecule(molecules[src_voxel.mol_index]);
          const unsigned tar_out_id(tar_voxel_species_id/out_id_);
          populate_outmolecule_in_voxel(tar_out_id,
                                        src_voxel_species_id,
                                        src_voxel.mol_index, tar_coord,
                                        tar_voxel);
          src_voxel.species_id = vacant_id_;
          src_molecule.coord = tar_coord;
        }
        else { //src is occupying an out_voxel and tar is a vacant out voxel:
          const unsigned src_out_id(src_voxel_species_id/out_id_);
          const unsigned tar_out_id(tar_voxel_species_id/out_id_);
          const unsigned oindex(src_voxel.mol_index);
          const unsigned mol_index(outmolecules_[src_out_id][oindex].mol_index);
          Molecule& src_molecule(molecules[mol_index]);
          if (src_out_id == tar_out_id) {
            tar_voxel = src_voxel;
            //make src voxel vacant:
            src_voxel.species_id = tar_voxel_species_id;
            outmolecules_[src_out_id][src_voxel.mol_index].coord = tar_coord;
          }
          else {
            const int src_species_id(src_voxel_species_id%out_id_);
            populate_outmolecule_in_voxel(tar_out_id,
                                          src_species_id,
                                          mol_index, tar_coord,
                                          tar_voxel);
            vacate_outmolecule_in_voxel(g, src_out_id, src_voxel);
          }
          src_molecule.coord = tar_coord;
        }
        return true;
      }
    }
    ++adj_index;
    if (adj_index > 11) {
      adj_index = 0;
    }
  }
  return false;
}


void Compartment::add_product(Species &s, Lattice &g, int coord) {
  Voxel& voxel(g.get_voxel(coord));
  add_molecule(s.get_id(), coord, get_new_mol_id(), voxel); 
}

double Compartment::get_reaction_propensity(Reaction &reaction) {
  const unsigned size(species_molecules_[reaction.getS0().get_id()].size()); 
  return reaction.get_k()*size;
}

//synchronous version of direct-method
//---start
double Compartment::get_local_propensity() {
  local_propensity_ = get_reaction_propensity(*direct_method_reactions_[0]);
  for (unsigned i(1); i < direct_method_reactions_.size(); ++i) {
    local_propensity_ += get_reaction_propensity(*direct_method_reactions_[i]);
  }
  return local_propensity_;
}
double Compartment::get_next_interval(ParallelEnvironment &pe,
                                      const double time_left) {
  double dt(std::numeric_limits<double>::infinity());
  const double old_propensity(global_propensity_);
  global_propensity_ = get_local_propensity();
  pe.getcart().Allreduce(MPI::IN_PLACE, &global_propensity_, 1,
                         MPI::DOUBLE, MPI::SUM);
  if (global_propensity_) {
    dt = old_propensity/global_propensity_*time_left;
  }
  return dt;
}
double Compartment::get_new_interval(ParallelEnvironment &pe) {
  double dt(std::numeric_limits<double>::infinity());
  global_propensity_ = get_local_propensity();
  pe.getcart().Allreduce(MPI::IN_PLACE, &global_propensity_, 1,
                         MPI::DOUBLE, MPI::SUM);
  if (global_propensity_) {
    if (!pe.getrank()) {
      dt = -log((*randdbl_)())/global_propensity_;
    }
    pe.getcart().Bcast(&dt, 1 , MPI_DOUBLE, 0);
  }
  return dt;
}
double Compartment::react_direct_method(Lattice &g, ParallelEnvironment &pe) {
  double local_propensities[pe.getsize()];
  pe.getcart().Allgather(&local_propensity_, 1, MPI::DOUBLE,
                         &local_propensities, 1, MPI::DOUBLE);
  double random((*global_randdbl_)());
  /*
  if (!pe.getrank()) {
    random = (*randdbl_)();
  }
  pe.getcart().Bcast(&random, 1 , MPI_DOUBLE, 0);
  */
  for (unsigned i(1); i < pe.getsize(); ++i) {
    local_propensities[i] += local_propensities[i-1];
  }
  double random_propensity(random*local_propensities[pe.getsize()-1]);
  if (local_propensities[pe.getrank()] >= random_propensity && 
      (!pe.getrank() || 
       local_propensities[pe.getrank()-1] < random_propensity)) {
    double accumulated_propensity(local_propensities[pe.getrank()]-
                                  local_propensity_);
    for (unsigned i(0); i < direct_method_reactions_.size() &&
         accumulated_propensity < random_propensity; ++i) {
      Reaction& reaction(*direct_method_reactions_[i]);
      accumulated_propensity += get_reaction_propensity(reaction);
      if(accumulated_propensity >= random_propensity) {
        do_direct_method_reaction(reaction, g, pe);
      }
    }
  }
  return get_new_interval(pe);
}
//---end

void Compartment::calculate_probability(Reaction &f, Lattice &g) {
  const int r0 = f.getR0();
  const int r1 = f.getR1();
  Species&  s0 = f.getS0();
  Species&  s1 = f.getS1();
  const double kAB = f.get_k();
  const double DA = s0.getD();
  const double DB = s1.getD();
  const double rv = g.getradius();
  const double SQR2 = 1.414213562;
  if(r0==1 && r1==1) {
    const double probability = kAB /(6.0*SQR2*rv*(DA+DB));
    f.set_probability(probability);
  }
}

void Compartment::calculate_max_probability(Species &s) {
  int m = influenced_reactions_.size(); 
  double maxPi = 0.0;
  for(int j=0; j<m; ++j) {
    Reaction& Rj(*influenced_reactions_[j]);
    Species& s0 = Rj.getS0();
    Species& s1 = Rj.getS1(); 
    if( &s==&s0 || &s==&s1 ) {
      const double Pi = Rj.get_probability();
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
    std::vector<Molecule>& molecules(species_molecules_[s.get_id()]);
    const unsigned begin_mol_size(molecules.size());
    unsigned index(0);
    std::vector<unsigned>& is_reactive(is_reactive_[s.get_id()]);
    std::vector<double>& reaction_probabilities(reaction_probabilities_[
                                                s.get_id()]);
    while (index < begin_mol_size && index < molecules.size()) { 
      Molecule& src_molecule(molecules[index]);
      const unsigned src_coord(src_molecule.coord);
      Voxel& src_voxel(g.get_voxel(src_coord)); 
      const unsigned adj_index((*adj_rand_)());
      const unsigned tar_coord(g.get_adjacent_coord(src_coord, adj_index));
      const int src_voxel_species_id(src_voxel.species_id);
      Voxel& tar_voxel(g.get_voxel(tar_coord));
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
        ghost_species_ids.push_back(s.get_id());
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
  const unsigned reaction_id(influenced_reaction_ids_[s.get_id()][tarID]);
  Reaction& reaction(*influenced_reactions_[reaction_id]);
  Species& s0(reaction.getS0());
  if (&s0 == &s) {
    do_reaction(g, reaction, src_voxel, tar_voxel, src_coord, tar_coord);
  }
  else {
    do_reaction(g, reaction, tar_voxel, src_voxel, tar_coord, src_coord);
  }
}

void Compartment::do_reaction(Lattice& g, Reaction& r, Voxel& v0, Voxel& v1,
                              const unsigned c0, const unsigned c1) {
  Species& s0(r.getS0());
  Species& s1(r.getS1());
  if (&s0 != &s1) {
    if (r.getP2()) {
      Species& p0(r.getS2());
      if (&p0 != &s0 && &p0 != &s1) {
        //must remove before adding since in the removed species, the replacing
        //molecule might be the one removed and it will be occupied
        //by the new molecule:
        remove_molecule(g, s0.get_id(), v0);
        add_molecule(p0.get_id(), c0, get_new_mol_id(), v0);
      }
    }
    else {
      remove_molecule(g, s0.get_id(), v0);
    }
    if (r.getP3()) {
      Species& p1(r.getS3());
      if (&p1 != &s1 && &p1 != &s0) {
        //must remove before adding since in the removed species, the replacing
        //molecule might be the one removed and it will be occupied
        //by the new molecule:
        remove_molecule(g, s1.get_id(), v1);
        add_molecule(p1.get_id(), c1, get_new_mol_id(), v1);
      }
    }
    else {
      remove_molecule(g, s1.get_id(), v1);
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
  Voxel& back_voxel(g.get_voxel(back_outmolecule.coord));
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
  Voxel& back_voxel(g.get_voxel(back_coord));
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


void Compartment::calculate_collision_time(Species& s,
                                           ParallelEnvironment& pe) {
  s.calcCollisionTime();
  const unsigned id0(s.get_id());
  for (unsigned i(0); i < is_reactive_[id0].size(); ++i) {
    if (is_reactive_[id0][i]) {
      const unsigned reaction_id(influenced_reaction_ids_[id0][i]);
      Reaction& reaction(*influenced_reactions_[reaction_id]);
      double probability(reaction.get_probability()*s.get_walk_probability());
      reaction_probabilities_[id0][i] = probability;
      if (!pe.getrank()) {
        std::cout << "reaction:" << reaction.get_name() << std::endl;
        std::cout << "\tspecies:" << s.get_name() << std::endl;
        std::cout << "\tp:" << probability << std::endl;
        std::cout << "\twalk probability:" << s.get_walk_probability() <<
          std::endl;
        std::cout << "\twalk interval:" << s.get_walk_interval() <<
          std::endl;
      }
    }
  }
}

void Compartment::add_diffusion_influenced_reaction(Reaction& r) {
  Species& s0(r.getS0());
  Species& s1(r.getS1());
  const unsigned reaction_id(influenced_reactions_.size());
  influenced_reactions_.push_back(&r);
  const unsigned id0(s0.get_id());
  const unsigned id1(s1.get_id());
  influenced_reaction_ids_[id0][id1] = reaction_id;
  influenced_reaction_ids_[id1][id0] = reaction_id;
  is_reactive_[id0][id1] = 1;
  is_reactive_[id1][id0] = 1;
}

unsigned Compartment::check_species_size(Lattice& g, ParallelEnvironment& pe,
                                         const unsigned sid,
                                         const double id) {
  Species& s(*species_list_[sid]);
  const std::vector<Molecule>& molecules(species_molecules_[s.get_id()]);
  unsigned size(molecules.size());
  pe.getcart().Allreduce(MPI::IN_PLACE, &size, 1, MPI::INT, MPI::SUM);
  fout << id << " " << s.get_name() << " size:" << size << std::endl;
  return size;
}


void Compartment::walk_on_ghost(std::vector<unsigned>& species_ids,
                                std::vector<unsigned>& src_coords,
                                std::vector<unsigned>& tar_coords, Lattice& g,
                                ParallelEnvironment& pe) { 
  //prepare sub-vector of molecules for each sub-volume 
  std::vector<unsigned> sub_indices[8];
  for (unsigned i(0); i < 8; ++i) {
    sub_indices[i].reserve(src_coords.size()/4);
  }

  for (unsigned i(0); i < src_coords.size(); ++i) {
    sub_indices[coord_to_subvolume(g, src_coords[i])].push_back(i);
  } 
  
  std::shuffle(order_.begin(), order_.end(), global_rng_);
  std::vector<unsigned> local_order(8);
  if (!pe.getrank()) {
    local_order = order_;
  }
  pe.getcart().Bcast(local_order.data(), 8, MPI::INT, 0); 
  for (unsigned i(0); i < 8; ++i) {
    if (local_order[i] != order_[i]) {
      std::cout << "error in global seed:" << pe.getrank() << std::endl;
    }
  }
  
  for(unsigned sv(0); sv < 8; ++sv) {
    g.jumpincoords.clear(); 
    std::vector<SpillMolecule> spill_coords;
    const unsigned N(order_[sv]);
    // load ghost cells from adjacent process using MPI 
    // this is the costliest operation, need to further optimize
    // by keeping a list of molecules that are in the in-compartment-ghost:
    g.load_ghost(pe, outmolecules_[N+1], outmolecules_[9], N);
    // loop over sub-moleculeVector
    for (unsigned i(0); i < sub_indices[N].size(); ++i) {
      const unsigned index(sub_indices[N][i]);
      const unsigned tar_coord(tar_coords[index]);
      const unsigned src_coord(src_coords[index]);
      Voxel& tar_voxel(g.get_voxel(tar_coord));
      int tar_voxel_species_id(tar_voxel.species_id);
      Voxel& src_voxel(g.get_voxel(src_coord));
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
  const unsigned reaction_id(influenced_reaction_ids_[s.get_id()][tarID]);
  Reaction& reaction(*influenced_reactions_[reaction_id]);
  Species& s0(reaction.getS0());
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
  Species& s0(r.getS0());
  Species& s1(r.getS1());
  if (&s0 != &s1) {
    if (r.getP2()) {
      Species& p0(r.getS2());
      if (&p0 != &s0 && &p0 != &s1) {
        //must remove before adding since in the removed species, the replacing
        //molecule might be the one removed and it will be occupied
        //by the new molecule:
        remove_molecule_on_ghost(g, s0.get_id(), c0, v0, spill_coords);
        add_molecule_on_ghost(p0.get_id(), c0, v0, spill_coords);
      }
    }
    else {
      remove_molecule_on_ghost(g, s0.get_id(), c0, v0, spill_coords);
    }
    if (r.getP3()) {
      Species& p1(r.getS3());
      if (&p1 != &s1 && &p1 != &s0) {
        //must remove before adding since in the removed species, the replacing
        //molecule might be the one removed and it will be occupied
        //by the new molecule:
        remove_molecule_on_ghost(g, s1.get_id(), c1, v1, spill_coords);
        add_molecule_on_ghost(p1.get_id(), c1, v1, spill_coords);
      }
    }
    else {
      remove_molecule_on_ghost(g, s1.get_id(), c1, v1, spill_coords);
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
    Voxel& v(g.get_voxel(coord));
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

void Compartment::add_coordinates_species(Species& species) {
  output_coord_species_.push_back(species);
}

void Compartment::add_number_species(Species& species) {
  output_number_species_.push_back(species);
}


/*
std::vector<float> get_logspace_vector(float a, float b, int k) {
  const auto exp_scale = (b - a) / (k - 1);
  std::vector<float> logspace;
  logspace.reserve(k);
  for (int i = 0; i < k; i++) {
    logspace.push_back(i * exp_scale);
  }
  std::for_each(logspace.begin(), logspace.end(), [](float &x) {
    x = pow(10, x); });
  return logspace;
}
*/

/*
std::vector<float> get_logspace_vector(float a, float b, int k) {
  std::vector<float> logspace;
  for (int i = 0; i < k; i++) {
    logspace.push_back(pow(10, i * (b - a) / (k - 1)));
  }
  return logspace;
}
*/

template<typename T>
class Logspace {
private:
    T curValue, base;

public:
    Logspace(T first, T base) : curValue(first), base(base) {}

    T operator()() {
        T retval = curValue;
        curValue *= base;
        return retval;
    }
};

std::vector<float> get_logspace_vector(float start, float stop, int num = 50,
                                       float base = 10) {
  float realStart = pow(base, start);
  float realBase = pow(base, (stop-start)/num); 
  std::vector<float> retval;
  retval.reserve(num);
  std::generate_n(std::back_inserter(retval), num,
                  Logspace<float>(realStart,realBase));
  return retval;
}

void Compartment::output_coordinates_header(Lattice& lattice,
                                          ParallelEnvironment &pe,
                                          const float dt,
                                          const float start_time,
                                          const float end_time,
                                          const unsigned n_logs) {
  if (dt < 0) {
    output_coords_dt_ = std::numeric_limits<double>::infinity();
    for (unsigned i(0); i < output_coord_species_.size(); ++i) {
      Species& s(output_coord_species_[i]);
      if (s.get_walk_interval() < output_coords_dt_) {
        output_coords_dt_ = s.get_walk_interval();
      }
    }
  }
  else {
    output_coords_dt_ = dt;
  }
  if (n_logs) {
    coords_logspace_ = get_logspace_vector(start_time, end_time, n_logs);
  }
  fout2 << "Time";
  for (unsigned i(0); i < output_coord_species_.size(); ++i) {
    Species& s(output_coord_species_[i]);
    fout2 << "," << s.get_name();
  }
  Vector<unsigned> gdims(pe.get_global_dimensions());
  fout2 << "," << lattice.getradius() << " " << pe.getsize() << " " << 
    gdims.x << " " << gdims.y << " " << gdims.z << endl;
}

float Compartment::output_coordinates(const double current_time) {
  fout2 << setprecision(15) << current_time;
  for (unsigned i(0); i < output_coord_species_.size(); ++i) {
    Species& s(output_coord_species_[i]);
    std::vector<Molecule>& molecules(species_molecules_[s.get_id()]);
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
  if (coords_logspace_.size()) {
    if (coords_logspace_cnt_+1 >= coords_logspace_.size()) {
      return std::numeric_limits<float>::infinity();
    }
    float curr_time(coords_logspace_[coords_logspace_cnt_]);
    float next_time(coords_logspace_[coords_logspace_cnt_+1]);
    coords_logspace_cnt_++;
    return next_time-curr_time;
  }
  return output_coords_dt_;
}


void Compartment::output_numbers_header(ParallelEnvironment& pe,
                                        const float dt, const float start_time,
                                        const float end_time,
                                        const unsigned n_logs) {
  if (dt < 0) {
    output_numbers_dt_ = std::numeric_limits<double>::infinity();
    for (unsigned i(0); i < output_coord_species_.size(); ++i) {
      Species& s(output_coord_species_[i]);
      if (s.get_walk_interval() < output_numbers_dt_) {
        output_numbers_dt_ = s.get_walk_interval();
      }
    }
  }
  else {
    output_numbers_dt_ = dt;
  }
  if (n_logs) {
    numbers_logspace_ = get_logspace_vector(start_time, end_time, n_logs);
  }
  fout3 << setprecision(15) << "Time " << pe.getsize();
  for (unsigned i(0); i < output_number_species_.size(); ++i) {
    Species& s(output_number_species_[i]);
    fout3 << "," << s.get_name();
  }
  fout3 << endl;
}


float Compartment::output_numbers(const double current_time) {
  fout3 << setprecision(15) << current_time;
  for (unsigned i(0); i < output_number_species_.size(); ++i) {
    Species& s(output_number_species_[i]);
    fout3 << "," << species_molecules_[s.get_id()].size();
  }
  fout3 << endl;
  if (numbers_logspace_.size()) {
    if (numbers_logspace_cnt_ >= numbers_logspace_.size()) {
      return std::numeric_limits<float>::infinity();
    }
    float curr_time(numbers_logspace_[numbers_logspace_cnt_]);
    float next_time(numbers_logspace_[numbers_logspace_cnt_+1]);
    numbers_logspace_cnt_++;
    return next_time-curr_time;
  }
  return output_numbers_dt_;
}

/*
//asynchronous version of direct-method
//---start
double Compartment::get_local_propensity() {
  double propensity(get_reaction_propensity(*direct_method_reactions_[0]));
  for (unsigned i(1); i < direct_method_reactions_.size(); ++i) {
    propensity += get_reaction_propensity(*direct_method_reactions_[i]);
  }
  return propensity;
}

double Compartment::get_next_interval(ParallelEnvironment &pe,
                                      const double time_left) {
  double dt(std::numeric_limits<double>::infinity());
  const double old_propensity(local_propensity_);
  local_propensity_ = get_local_propensity();
  if (local_propensity_) {
    dt = old_propensity/local_propensity_*time_left;
  }
  return dt;
}

double Compartment::get_new_interval(ParallelEnvironment &pe) {
  double dt(std::numeric_limits<double>::infinity());
  local_propensity_ = get_local_propensity();
  if (local_propensity_) {
    dt = -log((*randdbl_)())/local_propensity_;
  }
  return dt;
}

double Compartment::react_direct_method(Lattice &g, ParallelEnvironment &pe) {
  double random_propensity((*randdbl_)()*local_propensity_);
  double accumulated_propensity(0);
  for (unsigned i(0); i < direct_method_reactions_.size() &&
         accumulated_propensity < random_propensity; ++i) {
    Reaction& reaction(*direct_method_reactions_[i]);
    accumulated_propensity += get_reaction_propensity(reaction);
    if(accumulated_propensity >= random_propensity) {
      do_direct_method_reaction(reaction, g, pe);
    }
  }
  return get_new_interval(pe);
}
//---end
*/

