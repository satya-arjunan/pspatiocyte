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


#ifndef __COMPARTMENT_HPP
#define __COMPARTMENT_HPP

#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include "Common.hpp"
#include "Molecule.hpp"
#include "Species.hpp"
#include "Lattice.hpp"
#include "Reaction.hpp"
#include "ParallelEnvironment.hpp"
#include "Vector.hpp"
#include "SpatiocyteEvent.hpp"
#include "EventScheduler.hpp"
#include <boost/random.hpp>
#include <limits>
#include <random>

using namespace std;
using namespace boost;

class Compartment {
public:
  Compartment(string name, COMPARTMENT_TYPE type, 
              uint32_t seed,
              const int proc_size, const int proc_id, const int invalid_id,
              const int vacant_id, const int ghost_id,
              const bool is_force_search_vacant):
    is_parallel_(proc_size > 1),
    invalid_id_(invalid_id),
    vacant_id_(vacant_id),
    ghost_id_(ghost_id),
    is_force_search_vacant_(is_force_search_vacant),
    name_(name),
    type_(type),
    rng_(seed),
    gen_(seed),
    local_propensity_(0),
    mol_id_min_(std::numeric_limits<unsigned>::max()/proc_size*proc_id),
    mol_id_max_(mol_id_min_+std::numeric_limits<unsigned>::max()/proc_size),
    mol_id_(mol_id_min_) {
      randdbl_ = new boost::variate_generator<
        boost::mt19937&,boost::uniform_real<>>(gen_, uniform_real<>());
      adj_rand_ = new boost::variate_generator<
        boost::mt19937&, boost::uniform_int<>>(gen_, dist_);
  }

  ~Compartment() {
    delete randdbl_;
    delete adj_rand_;
  }

  std::string get_name() {
    return name_;
  }

  COMPARTMENT_TYPE get_type() {
    return type_;
  }

  int get_id() {
    return compartment_id_;
  }

  void set_id(int id) {
    compartment_id_ = id;
  }

  void set_out_id(const int out_id) {
    out_id_ = out_id;
  }
  
  unsigned get_new_mol_id() {
    ++mol_id_;
    if (mol_id_ > mol_id_max_) {
      mol_id_ = mol_id_min_;
    }
    return mol_id_;
  }

  unsigned coord_to_subvolume(Lattice& g, const unsigned coord);
  bool walk_molecule(Lattice& g, const int src_coord, const int curr_coord);
  void initialize(Lattice &g, ParallelEnvironment &pe,
                  std::vector<Species*>& species_list);
  void populate_molecules(Species& s, const unsigned size, Lattice &g,
                          ParallelEnvironment &pe);
  void populate_molecules(Species& s, const unsigned size, Lattice &g);
  void add_direct_method_reaction(Reaction& r);
  bool remove_reactant(Species &s, Lattice &g, bool binary, int &center,
                       int &neighbor);
  void add_product(Species &s, Lattice &g, int coord);
  double get_reaction_propensity(Reaction &r);
  double get_local_propensity();
  double get_new_interval(ParallelEnvironment &pe);
  double get_next_interval(ParallelEnvironment &pe, const double time_left);
  double react_direct_method(Lattice &g, ParallelEnvironment &pe);
  void add_diffusion_influenced_reaction(Reaction& r);
  void calculate_probability(Reaction &r, Lattice &g);
  void calculate_max_probability(Species &s);
  void calculate_collision_time(Species &s, ParallelEnvironment &pe);
  void walk(std::vector<Species*>& species_list, Lattice &g,
            ParallelEnvironment &pe);

  void walk_on_ghost(std::vector<unsigned>& species_ids,
                     std::vector<unsigned>& src_coords,
                     std::vector<unsigned>& tar_coords,
                     Lattice& g, ParallelEnvironment& pe);
  bool do_direct_method_reaction(Reaction &r, Lattice &g,
                                      ParallelEnvironment &pe);
  bool get_vacant_neighbor(Lattice &g, const int &coord, int &center,
                           int &neighbor);
  void add_coordinates_species(Species& species);
  void add_number_species(Species& species);
  void output_coordinates_header(Lattice& lattice, ParallelEnvironment &pe,
                                 const float dt, const float start_time,
                                 const float end_time, const unsigned n_logs);
  float output_coordinates(const double current_time);
  void output_numbers_header(ParallelEnvironment& pe, const float dt,
                             const float start_time, const float end_time,
                             const unsigned n_logs);
  float output_numbers(const double current_time);

  void react(Lattice& g, Species& s, const unsigned tarID, Voxel& src_voxel,
             Voxel& tar_voxel, const unsigned src_coord,
             const unsigned tar_coord);
  void react_on_ghost(Lattice& g, Species& s, const unsigned tarID,
                      Voxel& src_voxel, Voxel& tar_voxel,
                      const unsigned src_coord, const unsigned tar_coord, 
                      std::vector<SpillMolecule>& spill_coords);
  void do_reaction(Lattice& g, Reaction& r, Voxel& v0, Voxel& v1,
                   const unsigned c0, const unsigned c1);
  void do_reaction_on_ghost(Lattice& g, Reaction& r, Voxel& v0, Voxel& v1,
                            const unsigned c0, const unsigned c1,
                            std::vector<SpillMolecule>& spill_coords);
  void add_molecule(const int species_id, const unsigned coord,
                    const unsigned mol_id, Voxel& voxel);
  unsigned remove_molecule(Lattice& g, const int species_id, Voxel& voxel);
  void populate_outmolecule_in_voxel(const int out_id,
                                     const int species_id,
                                     const unsigned mol_index,
                                     const unsigned coord, Voxel& voxel);
  void vacate_outmolecule_in_voxel(Lattice& g, const int out_id, Voxel& voxel);

  unsigned remove_molecule_from_molecules(Lattice& g, const int species_id,
                                          const unsigned mol_index);
  void add_molecule_on_ghost(const int species_id, const unsigned coord,
                             Voxel& v,
                             std::vector<SpillMolecule>& spill_coords);
  void remove_molecule_on_ghost(Lattice& g, const int species_id,
                                const unsigned coord, Voxel& v,
                                std::vector<SpillMolecule>& spill_coords);
  void check_voxels(Lattice& g, const double id);
  unsigned check_species_size(Lattice& g, ParallelEnvironment& pe,
                              const unsigned sid, const double id);
  void spill_molecules(Lattice& g, const std::vector<SpillMolecule>& coord,
                      const Vector<unsigned>& ghost);
  void jumpin_molecules(Lattice& g, ParallelEnvironment& pe);

private:
  const int is_parallel_;
  const int invalid_id_;
  const int vacant_id_;
  const int ghost_id_;
  const bool is_force_search_vacant_;
  unsigned nx_;
  unsigned ny_;
  unsigned nz_;
  unsigned n_out_voxels_ = 0;
  unsigned n_out_ghost_voxels_ = 0;
  unsigned n_ghost_voxels_ = 0;
  unsigned out_cnts_[10];
  unsigned out_ghost_cnts_[10];
  int cntES = 0;
  int out_id_ = vacant_id_; //default value
  std::mt19937 rng_;
  std::vector<Species*> species_list_;
  std::vector<OutMolecule> outmolecules_[10];
  std::vector<std::vector<Molecule>> species_molecules_;
  string name_;                        // name of compartment
  COMPARTMENT_TYPE type_;              // volume or surface
  int compartment_id_;                  // compartment ID
  int n_voxels_;                  // number of voxels
  vector<int> voxel_coords_;            // linear coordinates of voxels
  vector<Species> output_coord_species_;
  vector<Species> output_number_species_;
  vector<Reaction*> direct_method_reactions_;  // list of reactions
  vector<Reaction*> influenced_reactions_;    // list of reactions

  double global_volume_;
  double local_volume_;
  double local_propensity_;
  double global_propensity_;

  boost::uniform_int<> dist_ = boost::uniform_int<>(0, 11);
  boost::mt19937 gen_;
  boost::variate_generator<boost::mt19937&, boost::uniform_int<>>* adj_rand_;
  boost::variate_generator<boost::mt19937&, boost::uniform_real<>>* randdbl_;

  void writeMolecule(vector<Molecule> &mv);
  std::vector<std::vector<unsigned>> is_reactive_;
  std::vector<std::vector<unsigned>> influenced_reaction_ids_;
  std::vector<std::vector<double>> reaction_probabilities_;
  Vector<unsigned> ghosts_[8];
  Vector<int> mid_span_;
  Vector<int> begin_[8];
  Vector<int> end_[8];
  Vector<int> bit_[8];
  const unsigned mol_id_min_ = 0;
  const unsigned mol_id_max_ = UINT_MAX;
  unsigned mol_id_;
  float output_coords_dt_;
  float output_numbers_dt_;
  std::vector<float> coords_logspace_;
  std::vector<float> numbers_logspace_;
  unsigned coords_logspace_cnt_ = 0;
  unsigned numbers_logspace_cnt_ = 0;
  std::vector<unsigned> order_ = {0, 1, 2, 3, 4, 5, 6, 7};
};

#endif /* __COMPARTMENT_HPP */

