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


#ifndef __WORLD_HPP
#define __WORLD_HPP

#include <string>
#include <vector>
#include <map>
#include <time.h>
#include "Common.hpp"
#include "Lattice.hpp"
#include "Compartment.hpp"
#include "Species.hpp"
#include "Reaction.hpp"
#include "SpatiocyteEvent.hpp"

class World {
public:
  World(int argc, char* argv[], const unsigned Nx, const int Ny,
        const unsigned Nz, const double rv,
        const float output_numbers_dt=0,
        const std::string dirname="output",
        const unsigned seed=time(NULL),
        const float output_coords_dt=0,
        const bool is_force_search_vacant=false);
  ~World() {} 
  int add_species(Species* newspec); 
  int add_independent_reaction(Reaction* reaction); 
  int add_influenced_reaction(Reaction* reaction);
  double get_current_time();
  void initialize();
  void run(const double end_time, const unsigned verbose);
  void add_output_numbers_species(Species& species);
  void add_output_coords_species(Species& species);
  void set_output_numbers_logspace(const float start_time,
                                   const float end_time, const unsigned n_logs);
  void set_output_coords_logspace(const float start_time, const float end_time,
                                  const unsigned n_logs);
private:
  ParallelEnvironment parallel_environment_;
  std::vector<Species*> species_list_;
  Species invalid_species_;
  Species vacant_species_;
  const int ghost_id_;
  float output_numbers_dt_;
  float output_coords_dt_;
  Lattice lattice_;
  EventScheduler<SpatiocyteEvent> scheduler_;
  Compartment compartment_;
  Species* out_species_;
  std::vector<Reaction*> independent_reactions_;
  std::vector<Reaction*> influenced_reactions_;
  std::vector<Species*> output_numbers_species_list_;
  std::vector<Species*> output_coords_species_list_;
  unsigned output_coords_n_logs_ = 0;
  unsigned output_numbers_n_logs_ = 0;
  float output_coords_end_;
  float output_numbers_end_;
  float output_coords_start_;
  float output_numbers_start_;
};

#endif /* __WORLD_HPP */
