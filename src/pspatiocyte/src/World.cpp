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

#include "World.hpp"

// parallel file streams
ofstream fout, fout2, fout3;

World::World(int argc, char* argv[], const unsigned Nx,
             const int Ny, const unsigned Nz, const double rv,
             const float output_numbers_dt,
             const std::string dirname,
             const unsigned seed,
             const float output_coords_dt,
             const bool is_force_search_vacant):
  invalid_species_("Invalid", 0, 0, *this),
  vacant_species_("Vacant", 0, 0, *this),
  ghost_id_(-1),
  output_numbers_dt_(output_numbers_dt),
  output_coords_dt_(output_coords_dt),
  parallel_environment_(argc, argv, Nx, Ny, Nz, dirname),
  lattice_("SchisMatrix", rv, parallel_environment_, invalid_species_.get_id(),
               vacant_species_.get_id(), ghost_id_),
  compartment_("Cell", VOLUME, //rand(),
               scheduler_,
               seed*parallel_environment_.getsize()+
               parallel_environment_.getrank(),
               parallel_environment_.getsize(), parallel_environment_.getrank(),
               invalid_species_.get_id(), vacant_species_.get_id(), ghost_id_,
               is_force_search_vacant) {
  }

int DistributeMolecule(int Total, int Nproc, int rank) {
  const int Ndist = (int)(Total/Nproc); 
  // distribution of the remainder
  if( rank < (Total%Nproc) ) {
    return (Ndist+1);
  }
  else {
    return Ndist;
  }
}

double World::get_current_time() {
  return scheduler_.get_time();
}

void World::add_output_numbers_species(Species& species) {
  output_numbers_species_list_.push_back(&species);
}

void World::add_output_coords_species(Species& species) {
  output_coords_species_list_.push_back(&species);
}

void World::set_output_numbers_logspace(const float start_time,
                                        const float end_time,
                                        const unsigned n_logs) {
  output_numbers_start_ = start_time;
  output_numbers_end_ = end_time;
  output_numbers_n_logs_ = n_logs;
  if (output_numbers_dt_ <= 0) {
    output_numbers_dt_ = 1e-7;
  }
}

void World::set_output_coords_logspace(const float start_time,
                                       const float end_time,
                                       const unsigned n_logs) {
  output_coords_start_ = start_time;
  output_coords_end_ = end_time;
  output_coords_n_logs_ = n_logs;
  if (output_coords_dt_ <= 0) {
    output_coords_dt_ = 1e-7;
  }
}

void World::initialize() {
  out_species_= new Species("Out", 0, 0, *this),
  compartment_.set_out_id(out_species_->get_id());
  lattice_.set_out_id(out_species_->get_id());
  compartment_.initialize(lattice_, parallel_environment_, species_list_);
  const int nproc(parallel_environment_.getsize());
  const int myrank(parallel_environment_.getrank());
  for (unsigned i(0); i < species_list_.size(); ++i) { 
    Species& species(*species_list_[i]);
    unsigned init_size(species.get_init_size());
    if (init_size) {
      init_size = DistributeMolecule(init_size, nproc, myrank);
    }
    compartment_.populate_molecules(species, init_size, lattice_,
                                    parallel_environment_);
  }

  for (unsigned i(0); i < independent_reactions_.size(); ++i) {
    compartment_.add_direct_method_reaction(*independent_reactions_[i]);
  }

  std::vector<Species*> reactants;
  for (unsigned i(0); i < influenced_reactions_.size(); ++i) {
    Reaction& reaction(*influenced_reactions_[i]);
    compartment_.add_diffusion_influenced_reaction(reaction);
    compartment_.calculate_probability(reaction, lattice_);
    if (std::find(reactants.begin(), reactants.end(), reaction.getS0()) == 
        reactants.end()) {
      reactants.push_back(reaction.getS0());
    }
    if (std::find(reactants.begin(), reactants.end(), reaction.getS1()) == 
        reactants.end()) {
      reactants.push_back(reaction.getS1());
    }
  }

  for (unsigned i(0); i < reactants.size(); ++i) {
    compartment_.calculate_max_probability(*reactants[i]);
  }

  for (unsigned i(0); i < reactants.size(); ++i) {
    compartment_.calculate_collision_time(*reactants[i], parallel_environment_);
  }

  if (output_numbers_dt_) {
    if (!output_numbers_species_list_.size()) {
      output_numbers_species_list_ = species_list_;
    }
    for (unsigned i(0); i < output_numbers_species_list_.size(); ++i) { 
      Species& species(*output_numbers_species_list_[i]);
      if (&species != &vacant_species_ && &species != &invalid_species_ &&
          &species != out_species_) {
        compartment_.add_number_species(species);
      }
    }
    compartment_.output_numbers_header(parallel_environment_,
                                       output_numbers_dt_, 
                                       output_numbers_start_,
                                       output_numbers_end_,
                                       output_numbers_n_logs_);
    float dt(compartment_.output_numbers(scheduler_.get_time()));
    scheduler_.add_event(SpatiocyteEvent(lattice_, compartment_,
                                 parallel_environment_, OUTPUT_NUMBERS, dt));
  }

  if (output_coords_dt_) {
    if (!output_coords_species_list_.size()) {
      output_coords_species_list_ = species_list_;
    }
    for (unsigned i(0); i < output_coords_species_list_.size(); ++i) { 
      Species& species(*output_coords_species_list_[i]);
      if (&species != &vacant_species_ && &species != &invalid_species_ &&
          &species != out_species_) {
        compartment_.add_coordinates_species(species);
      }
    }
    compartment_.output_coordinates_header(lattice_, parallel_environment_,
                                           output_coords_dt_,
                                           output_coords_start_,
                                           output_coords_end_,
                                           output_coords_n_logs_);
    float dt(compartment_.output_coordinates(scheduler_.get_time()));
    scheduler_.add_event(SpatiocyteEvent(lattice_, compartment_,
                               parallel_environment_, OUTPUT_COORDINATES, dt));
  }

  for (unsigned i(0); i < species_list_.size(); ++i) {
    Species& species(*species_list_[i]);
    if (species.getD()) {
      scheduler_.add_event(SpatiocyteEvent(lattice_, compartment_,
                                           parallel_environment_, species));
    }
  }

  if (independent_reactions_.size()) {
    const int independent_event_index_(
      scheduler_.add_event(SpatiocyteEvent(lattice_, compartment_,
                                          parallel_environment_)));

    compartment_.set_independent_event_index(independent_event_index_);
  }
  //queue events should be executed after pushing all the events in the
  //scheduler, otherwise the event pointers in queue will become invalid:
  scheduler_.queue_events();
}

void World::run(const double end_time, const unsigned verbose) {
  double interval(end_time/10);
  double prev_time(scheduler_.get_time());
  unsigned n(0);
  parallel_environment_.starttimer();
  if(!parallel_environment_.getrank() && verbose) {
    std::cout << std::endl;
  }
  while(scheduler_.get_time() < end_time) {
    /*
    if (!parallel_environment_.getrank() && scheduler_.get_time() > 62) {
      fout << "step: time:" << scheduler_.get_time() << " top time:" <<
        scheduler_.get_top_time() << " name:" << 
        scheduler_.get_top_event().get_name() << std::endl;
    }
    */
    scheduler_.step();
    if(scheduler_.get_time()-prev_time > interval) {
      prev_time = scheduler_.get_time(); 
      if(!parallel_environment_.getrank() && verbose) {
        std::cout << "t/duration: " << scheduler_.get_time() << "/" <<
          end_time << std::endl;
      }
      ++n;
    }
  }
  parallel_environment_.stoptimer();
}

int World::add_species(Species* newspec) {
  species_list_.push_back(newspec); 
  return species_list_.size()-1;
}

int World::add_independent_reaction(Reaction* reaction) {
  independent_reactions_.push_back(reaction);
  return independent_reactions_.size()-1;
}

int World::add_influenced_reaction(Reaction* reaction) {
  influenced_reactions_.push_back(reaction);
  return influenced_reactions_.size()-1;
}

