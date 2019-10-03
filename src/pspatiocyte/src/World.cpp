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
  lattice_("SchisMatrix", rv, parallel_environment_, invalid_species_.getID(),
               vacant_species_.getID(), ghost_id_),
  compartment_("Cell", VOLUME, //rand(),
               seed*parallel_environment_.getsize()+
               parallel_environment_.getrank(),
               parallel_environment_.getsize(), parallel_environment_.getrank(),
               invalid_species_.getID(), vacant_species_.getID(), ghost_id_,
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
  return scheduler_.getTime();
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
  if (!output_numbers_dt_) {
    output_numbers_dt_ = 1e-7;
  }
}

void World::set_output_coords_logspace(const float start_time,
                                       const float end_time,
                                       const unsigned n_logs) {
  output_coords_start_ = start_time;
  output_coords_end_ = end_time;
  output_coords_n_logs_ = n_logs;
  if (!output_coords_dt_) {
    output_coords_dt_ = 1e-7;
  }
}

void World::initialize() {
  out_species_= new Species("Out", 0, 0, *this),
  compartment_.set_out_id(out_species_->getID());
  lattice_.set_out_id(out_species_->getID());
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

  for (unsigned i(0); i < species_list_.size(); ++i) {
    Species& species(*species_list_[i]);
    if (species.getD()) {
      scheduler_.addEvent(SpatiocyteEvent(&lattice_, &compartment_,  &species));
    }
  }

  if (independent_reactions_.size()) {
    independent_event_id_ =
      scheduler_.addEvent(SpatiocyteEvent(&lattice_, &compartment_,
                                          &parallel_environment_));
  }

  if (output_numbers_dt_ > 0) {
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
    float dt(compartment_.output_numbers(scheduler_.getTime()));
    scheduler_.addEvent(SpatiocyteEvent(&lattice_, &compartment_,
                                  &parallel_environment_, OUTPUT_NUMBERS, dt));
  }

  if (output_coords_dt_ > 0) {
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
    float dt(compartment_.output_coordinates(scheduler_.getTime()));
    scheduler_.addEvent(SpatiocyteEvent(&lattice_, &compartment_,
                              &parallel_environment_, OUTPUT_COORDINATES, dt));
  }

}

void World::run(const double end_time, const unsigned verbose) {
  double interval(end_time/10);
  double prev_time(scheduler_.getTime());
  unsigned n(0);
  parallel_environment_.starttimer();
  while(scheduler_.getTime() < end_time) {
    /*
    std::cout << "step: time:" << scheduler_.getTime() << " top time:" <<
      scheduler_.getTopTime() << " name:" << 
      scheduler_.getTopEvent().getName() << std::endl;
      */
    if(independent_event_id_ >= 0 && 
       (scheduler_.getTime() < scheduler_.getTopTime()) ||
       (independent_event_id_ == scheduler_.getTopID() &&
        scheduler_.getTime() == scheduler_.getTopTime())) {
      scheduler_.updateEventTime(independent_event_id_,
       compartment_.get_next_time(parallel_environment_, scheduler_.getTime()));
    }
    scheduler_.step(parallel_environment_);
    if(scheduler_.getTime()-prev_time > interval) {
      prev_time = scheduler_.getTime(); 
      if(!parallel_environment_.getrank() && verbose) {
        std::cout << "current t/duration: " << scheduler_.getTime() << "/" << end_time <<
          std::endl;
      }
      ++n;
    }
  }
  const double last_time(scheduler_.getTime());
  if(!parallel_environment_.getrank() && verbose) {
    std::cout << "current t:" << last_time << std::endl;
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

