#include "World.hpp"

using namespace std;

// parallel file streams
ofstream fout, fout2, fout3;

//distribute molecules thrown in for each process
World::World(int argc, char* argv[], const unsigned Nx,
             const int Ny, const unsigned Nz, const double rv):
  invalid_species_("Invalid", 0, 0, *this),
  vacant_species_("Vacant", 0, 0, *this),
  ghost_id_(-1),
  parallel_environment_(argc, argv, Nx, Ny, Nz),
  lattice_("SchisMatrix", rv, parallel_environment_, invalid_species_.getID(),
               vacant_species_.getID(), ghost_id_),
  compartment_("Cell", VOLUME, rand(), parallel_environment_.getsize(),
               parallel_environment_.getrank(), invalid_species_.getID(),
               vacant_species_.getID(), ghost_id_) {
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
    compartment_.throwinMolecules(species, init_size, lattice_);
  }

  for (unsigned i(0); i < independent_reactions_.size(); ++i) {
    compartment_.attachIndependentReaction(*independent_reactions_[i]);
  }

  std::vector<Species*> reactants;
  for (unsigned i(0); i < influenced_reactions_.size(); ++i) {
    Reaction& reaction(*influenced_reactions_[i]);
    compartment_.attachInfluencedReaction(reaction);
    compartment_.calculateProbability(reaction, lattice_);
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
    compartment_.findMaxProbability(*reactants[i]);
  }

  for (unsigned i(0); i < reactants.size(); ++i) {
    compartment_.calculateCollisionTime(*reactants[i]);
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
}

void World::run(const double log_interval, const double end_time,
                const unsigned verbose) {
  for (unsigned i(0); i < species_list_.size(); ++i) { 
    Species& species(*species_list_[i]);
    if (&species != &vacant_species_ && &species != &invalid_species_ &&
        &species != out_species_) {
      compartment_.addCoordinatesSpecies(species);
      compartment_.addNumberSpecies(species);
    }
  }
  double current_time(scheduler_.getTime());
  compartment_.outputCoordinatesHeader(lattice_, parallel_environment_);
  compartment_.outputNumbersHeader();
  double prev_time(scheduler_.getTime());
  compartment_.outputCoordinates(prev_time);
  compartment_.outputNumbers(prev_time);
  unsigned n(0);
  parallel_environment_.starttimer();
  while(scheduler_.getTime() < end_time) {
    /*
    std::cout << "step: time:" << scheduler_.getTime() << " top time:" <<
      scheduler_.getTopTime() << " name:" << 
      scheduler_.getTopEvent().getName() << std::endl;
      */
    if(scheduler_.getTime() < scheduler_.getTopTime() ||
       (independent_event_id_ == scheduler_.getTopID() &&
        scheduler_.getTime() == scheduler_.getTopTime())) {
      scheduler_.updateEventTime(independent_event_id_,
         compartment_.getNextTime(parallel_environment_, scheduler_.getTime()));
    }
    scheduler_.step(parallel_environment_);
    if(scheduler_.getTime()-prev_time > log_interval) {
      prev_time = scheduler_.getTime(); 
      //compartment_.outputCoordinates(prev_time);
      compartment_.outputNumbers(prev_time);
      if(!parallel_environment_.getrank() && verbose) {
        std::cout << "n:" << n << " t:" << scheduler_.getTime() << std::endl;
      }
      ++n;
    }
  }
  const double last_time(scheduler_.getTime());
  if(!parallel_environment_.getrank() && verbose) {
    std::cout << "n:" << n << " t:" << last_time << std::endl;
  }
  //compartment_.outputCoordinates(last_time);
  compartment_.outputNumbers(last_time);
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

