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
#include "Event.hpp"

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
  Compartment compartment_;
  Species* out_species_;
  EventScheduler<SpatiocyteEvent> scheduler_;
  std::vector<Reaction*> independent_reactions_;
  std::vector<Reaction*> influenced_reactions_;
  int independent_event_index_ = -1;
  unsigned independent_event_id_;
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
