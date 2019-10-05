#include "Species.hpp"
#include "World.hpp"

Species::Species(string n, double D,  unsigned init_size,
                 World& world, double P):
  name_(n),
  D_(D),
  P_(std::min(P, 1.0)),
  init_size_(init_size),
  walk_interval_(std::numeric_limits<double>::infinity()),
  species_id_(world.add_species(this)) {}

Species::Species(const Species &s) {
  name_ = s.name_;
  D_ = s.D_;
  walk_interval_ = s.walk_interval_;
  P_ = s.P_;
  rho_ = s.rho_;
  walk_probability_ = s.walk_probability_;
  species_id_ = s.species_id_;
}


