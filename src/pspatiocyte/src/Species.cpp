#include "Species.hpp"
#include "World.hpp"

Species::Species(string n, double D,  unsigned init_size,
                 World& world, double P):
  name_(n),
  D_(D),
  P_(std::min(P, 1.0)),
  init_size_(init_size),
  dt_(std::numeric_limits<double>::infinity()),
  species_id_(world.add_species(this)) {}

Species::Species(const Species &s) {
  name_ = s.name_;
  D_ = s.D_;
  dt_ = s.dt_;
  P_ = s.P_;
  rho_ = s.rho_;
  walk_probability_ = s.walk_probability_;
  species_id_ = s.species_id_;
}

void Species::diagnostics() {
  fout << "Species:" << endl;
  fout << "   name = " <<  name_ << endl;
  fout << "   D    = " <<  D_    << endl;
  fout << "   dt   = " <<  dt_   << endl;
  fout << "   ID   = " <<  species_id_   << endl;
}
