#ifndef __SPECIES_HPP
#define __SPECIES_HPP

#include <iostream>
#include <string>
#include <limits>
#include "Common.hpp"
using namespace std;

class Species {
public:
  Species() {}
  Species(string n, double D, unsigned init_size, World& world, double P=1);
  Species(const Species &s);
  ~Species() {} 
  string getName() {
    return name_;
  } 
  SPECIES_TYPE getType() { 
    return type_; 
  } 
  double getD() {
    return D_;
  }
  double getDt() {
    return dt_*walk_probability_;
  } 
  int getID() {
    return species_id_;
  } 
  void setType(SPECIES_TYPE type) { 
    type_ = type; 
  } 
  void setVolumeDt(double rv) { 
    const double alpha = 2.0/3.0;
    dt_ = ( type_==DIFFUSIVE ) ? alpha*rv*rv/D_ :
      std::numeric_limits<double>::infinity();
  } 
  void setSurfaceDt(double rv) {
    const double SQR2  = 1.414213562373;  
    const double SQR3  = 1.732050807569; 
    const double SQR6  = 2.449489742783;
    const double SQR22 = 4.690415759823;
    const double sqra = ( 2*SQR2 + 4*SQR3 + 3*SQR6 + SQR22 )/
                        ( 6*SQR2 + 4*SQR3 + 3*SQR6 );
    dt_ = ( type_==DIFFUSIVE ) ? sqra*sqra*rv*rv/D_ :
      std::numeric_limits<double>::infinity(); 
  } 
  void set_rho(double rho) {
    rho_ = rho;
  } 
  void calcCollisionTime() { 
    if (rho_ > P_) {
      walk_probability_ = P_/rho_;
    }
  } 
  float get_walk_probability() {
    return walk_probability_;
  } 
  bool operator==(const Species &s) const {
    return name_==s.name_ &&
      type_==s.type_ &&
      species_id_==s.species_id_;
  } 
  unsigned get_init_size() {
    return init_size_;
  }
  void diagnostics();

private:
  string name_;        // name of species
  double D_ = 0;     // diffusion coefficient
  double P_ = 1; // user-defined upper limit of reaction probability
  unsigned init_size_ = 0;
  double walk_probability_ = 1;
  SPECIES_TYPE type_;  // diffusive, immobile or HD
  double dt_;          // diffusion interval time
  double rho_;         // max reaction probability [max(p_j)] of all reactions
  int species_id_;
};

#endif /* __SPECIES_HPP */
