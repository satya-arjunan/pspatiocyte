#ifndef __SPATIOEVENT_HPP
#define __SPATIOEVENT_HPP

#include <iostream>
#include <limits>
#include "Compartment.hpp"
#include "Lattice.hpp"
#include "Species.hpp"
#include "Reaction.hpp"
#include "Common.hpp"
#include "EventScheduler.hpp"
#include "EventBase.hpp"

using namespace std;

struct SpatiocyteEvent: public EventBase {
  SpatiocyteEvent(Lattice* latt, Compartment* comp, Species* spec) {
    etype_ = DIFFUSION;
    dt_ = spec->getDt();
    pspec_ = spec;
    species_list_.push_back(spec);
    pcomp_ = comp;
    platt_ = latt;
    name_ = spec->getName();
    setTime(dt_);
  }
  
  SpatiocyteEvent(Lattice* latt, Compartment* comp, ParallelEnvironment* pe) {
    etype_ = INDEPENDENT_REACTION;
    dt_ = std::numeric_limits<double>::infinity();
    pcomp_ = comp;
    platt_ = latt;
    setTime(dt_);
  } 

  SpatiocyteEvent(Lattice* latt, Compartment* comp, ParallelEnvironment* pe,
                  EVENT_TYPE etype, double dt) {
    etype_ = etype;
    dt_ = dt;
    pcomp_ = comp;
    platt_ = latt;
    setTime(dt_);
  } 

  int getID() {
    return id_;
  } 
  
  void setID(int id) {
    id_ = id;
  }
  
  const string getName() const {
    return name_;
  }
  
  EVENT_TYPE get_event_type() {
    return etype_;
  }
  
  double get_interval() const {
    return dt_;
  } 
  
  Species* get_species() const {
    return pspec_;
  } 
  
  void add_species(Species* species) {
    species_list_.push_back(species);
  } 
  
  void fire(ParallelEnvironment &pe) {
    switch(etype_) {
    case DIFFUSION:
      pcomp_->walk(species_list_, *platt_, pe);
      break;
    case INDEPENDENT_REACTION:
      dt_ = pcomp_->react_direct_method(*platt_, pe, getTime());
      break;
    case OUTPUT_NUMBERS:
      pcomp_->output_numbers(getTime());
      break;
    case OUTPUT_COORDINATES:
      pcomp_->output_coordinates(getTime());
      break;
    }
    setTime(getTime() + dt_);
  }

private:
  int id_;
  EVENT_TYPE etype_;
  double dt_;
  Species* pspec_;
  Compartment* pcomp_;
  Lattice* platt_;
  string name_;
  std::vector<Species*> species_list_;
};

#endif /* __SPATIOEVENT_HPP */
