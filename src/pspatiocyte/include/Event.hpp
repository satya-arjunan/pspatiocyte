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

struct SpatiocyteEvent: public EventBase {
  SpatiocyteEvent(Lattice& lattice, Compartment& compartment,
                  ParallelEnvironment& parallel_environment, Species& species):
    lattice_(lattice),
    compartment_(compartment),
    parallel_environment_(parallel_environment),
    species_(&species),
    name_(species.get_name()),
    type_(DIFFUSION),
    interval_(species.get_walk_interval()) {
      if (!parallel_environment_.getrank()) { 
        std::cout << "Species " << species.get_name() << 
          " walk interval:" << species.get_walk_interval() << std::endl;
      }
      species_list_.push_back(&species);
      set_time(interval_);
  }

  SpatiocyteEvent(Lattice& lattice, Compartment& compartment,
                  ParallelEnvironment& parallel_environment):
    lattice_(lattice),
    compartment_(compartment),
    parallel_environment_(parallel_environment),
    name_("IndependentReaction"),
    type_(INDEPENDENT_REACTION),
    interval_(std::numeric_limits<float>::infinity()) {
      set_time(interval_);
  }

  SpatiocyteEvent(Lattice& lattice, Compartment& compartment,
                  ParallelEnvironment& parallel_environment,
                  const EVENT_TYPE type, const float interval):
    lattice_(lattice),
    compartment_(compartment),
    parallel_environment_(parallel_environment),
    name_("Logger"),
    type_(type),
    interval_(interval) {
      set_time(interval_);
      set_priority(-5);
  }
  
  int get_id() const {
    return id_;
  } 
  
  void set_id(int id) {
    id_ = id;
  }
  
  const std::string get_name() const {
    return name_;
  }
  
  EVENT_TYPE get_event_type() {
    return type_;
  }
  
  float get_interval() const {
    return interval_;
  } 
  
  Species& get_species() const {
    return *species_;
  } 
  
  void add_species(Species& species) {
    species_list_.push_back(&species);
  }
  
  void fire() {
    switch(type_) {
    case DIFFUSION:
      compartment_.walk(species_list_, lattice_, parallel_environment_);
      break;
    case INDEPENDENT_REACTION:
      interval_ = compartment_.react_direct_method(lattice_,
                                           parallel_environment_, get_time());
      break;
    case OUTPUT_NUMBERS:
      interval_ = compartment_.output_numbers(get_time());
      break;
    case OUTPUT_COORDINATES:
      interval_ = compartment_.output_coordinates(get_time());
      break;
    }
    set_time(get_time() + interval_);
  }

private:
  Lattice& lattice_;
  Compartment& compartment_;
  ParallelEnvironment& parallel_environment_;
  Species* species_;
  std::string name_;
  EVENT_TYPE type_;
  float interval_;
  std::vector<Species*> species_list_;
  int id_;
};

#endif /* __SPATIOEVENT_HPP */
