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

#ifndef __SpatiocyteEvent_HPP
#define __SpatiocyteEvent_HPP

#include <iostream>
#include <limits>
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
    name_(std::string("Diffusion:")+species.get_name()),
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
    interval_(std::numeric_limits<double>::infinity()) {
      set_time(interval_);
  }

  SpatiocyteEvent(Lattice& lattice, Compartment& compartment,
                  ParallelEnvironment& parallel_environment,
                  const EVENT_TYPE type, const double interval):
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
  
  double get_interval() const {
    return interval_;
  } 
  
  Species& get_species() const {
    return *species_;
  } 
  
  void add_species(Species& species) {
    species_list_.push_back(&species);
  }
  void fire();

private:
  Lattice& lattice_;
  Compartment& compartment_;
  ParallelEnvironment& parallel_environment_;
  Species* species_;
  std::string name_;
  EVENT_TYPE type_;
  double interval_;
  std::vector<Species*> species_list_;
  int id_;
};

#endif /* __SpatiocyteEvent_HPP */
