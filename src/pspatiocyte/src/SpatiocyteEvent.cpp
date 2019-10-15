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

#include "SpatiocyteEvent.hpp"
#include "Compartment.hpp"

SpatiocyteEvent::SpatiocyteEvent(Lattice& lattice, Compartment& compartment,
                                 ParallelEnvironment& parallel_environment,
                                 EventScheduler<SpatiocyteEvent>& scheduler,
                                 Species& species):
  lattice_(lattice),
  compartment_(compartment),
  parallel_environment_(parallel_environment),
  scheduler_(scheduler),
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

SpatiocyteEvent::SpatiocyteEvent(Lattice& lattice, Compartment& compartment,
                                 ParallelEnvironment& parallel_environment,
                                 EventScheduler<SpatiocyteEvent>& scheduler):
  lattice_(lattice),
  compartment_(compartment),
  parallel_environment_(parallel_environment),
  scheduler_(scheduler),
  name_("IndependentReaction"),
  type_(INDEPENDENT_REACTION),
  interval_(compartment.get_new_interval(parallel_environment)) {
    set_time(interval_);
}

SpatiocyteEvent::SpatiocyteEvent(Lattice& lattice, Compartment& compartment,
                                ParallelEnvironment& parallel_environment, 
                                EventScheduler<SpatiocyteEvent>& scheduler,
                                const EVENT_TYPE type, const double interval):
  lattice_(lattice),
  compartment_(compartment),
  parallel_environment_(parallel_environment),
  scheduler_(scheduler),
  name_("Logger"),
  type_(type),
  interval_(interval) {
    set_time(interval_);
    set_priority(-5);
}

//for direct method event:
void SpatiocyteEvent::update_next_time() {
  double next_time(scheduler_.get_time());
  const double old_next_time(get_time());
  if (old_next_time < std::numeric_limits<double>::infinity()) {
    next_time += compartment_.get_next_interval(parallel_environment_,
                                              old_next_time-next_time);
  }
  else {
    next_time += compartment_.get_new_interval(parallel_environment_);
  }
  set_time(next_time);
  if(next_time >= old_next_time) {
    scheduler_.get_queue().move_down(get_id());
  }
  else {
    scheduler_.get_queue().move_up(get_id());
  }
}

void SpatiocyteEvent::fire() {
  const double time(get_time());
  switch(type_) {
  case DIFFUSION:
    parallel_environment_.getcart().Barrier();
    compartment_.walk(species_list_, lattice_, parallel_environment_);
    parallel_environment_.getcart().Barrier();
    //If independent reaction event exists, update its time and reschedule
    //it in the execution queue:
    if (direct_method_event_) {
      direct_method_event_->update_next_time();
    }
    set_time(time + interval_);
    return;
  case INDEPENDENT_REACTION:
    set_time(time + compartment_.react_direct_method(lattice_,
                                                     parallel_environment_));
    return;
  case OUTPUT_NUMBERS:
    parallel_environment_.getcart().Barrier();
    set_time(time + compartment_.output_numbers(time));
    parallel_environment_.getcart().Barrier();
    return;
  case OUTPUT_COORDINATES:
    parallel_environment_.getcart().Barrier();
    set_time(time + compartment_.output_coordinates(time));
    parallel_environment_.getcart().Barrier();
  }
}
