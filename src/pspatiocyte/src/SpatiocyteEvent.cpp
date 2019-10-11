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
                                 Species& species):
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

SpatiocyteEvent::SpatiocyteEvent(Lattice& lattice, Compartment& compartment,
                                 ParallelEnvironment& parallel_environment):
  lattice_(lattice),
  compartment_(compartment),
  parallel_environment_(parallel_environment),
  name_("IndependentReaction"),
  type_(INDEPENDENT_REACTION),
  interval_(compartment.get_next_time(parallel_environment, 0)) {
    set_time(interval_);
}

SpatiocyteEvent::SpatiocyteEvent(Lattice& lattice, Compartment& compartment,
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

void SpatiocyteEvent::fire() {
  switch(type_) {
  case DIFFUSION:
    compartment_.walk(species_list_, lattice_, parallel_environment_,
                      get_time());
    set_time(get_time() + interval_);
    return;
  case INDEPENDENT_REACTION:
    set_time(compartment_.react_direct_method(lattice_,
                                        parallel_environment_, get_time()));
    return;
  case OUTPUT_NUMBERS:
    interval_ = compartment_.output_numbers(get_time());
    set_time(get_time() + interval_);
    return;
  case OUTPUT_COORDINATES:
    interval_ = compartment_.output_coordinates(get_time());
    set_time(get_time() + interval_);
  }
}
