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
