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

#ifndef __EventScheduler_hpp
#define __EventScheduler_hpp

#include "PriorityQueue.hpp"
#include "ParallelEnvironment.hpp"
#include "Common.hpp"

template <class Event>
class EventScheduler { 
public: 
  typedef PriorityQueue<Event*> EventPriorityQueue; 
  typedef typename EventPriorityQueue::ID EventID; 
  EventScheduler() {
    queue_.clear();
  }
  ~EventScheduler() {} 
  double get_time() {
    return time_;
  } 
  double get_top_time() {
    return get_top_event().get_time();
  } 
  EventID get_top_id() const {
    return queue_.get_top_id();
  } 
  Event& get_top_event() {
    return *queue_.get_top();
  } 

  unsigned add_event(Event event) {
    for (unsigned i(0); i < events_.size(); ++i) {
      if (event.get_event_type() == DIFFUSION &&
          events_[i].get_event_type() == DIFFUSION &&
          event.get_interval() == events_[i].get_interval()) {
        events_[i].add_species(event.get_species());
        return i;
      }
    }
    events_.push_back(event);
    return events_.size()-1;
  }

  void queue_events() {
    for (unsigned i(0); i < events_.size(); ++i) {
      Event& event(events_[i]);
      event.set_id(queue_.push(&event));
    }
  }

  const std::vector<Event>& get_events() const {
    return events_;
  }

  void update_event_time(const unsigned index, const double next_time) {
    Event& event(events_[index]);
    const double old_time(event.get_time());
    event.set_time(next_time);
    if(next_time >= old_time) {
      queue_.move_down(event.get_id());
    }
    else {
      queue_.move_up(event.get_id());
    }
  }

  void step() {
    time_ = queue_.get_top()->get_time();
    queue_.get_top()->fire();
    queue_.move_top();
  }

  /*
  Event& get_event(const EventID id) {
    return *queue_.get(id);
  } 

  void step(ParallelEnvironment &pe) {
    long long unsigned int idw[2];
    if(!pe.getrank()) {
      Event& topEvent(get_top_event());
      EventID ID(queue_.get_top_id());
      time_ = topEvent.get_time();
      MPI::Cartcomm cart = pe.getcart();
      idw[0] = ID;                                
      cart.Bcast( idw, 2, MPI::LONG_LONG, 0 );   
      ID = idw[0];                              
      topEvent.fire();
      queue_.move_top();
    } else {
      EventID ID;
      MPI::Cartcomm cart = pe.getcart();
      idw[0] = ID;                                 
      cart.Bcast( idw, 2, MPI::LONG_LONG, 0 );    
      ID = idw[0];                               
      Event& event(get_event(ID));
      time_ = event.get_time();
      event.fire();
      queue_.move_down(ID);
    }
  }
  */
private:
  double time_ = 0;
  EventPriorityQueue queue_;
  std::vector<Event> events_;
};

#endif /* __EventScheduler_hpp */
