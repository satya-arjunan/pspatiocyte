#ifndef __EVENTSCHEDULER_HPP
#define __EVENTSCHEDULER_HPP

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
    const float curr_time(event.get_time());
    event.set_time(next_time);
    if(next_time < curr_time) {
      queue_.move_up(event.get_id());
    }
    else {
      queue_.move_down(event.get_id());
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
      topEvent.fire(pe);
      queue_.move_top();
    } else {
      EventID ID;
      MPI::Cartcomm cart = pe.getcart();
      idw[0] = ID;                                 
      cart.Bcast( idw, 2, MPI::LONG_LONG, 0 );    
      ID = idw[0];                               
      Event& event(get_event(ID));
      time_ = event.get_time();
      event.fire(pe);
      queue_.move_down(ID);
    }
  }
  */
private:
  double time_ = 0;
  EventPriorityQueue queue_;
  std::vector<Event> events_;
};

#endif /* __EVENTSCHEDULER_HPP */
