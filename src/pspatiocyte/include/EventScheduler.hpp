#ifndef __EVENTSCHEDULER_HPP
#define __EVENTSCHEDULER_HPP

#include "DynamicPriorityQueue.hpp"
#include "ParallelEnvironment.hpp"
#include "Common.hpp"

//#undef DIAGNOSTICS

/**
   Event scheduler.

   This class works as a sequential
   event scheduler with a heap-tree based priority queue.

*/

template <class Event_>
class EventScheduler
{
  
public:

    typedef Event_ Event;
    typedef DynamicPriorityQueue<Event> EventPriorityQueue;

    typedef typename DynamicPriorityQueue<Event>::Index EventIndex;
    typedef typename DynamicPriorityQueue<Event>::ID EventID;

    //    typedef std::vector<EventIndex> EventIndexVector;
    //    typedef std::vector<EventIndexVector> EventIndexVectorVector;


    EventScheduler(): time( 0.0 )
    {
        ; // do nothing
    }

    ~EventScheduler()
    {
        ; // do nothing
    }


    const double getTime() const
    {
        return time;
    }

    const double getTopTime() const
    {
        assert( this->getSize() != 0 );  // FIXME: use exception

        return getTopEvent().getTime();
    }

    const EventIndex getSize() const
    {
        return this->eventPriorityQueue.getSize();
    }

    const Event& getTopEvent() const
    {
        return this->eventPriorityQueue.getTop();
    }

    EventID getTopID() const
    {
        return this->eventPriorityQueue.getTopID();
    }

    const Event& getEvent( const EventID id ) const
    {
        return this->eventPriorityQueue.get( id );
    }

    const Event& getEventByIndex( const EventIndex index ) const
    {
        return this->eventPriorityQueue.getByIndex( index );
    }

    void step_async(ParallelEnvironment &pe)
    {
        this->time = this->eventPriorityQueue.getTop().getTime();
        this->eventPriorityQueue.getTop().fire(pe);
        this->eventPriorityQueue.moveTop();
    }

    /*
    void step(ParallelEnvironment &pe)
    {
        MPI::Cartcomm cart = pe.getcart();
        cart.Barrier();
        this->time = this->eventPriorityQueue.getTop().getTime();
        this->eventPriorityQueue.getTop().fire(pe);
        this->eventPriorityQueue.moveTop();
    }
    /*

    void step(ParallelEnvironment &pe)
    {
        MPI::Cartcomm cart = pe.getcart();
        cart.Barrier();
        Event topEvent( getTopEvent() );
        EventID ID( this->eventPriorityQueue.getTopID() );
        this->time = topEvent.getTime();
        topEvent.fire(pe);
        if( topEvent.getTime() >= getTime() )
        {
            this->eventPriorityQueue.replace( ID, topEvent );
        }
        else
        {
            this->eventPriorityQueue.pop( ID );
        }
    }*/


    void step(ParallelEnvironment &pe)
    {
        // master process
        long long unsigned int idw[2];                      // temporal tuning for K
        if( pe.getrank()==0 )
        {
            // Here I copy construct the top event and use its event
            // ID to reschedule it.  This is necessary if events can
            // be created or deleted within fire() and the dynamic
            // priority queue can reallocate internal data structures.
            // Most of the cost of using this is optimized away when
            // the dynamic priority queue has a VolatileIDPolicy.
            // 
            // const EventID ID( this->eventPriorityQueue.getTopID() );
            Event topEvent( getTopEvent() );
            EventID ID( this->eventPriorityQueue.getTopID() );
            this->time = topEvent.getTime();
    
            // broadcast eventID 
            MPI::Cartcomm cart = pe.getcart();
            // cart.Bcast( &ID, 1, MPI::LONG_LONG, 0 );     // temporal tuning for K
            idw[0] = ID;                                    //
            cart.Bcast( idw, 2, MPI::LONG_LONG, 0 );        //
            ID = idw[0];                                    //

            #ifdef DIAGNOSTICS
            // fout << "DEBUG: EventID = " << ID << endl;
            #endif
    
            // Fire top
            topEvent.fire(pe);
    
            // If the event is rescheduled into the past, remove it.
            // Otherwise, reuse the event.
            if( topEvent.getTime() >= getTime() )
            {
                this->eventPriorityQueue.replace( ID, topEvent );
            }
            else
            {
                this->eventPriorityQueue.pop( ID );
            }
        // worker processes
        } else {
            // broadcast eventID 
            EventID ID;
            MPI::Cartcomm cart = pe.getcart();
            // cart.Bcast( &ID, 1, MPI::LONG_LONG, 0 );     // temporal tuning for K
            idw[0] = ID;                                  //
            cart.Bcast( idw, 2, MPI::LONG_LONG, 0 );        //
            ID = idw[0];                                    //

            #ifdef DIAGNOSTICS
            // fout << "DEBUG: EventID = " << ID << endl;
            #endif

            // identify event to be fired
            Event topEvent = getEvent(ID);
            this->time = topEvent.getTime();

            // Fire identified event
            topEvent.fire(pe);

            // Actually, PriorityQue of slave processes is not used,
            // exists only to hold SpatiocytoEvent.
            // However following code fragments are necessary to get
            // slave-timers to be consistent to master-timer.
            if( topEvent.getTime() >= getTime() )
            {
                this->eventPriorityQueue.replace( ID, topEvent );
            }
            else
            {
                this->eventPriorityQueue.pop( ID );
            }
        }
    }

    void clear()
    {
        time = 0.0;
        this->eventPriorityQueue.clear();
    }


    /*
    EventID addEvent( const Event& event )
    {
        EventID id(this->eventPriorityQueue.push( event ));
        Event& newEvent(this->eventPriorityQueue.get(id));
        newEvent.setID(id);
        return id;
    }
    */

    EventID addEvent(Event event) {
      for (unsigned i(0); i < events_.size(); ++i) {
        if (event.get_event_type() == DIFFUSION &&
            event.get_event_type() == events_[i]->get_event_type() &&
            event.get_interval() == events_[i]->get_interval()) {
          Event& queued_event(
              this->eventPriorityQueue.get(events_[i]->getID()));
          queued_event.add_species(event.get_species());
          fout << std::endl;
          if (event.get_event_type() == DIFFUSION) {
            fout << "diffusion: " << events_[i]->get_interval() << 
              std::endl;
          }
          fout << "adding:" << event.get_species()->getName() <<
            " parent:" << queued_event.get_species()->getName() << std::endl;
          return 0;
        }
      }
      EventID id(this->eventPriorityQueue.push(event));
      Event& newEvent(this->eventPriorityQueue.get(id));
      newEvent.setID(id);
      events_.push_back(&newEvent);
      return id;
    }

    void removeEvent( const EventID id )
    {
        this->eventPriorityQueue.pop( id );
    }


    void updateEventTime( const EventID id, const double t )
    {
        const EventIndex index( this->eventPriorityQueue.getIndex( id ) );
        Event& event( this->eventPriorityQueue.getByIndex( index ) );

        event.setTime( t );
        this->eventPriorityQueue.move( index );
    }

    const bool check() const
    {
        return this->eventPriorityQueue.check();
    }

private:
    EventPriorityQueue       eventPriorityQueue;
    double                   time;
    std::vector<Event*> events_;
};

#endif /* __EVENTSCHEDULER_HPP */
