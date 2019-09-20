#ifndef EVENTBASE_HPP
#define EVENTBASE_HPP

/**
   EventBase
   
   A subclass must define three customization points;

   void fire()
   {
   (1) do what this event is supposed to do.
   (2) setTime( next scheduled time of this event );
   }

   void update( const Event& anEvent )
   {
   Given the last fired Event (anEvent) that this Event
   depends on,

   (1) recalculate scheduled time (if necessary).
   (2) setTime( new scheduled time ).
   }

   const bool isDependentOn( const Event& anEvent )
   {
   Return true if this Event must be updated when the
   given Event (anEvent) fired.  Otherwise return false;
   }
*/

class EventBase
{
public:
    EventBase( const double& time )
        : time( time )
    {
        ; // do nothing
    }

    const double getTime() const
    {
        return this->time;
    }

    void setTime(const double& time)
    {
        this->time = time;
    }

    const bool operator<= ( const EventBase& rhs ) const
    {
        return getTime() <= rhs.getTime();
    }

    const bool operator< ( const EventBase& rhs ) const
    {
        return getTime() < rhs.getTime();
    }

    const bool operator== ( const EventBase& rhs ) const
    {
        return getTime() == rhs.getTime();
    }

    const bool operator!= ( const EventBase& rhs ) const
    {
        return ! this->operator==( rhs );
    }

    // dummy, because DynamicPriorityQueue requires this. better without.
    EventBase(): time( -1.0 )
    {
        ; // do nothing
    }
private:
    double  time;
};

#endif /* EVENTBASE_HPP */
