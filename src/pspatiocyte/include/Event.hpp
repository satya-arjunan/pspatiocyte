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

//#undef DIAGNOSTICS

struct SpatiocyteEvent: public EventBase
{
    // diffusion event
    SpatiocyteEvent(Lattice* latt, Compartment* comp, Species* spec)
    {
        etype_ = DIFFUSION;
        dt_ = spec->getDt();
        std::cout << spec->getName() << " diffusion interval:" << dt_ << 
          std::endl;
        pspec_ = spec;
        species_list_.push_back(spec);
        pcomp_ = comp;
        platt_ = latt;
        name_ = spec->getName();
        setTime(dt_);                 // dt as initial time in EventBase
    }

    // independent-reaction event
    SpatiocyteEvent(Lattice* latt, Compartment* comp, ParallelEnvironment* pe)
    {
        //fout << endl << "independent event" << std::endl;
        etype_ = INDEPENDENT_REACTION;
        dt_ = std::numeric_limits<double>::infinity();
        //fout << "init independent event dt:" << dt_ << std::endl;
        pcomp_ = comp;
        platt_ = latt;
        setTime(dt_);                 // dt as initial time in EventBase
    }

    int getID()
    {
        return id_;
    }

    void setID(int id)
    {
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

    void fire(ParallelEnvironment &pe)
    {
      /*
        fout << endl << "Fire!  time = "
             << setw(8) << setprecision(4) << showpoint << left
             << getTime();
        fout << " : dt = "   << dt_ ;
        fout << " : id = "   << id_ 
             << " : type = " << etype_ << std::endl;;
             */
        switch(etype_)
        {
            case DIFFUSION:
                 #ifdef DIAGNOSTICS
                 fout << " : diffusion "<< name_ << " id:" << id_ << endl; 
                 fout << "\t>>>> walk invoked <<<<" << endl;
                 #endif
                 pcomp_->walk(species_list_, *platt_, pe);
                 break;
            case INDEPENDENT_REACTION:
                 /*
                 fout << " : independent-reaction " << name_ << " id:" << id_ << endl;
                 fout << "\t>>>> independentReaction invoked <<<<" << endl;
                 */
                 dt_ = pcomp_->executeDirectMethodReaction(*platt_, pe, getTime());
                 break;
        }
        //fout << "fire current time:" << getTime() << std::endl;
        setTime( getTime() + dt_ );   // update time in EventBase
        //fout << "fire next time:" << getTime() << std::endl;
    }

private:
    int id_;
    EVENT_TYPE etype_;      // DIFFUSION, COLLISION, REACTION
    double dt_;
    Species* pspec_;
    Compartment* pcomp_;
    Lattice* platt_;
    string name_;
    std::vector<Species*> species_list_;
};

//theTime: next_react_time_
//

#endif /* __SPATIOEVENT_HPP */
