#ifndef __COMMON_HPP
#define __COMMON_HPP

#include <fstream>
#include <boost/random.hpp>
using namespace std;
using namespace boost;

//===== print control in general
//#define DIAGNOSTICS          //***

//===== print control in compartment.cpp
//#define PRINT_ALL_MOLECULES  //***
//#define SHOW_SUBVOLUME
//#define CHECK_INSIDE
//#define CONFIRM_LOCATION

//===== print control in lattice.cpp
//#define PRINT_DETAIL         //***

//===== print control in parallelenvironment.cpp
#define PRINT_NODETIMES

extern ofstream fout, fout2, fout3;          // for parallel files

/*
 *  random number generator of double precision between 0 and 1.
 *
 *  usage: Dice< uniform_real<>, double > random(seed);
 */
template <class X, class Y>
class Dice
{
public:
   Dice(uint32_t seed) : urand(mt19937(seed), X()) {}

   Y operator()() { return urand(); }

private:
   variate_generator< mt19937, X > urand;
};

// species type
enum SPECIES_TYPE
{
    DIFFUSIVE = 1,
    IMMOBILE = 2,
    HD = 3
};

// compartment type
enum COMPARTMENT_TYPE
{
    VOLUME = 1,
    SURFACE = 2,
};

// event type
enum EVENT_TYPE
{
    DIFFUSION = 1,
    INDEPENDENT_REACTION =2,
    INFLUENCED_REACTION = 3,
    COLLISION = 4,
    REACTION = 5
};

// reaction type
enum REACTION_TYPE
{
    INDEPENDENT = 1,
    INFLUENCED = 2
};

// cartesian coordinates
struct Coordinate {
    int ci;        // x component
    int cj;        // y component
    int ck;        // z component
};

// voxel
struct Voxel {
  int species_id;             // species-ID
  int mol_index;             // molecular-ID
};

class Species;
class World;
class Compartment;

// ghost size
#define GHOST_SIZE  1
#endif /* __COMMON_HPP */
