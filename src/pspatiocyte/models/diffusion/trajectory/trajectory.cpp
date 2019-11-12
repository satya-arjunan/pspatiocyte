#include <time.h>
#include "World.hpp"
#include "Species.hpp"
#include "Reaction.hpp"

int main( int argc, char *argv[] ) {
  std::string dirname("output");
  double D(0.06); // um2/s
  const double molecule_radius(0.0025); // um
  const int nx(476);
  const int ny(476);
  const int nz(476);
  const double rv(1.0208582*molecule_radius); // see Chew et al. PRE 2018.
  double volume(rv*rv*rv*(nx-2)*(ny-2)*(nz-2)*2.0*sqrt(3.0)*
                      sqrt(8.0/3.0));
  const bool verbose(true);
  const double duration(10); // s
  const int NA(1);
  const int NB(3);
  const int NC(1);
  const int nlogs(1000);
  const unsigned seed(time(NULL));
  const double numbers_log_interval(0); //don't log numbers
  const double coords_log_interval(duration/nlogs);
  
  //set log_interval to record trajectory (coordinates and id of molecules):
  World world(argc, argv, nx, ny, nz, rv, numbers_log_interval, dirname, seed,
              coords_log_interval);

  Species A("A", D, NA, world);
  Species B("B", D, NB, world);
  Species C("C", D, NC, world);

  world.initialize();
  world.run(duration, verbose);
  return 0;
}


