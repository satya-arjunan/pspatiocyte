#include "World.hpp"
#include "Species.hpp"
#include "Reaction.hpp"

int main(int argc, char *argv[]) {
  //for unforced search up to 75% crowding is ok 
  //for forced search up to 95% crowding is good
  std::string output_dirname("output");
  bool force_search_vacant(false);
  float crowd_fraction(0.75);
  if (argc >= 2) {
    output_dirname = argv[1];
  }
  if (argc >= 3) {
    force_search_vacant = bool(atoi(argv[2]));
  }
  if (argc >= 4) {
    crowd_fraction = atof(argv[3]);
  }

  const int nx(64); // lattice dimensions
  const int ny(64);
  const int nz(64);
  const bool verbose(true);
  const double duration(11); // s
  const double volume(909e-18); // m^3
  const double rv(pow(volume/(((nz-2)*2.0)*((ny-2)*sqrt(3.0))*
                              ((nx-2)*sqrt(8.0/3.0))), 1.0/3.0)); //voxel radius
  const int num_voxels((nx-2)*(ny-2)*(nz-2));
  const int num_A(50000);
  const int num_crowd(std::max(int(crowd_fraction*num_voxels-num_A*2), 0));
  const double D(10e-12); // m^2/s (diffusion coefficient)
  const int nlogs(20000);
  const double log_interval(duration/nlogs);
  const bool log_coordinates(false);
  World world(argc, argv, nx, ny, nz, rv, log_interval, output_dirname,
              0, force_search_vacant);
  Species A("A", D, num_A, world);
  Species B("B", D, 0, world);
  Species C("C", D, 0, world);
  Species E("E", D, num_crowd, world);

  const double k(1); //dissociation rate

  Reaction product("A -> B + C", A, k, B, C, world);

  world.initialize();
  world.run(duration, verbose);
  return 0;
}


