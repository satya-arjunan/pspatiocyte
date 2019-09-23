#include "World.hpp"
#include "Species.hpp"
#include "Reaction.hpp"

int main( int argc, char *argv[] ) {
  std::string dirname("output");
  double D(0.06);
  double ratio(0.11103363181676379); //ratios[2]

  if (argc >= 2) {
    dirname = argv[1];
  }
  if (argc >= 3) {
    D = atof(argv[2]);
  }
  if (argc >= 4) {
    ratio = atof(argv[3]);
  }

  const double molecule_radius(0.0025);
  const int nx(476); // lattice dimensions for volume = 10
  const int ny(476);
  const int nz(476);
  const double rv(1.0208582*molecule_radius); // see Chew et al. PRE 2018.
  const double volume(rv*rv*rv*(nx-2)*(ny-2)*(nz-2)*2.0*sqrt(3.0)*
                      sqrt(8.0/3.0));
  const bool verbose(true);
  const double duration(400); // s
  const double ka1(0.04483455086786913);
  const double kd1(1.35);
  const double kcat1(1.5);
  const double ka2(0.09299017957780264);
  const double kd2(1.73);
  const double kcat2(15.0);
  const double trel(1e-6);
  const double k7(log(2.)/trel);
  const int NKT(120*volume); // total K
  const int NPP(rint(60*volume/(ratio+1)));
  const int NKK(60*volume-NPP);
  const double dt(rv*rv*2/(D*3)); // s (diffusion interval)
  const int nlogs(100);
  const double log_interval(duration/nlogs);

  World world(argc, argv, nx, ny, nz, rv, dirname);

  Species KK("KK", D, NKK, world);
  Species Kpp("Kpp", D, 0, world);
  Species PP("PP", D, NPP, world);
  Species K("K", D, NKT, world);
  Species Kp("Kp", D, 0, world);
  Species K_KK("K_KK", D, 0, world);
  Species Kp_KK("Kp_KK", D, 0, world);
  Species Kpp_PP("Kpp_PP", D, 0, world);
  Species Kp_PP("Kp_PP", D, 0, world);
  Species KKa("KKa", D, 0, world);
  Species PPa("PPa", D, 0, world);

  Reaction r1("K + KK -> K_KK", K, KK, ka1, K_KK, world);
  Reaction r2("K_KK -> K + KK", K_KK, kd1, K, KK, world);
  Reaction r3("K_KK -> Kp + KKa", K_KK, kcat1, Kp, KKa, world);
  Reaction r4("Kp + KK -> Kp_KK", Kp, KK, ka2, Kp_KK, world);
  Reaction r5("Kp_KK -> Kp + KK", Kp_KK, kd2, Kp, KK, world);
  Reaction r6("Kp_KK -> Kpp + KKa", Kp_KK, kcat2, Kpp, KKa, world);
  Reaction r7("KKa -> KK", KKa, k7, KK, world);
  Reaction r8("Kpp + PP -> Kpp_PP", Kpp, PP, ka1, Kpp_PP, world);
  Reaction r9("Kpp_PP -> Kpp + PP", Kpp_PP, kd1, Kpp, PP, world);
  Reaction r10("Kpp_PP -> Kp + PPa", Kpp_PP, kcat1, Kp, PPa, world);
  Reaction r11("Kp + PP -> Kp_PP", Kp, PP, ka2, Kp_PP, world);
  Reaction r12("Kp_PP -> Kp + PP", Kp_PP, kd2, Kp, PP, world);
  Reaction r13("Kp_PP -> K + PPa", Kp_PP, kcat2, K, PPa, world);
  Reaction r14("PPa -> PP", PPa, k7, PP, world);

  world.initialize();
  world.run(log_interval, duration, verbose);
  return 0;
}


