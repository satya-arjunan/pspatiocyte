#include "World.hpp"
#include "Species.hpp"
#include "Reaction.hpp"

int main( int argc, char *argv[] ) {
  std::string dirname("output");
  double D(4);
  double ratio(1);

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
  const double duration(10); // s
  const double ka1(0.04483455086786913);
  const double kd1(1.35);
  const double kcat1(1.5);
  const double ka2(0.09299017957780264);
  const double kd2(1.73);
  const double kcat2(15.0);
  const double trel(1e-6);
  const double k7(log(2.)/trel);
  const int NKT(120*volume); // total K
  const int NP(rint(60*volume/(ratio+1)));
  const int NKK(60*volume-NP);
  const double dt(rv*rv*2/(D*3)); // s (diffusion interval)
  const int nlogs(duration*2);
  const double log_interval(duration/nlogs);
  //std::cout << "NKT:" << NKT << " NP:" << NP << " NKK:" << NKK << std::endl;

  World world(argc, argv, nx, ny, nz, rv, log_interval, dirname);

  Species KK("KK", D, NKK, world);
  Species Kpp("Kpp", D, 0, world);
  Species P("P", D, NP, world);
  Species K("K", D, NKT, world);
  Species Kp("Kp", D, 0, world);
  Species K_KK("K_KK", D, 0, world); //KK-K in Takahashi et al. PNAS 2010
  Species Kp_KK("Kp_KK", D, 0, world); //KK-Kp
  Species Kpp_P("Kpp_P", D, 0, world); //P-Kpp
  Species Kp_P("Kp_P", D, 0, world); //P-Kp
  Species KKa("KKa", D, 0, world);
  Species Pa("Pa", D, 0, world);

  Reaction r1("K + KK -> K_KK", K, KK, ka1, K_KK, world);
  Reaction r2("K_KK -> K + KK", K_KK, kd1, K, KK, world);
  Reaction r3("K_KK -> Kp + KKa", K_KK, kcat1, Kp, KKa, world);

  Reaction r4("Kp + KK -> Kp_KK", Kp, KK, ka2, Kp_KK, world);
  Reaction r5("Kp_KK -> Kp + KK", Kp_KK, kd2, Kp, KK, world);
  Reaction r6("Kp_KK -> Kpp + KKa", Kp_KK, kcat2, Kpp, KKa, world);

  Reaction r7("KKa -> KK", KKa, k7, KK, world);

  Reaction r8("Kpp + P -> Kpp_P", Kpp, P, ka1, Kpp_P, world);
  Reaction r9("Kpp_P -> Kpp + P", Kpp_P, kd1, Kpp, P, world);
  Reaction r10("Kpp_P -> Kp + Pa", Kpp_P, kcat1, Kp, Pa, world);

  Reaction r11("Kp + P -> Kp_P", Kp, P, ka2, Kp_P, world);
  Reaction r12("Kp_P -> Kp + P", Kp_P, kd2, Kp, P, world);
  Reaction r13("Kp_P -> K + Pa", Kp_P, kcat2, K, Pa, world);

  Reaction r14("Pa -> P", Pa, k7, P, world);

  world.initialize();
  world.run(duration, verbose);
  return 0;
}


