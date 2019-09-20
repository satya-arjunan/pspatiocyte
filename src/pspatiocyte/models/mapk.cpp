#include "world.hpp"
#include "species.hpp"
#include "reaction.hpp"

int main( int argc, char *argv[] ) {
  const int nx(220); // lattice dimensions
  const int ny(220);
  const int nz(220);
  const double rv(0.0255215);
  const double D(0.06);
  const double volume(rv*rv*rv*nx*ny*nz*2.0*sqrt(3.0)*sqrt(8.0/3.0));
  const bool verbose(true);
  const double duration(20); // s
  const double ka1(0.04483455086786913);
  const double kd1(1.35);
  const double kcat1(1.5);
  const double ka2(0.09299017957780264);
  const double kd2(1.73);
  const double kcat2(15.0);
  const double trel(1e-6);
  const double k7(log(2.)/trel);
  const double ratio(1.3688745095370807); //ratios[6]
  const int NKT(120*volume); // total K
  const int NPP(rint(60*volume/(ratio+1)));
  const int NKK(60*volume-NPP);
  const double dt(rv*rv*2/(D*3)); // s (diffusion interval)
  std::cout << "dt:" << dt << "rv:" << rv << "volume:" << volume << std::endl;
  const int nlogs(100);
  const double log_interval(duration/nlogs);

  World world(argc, argv, nx, ny, nz, rv);

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


