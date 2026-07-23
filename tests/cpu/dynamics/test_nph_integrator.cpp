#include <gtest/gtest.h>
#include <cmath>
#include <cstdio>
#include "atomiverse/core/units.h"
#include "atomiverse/cpu/dynamics/init.h"
#include "atomiverse/cpu/dynamics/md.h"
#include "atomiverse/cpu/dynamics/nose_hoover_npt.h"
#include "atomiverse/cpu/dynamics/parrinello_rahman_npt.h"
#include "atomiverse/cpu/dynamics/nph_integrator.h"
#include "atomiverse/cpu/neighbor/neighbor_list.h"
#include "atomiverse/cpu/potentials/pair_lj_cut.h"
using namespace ATV;
static cpu::Frame mk() {
  int ns=4,n=ns*ns*ns; double L=cbrt(n/0.02071),sp=L/ns;
  cpu::Frame f(n); int i=0;
  for(int x=0;x<ns;x++) for(int y=0;y<ns;y++) for(int z=0;z<ns;z++){
    f.span(kAtom_x)[i]=(x+0.5)*sp; f.span(kAtom_y)[i]=(y+0.5)*sp; f.span(kAtom_z)[i]=(z+0.5)*sp;
    f.span(kAtom_mass)[i]=39.948; f.span(kAtom_type)[i]=0; i++;
  }
  f.simbox=cpu::Box::cubic(L); return f;
}
TEST(Diag,Long) {
  const double eps=120*units::kB_ev_per_k,sigma=3.4,cut=2.5*sigma,T=94.4,dt=2;
  auto run=[&](const char* tag, auto seti){
    auto f=mk(); double V0=f.simbox.volume();
    cpu::md::MD sim; sim.add_potential(cpu::PairLJCut{eps,sigma,cut});
    sim.set_neighbor(cpu::md::NeighborList<>{cut,0.5}); seti(sim);
    sim.add(cpu::md::InitVelocity{T,42u}); sim.enable_timing(false);
    sim.build();
    sim.run(f,200); // warm
    double H0=f.meta.get(kConserved); double sumrho=0; int ns=0;
    double worst=0;
    for(int s=0;s<20;s++){
      sim.run(f,20);
      double H=f.meta.get(kConserved);
      double rel=std::fabs(H-H0)/std::fabs(H0);
      if(rel>worst) worst=rel;
      sumrho+=64.0/f.simbox.volume(); ns++;
    }
    std::printf("%s V0=%.1f Vf=%.1f meanrho=%.5f worstHrel=%.3e H0=%.3f\n",
      tag,V0,f.simbox.volume(),sumrho/ns,worst,H0);
  };
  run("MTK",[](cpu::md::MD&s){s.set_integrator(cpu::md::NoseHooverNPT{{.T_target=94.4,.tau_T=50,.P_target=0,.tau_P=500,.dt=2}});});
  run("PR ",[](cpu::md::MD&s){s.set_integrator(cpu::md::ParrinelloRahmanNPT{94.4,0,50,500,2});});
  run("NPH",[](cpu::md::MD&s){s.set_integrator(cpu::md::NphIntegrator{0,500,2});});
  SUCCEED();
}
