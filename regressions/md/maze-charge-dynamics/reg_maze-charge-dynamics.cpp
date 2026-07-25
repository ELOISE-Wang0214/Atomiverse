
#include <cmath>
#include <cstdio>
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
using namespace ATV;

#include "atomiverse/cpu/dynamics/const_pot_maze.h"
int main() {
  cpu::Frame f(4);
  double zz[4]={8,12,9,11}; int ty[4]={1,2,3,3}; double qq[4]={0,0,0.5,-0.5};
  for(int i=0;i<4;++i){
    f.span(kAtom_x)[i]=f.span(kAtom_y)[i]=0; f.span(kAtom_z)[i]=zz[i];
    f.span(kAtom_type)[i]=ty[i]; f.span(kAtom_charge)[i]=qq[i];
    f.span(kAtom_mass)[i]=1; f.span(kAtom_tag)[i]=i+1;
  }
  f.simbox=cpu::Box::orthorhombic(1,1,20);
  cpu::md::ConstPotMaZe maze(2.0);
  maze.add_electrode(1,0.0); maze.add_electrode(2,1.0);
  maze.solver().xy_image_range=50;
  cpu::Context ctx; maze.build({ctx,f});
  maze.launch({ctx,f,0,0});
  double q0=maze.solver().total_charge(1);
  f.span(kAtom_z)[2]+=0.1;
  maze.launch({ctx,f,0,1});
  double dQ=maze.max_dQ();
  double q1=maze.solver().total_charge(1);
  std::printf("maze-charge-dynamics Q0=%.6e Q1=%.6e max_dQ=%.3e\n", q0,q1,dQ);
  int fail=0;
  if (!(dQ>=0)) ++fail;
  if (maze.solver().neutrality_residual()>1e-8) ++fail;
  return fail?1:0;
}
