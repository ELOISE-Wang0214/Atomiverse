
#include <cmath>
#include <cstdio>
#include "atomiverse/core/frame_schema.h"
#include "atomiverse/cpu/component/context.h"
#include "atomiverse/cpu/frame/frame.h"
using namespace ATV;

#include "atomiverse/cpu/dynamics/const_pot_bo.h"
#include "atomiverse/cpu/observe/elastance_observer.h"
int main() {
  cpu::Frame f(2);
  f.span(kAtom_z)[0]=5; f.span(kAtom_z)[1]=15;
  for(int i=0;i<2;++i){
    f.span(kAtom_x)[i]=f.span(kAtom_y)[i]=0; f.span(kAtom_type)[i]=i+1;
    f.span(kAtom_charge)[i]=0; f.span(kAtom_mass)[i]=1; f.span(kAtom_tag)[i]=i+1;
  }
  f.simbox=cpu::Box::orthorhombic(5,5,20);
  cpu::md::ConstPotBO bo(2.0);
  bo.add_electrode(1,0); bo.add_electrode(2,1); bo.xy_image_range=80;
  cpu::ElastanceObserver obs(&bo);
  cpu::Context ctx; bo.build({ctx,f}); bo.launch({ctx,f,0,0}); obs.launch({ctx,f,0,0});
  std::printf("elastance-condition-observer cond=%.6e C=%.6e\n",
         obs.last_condition(), obs.last_capacitance());
  return (obs.last_condition()>=1.0 && std::isfinite(obs.last_capacitance())
          && obs.last_capacitance()>0.0)?0:1;
}
