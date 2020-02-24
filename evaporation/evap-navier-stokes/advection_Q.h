#include "run.h"
#include "timestep.h"

face vector uf[];
#define u uf

double (* gradient) (double, double, double) = NULL;

extern scalar * tracers;

event defaults (i = 0) {
  for (scalar f in tracers)
    f.gradient = gradient;
}

event init (i = 0) {
  boundary ((scalar *){u});
  boundary (tracers);
}

#define velocity stability

event stability (i++, last) {
  dt = dtnext (timestep (u, DT));
}

event vof (i++, last);
event tracer_advection (i++, last);
event tracer_diffusion (i++, last);

#if TREE
event adapt (i++,last);
#endif
#include "tracer.h"


