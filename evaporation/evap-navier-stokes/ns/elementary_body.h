#ifndef BASILISK_HEADER_14
#define BASILISK_HEADER_14
#line 1 "./../../elementary_body.h"
#ifndef F_ERR
  #define F_ERR 1e-10
#endif

#include "vof.h"
#include "diffusion.h"
#include "curvature.h"
#include "my_functions.h"

attribute {
  double D;
  double tr_eq;
  double peclet;
}

void phase_change_velocity (scalar f, scalar tr, face vector v_pc) {

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  
  face vector gtr[];
  foreach_face()
    gtr.x[] = (tr[] - tr[-1])/Delta;
  boundary((scalar*){gtr});
  
  vector n[];
  compute_normal (f, n);

  foreach_face() {
    v_pc.x[] = 0.;
    if (interfacial(point, f) || interfacial(neighborp(-1), f)) {
      coord nf;
      foreach_dimension()
        nf.x = 0.;
      if (interfacial(point, f)) {
        foreach_dimension()
          nf.x += n.x[];
      }
      if (interfacial(neighborp(-1), f)) {
        nf.x += n.x[-1];
        nf.y += n.y[-1];
      }
      double norm = 0.;
      foreach_dimension()
        norm += fabs(nf.x);
      foreach_dimension()
        nf.x /= (tr.inverse ? norm : - norm);
      if (nf.x > 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[1, 0] 
                    + fabs(nf.y)*(nf.y > 0. ? gtr.x[1, 1] : gtr.x[1, -1]));
      }
      else if (nf.x < 0.) {
        v_pc.x[] = (fabs(nf.x)*gtr.x[-1, 0]
                    + fabs(nf.y)*(nf.y > 0. ? gtr.x[-1, 1] : gtr.x[-1, -1]));
      }
      v_pc.x[] *= fm.x[]*tr.D*tr.peclet;
    }
  }
  boundary((scalar *){v_pc});
}

struct Dirichlet_Diffusion {
  // mandatory
  scalar tr;
  scalar f;
  int max_level;
  double dt;
  double time_factor;
  // optional
  scalar tr_op; // default uniform 1.
};

mgstats dirichlet_diffusion (struct Dirichlet_Diffusion p) {
  
  scalar tr = p.tr, f = p.f, tr_op = automatic (p.tr_op);
  int max_level = p.max_level;
  double time_factor = p.time_factor;

  foreach()
    f[] = clamp(f[], 0., 1.);
  boundary ({f, tr});
  if (p.tr_op.i)
    boundary ({tr_op});

  scalar volumic_metric[], dirichlet_source_term[], dirichlet_feedback_term[];
  face vector diffusion_coef[];

  foreach() {
    volumic_metric[] = cm[];
    if (p.tr_op.i)  
      dirichlet_source_term[] = cm[]*tr.tr_eq*tr_op[]*tr.D*time_factor*(tr.inverse ?
                                  f[] : 1. - f[])*sq((1<<max_level)/L0);
    else
      dirichlet_source_term[] = cm[]*tr.tr_eq*tr.D*time_factor*(tr.inverse ?
                                  f[] : 1. - f[])*sq((1<<max_level)/L0);
    dirichlet_feedback_term[] = - cm[]*tr.D*time_factor*(tr.inverse ?
                                  f[] : 1. - f[])*sq((1<<max_level)/L0);
  }
  foreach_face()
    diffusion_coef.x[] = fm.x[]*tr.D;

  boundary({volumic_metric, dirichlet_source_term, dirichlet_feedback_term, diffusion_coef});

  return diffusion (tr, dt, D = diffusion_coef, r = dirichlet_source_term, beta = dirichlet_feedback_term, theta = volumic_metric);
  
}



#endif
