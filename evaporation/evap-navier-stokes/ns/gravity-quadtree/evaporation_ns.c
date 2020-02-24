/* Compile with:
 * 
 * 3D with MPI and graphics: CC99='mpicc -std=c99' qcc -Wall -O2 -grid=octree -D_MPI=1 evaporation_ns.c -o evaporation_ns -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
 *      
 * 2D with MPI and graphics: CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 evaporation_ns.c -o evaporation_ns -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
 * 
 * 2D without MPI: qcc -Wall -O2 droplet.c -o droplet -lm
 *  
 * Run with:
 * 
 * with MPI: mpirun -np 4 ./evaporation_ns
 * without MPI: ./droplet
 * 
 * Optional:
 * 
 * qcc -source droplet.c
 * mpicc -O2 -Wall -std=c99 -D_MPI=1 _droplet.c -o droplet -lm
 * mpirun -np 2 ./droplet
 * 
 */


#define R0 1
#define L 50. // size of the box

#define MIN_LEVEL 5
#define LEVEL 5
#define MAX_LEVEL 11
#define dR_refine (2.*L0/(1 << LEVEL))

#define F_ERR 1e-10

#define T_END 20 //600.
//#define DT_MAX 10.
#define DELTA_T 1 //10. // for measurements and videos


//#include "../advection_Q.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "two-phase.h"
#include "tension.h"


#include "../../elementary_body.h"
#include "view.h"
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

#define vapor_peclet 1e-3
#define D_V 100.
#define vcs 1.
#define cinf 0.2
#define dirichlet_time_factor 10.

#define MU_L 1.e-3
#define MU_G 1.e-6
#define RHO_L 10.
#define RHO_G 1.
#define SIGMA 0.0728
#define PRESSURE 101000

// Gravity acceleration
#define GRAVITY 0.1


scalar vapor[];
scalar * tracers = {vapor};
face vector uf_save[];


vapor[left] = dirichlet(cinf);
u.n[left] = neumann(0);
p[left] = dirichlet(PRESSURE);

vapor[right] = dirichlet(cinf);
u.n[right] = neumann(0);
p[right] = dirichlet(PRESSURE);

vapor[top]   = dirichlet(cinf);
u.n[top] = neumann(0);
p[top] = dirichlet(PRESSURE);

vapor[bottom]   = dirichlet(cinf);
u.n[bottom] = neumann(0);
p[bottom] = dirichlet(PRESSURE);


int main() {
  
  size (L);
  origin (-L/2, -2*L/3);
  
  rho1 = RHO_L;
  rho2 = RHO_G;    

  mu1 = MU_L;
  mu2 = MU_G;
  
  f.sigma = SIGMA;

  // Decreasing the tolerance on the Poisson solve improves the results 
  TOLERANCE = 1e-4;
  
  N = 1 << MAX_LEVEL; //LEVEL;
  init_grid (N);
  //DT = DT_MAX;
  
  run();
  
}


// Add the acceleration of gravity in the downward (-y) direction
event acceleration (i++) {  
  face vector av = a;
  foreach_face(y)
    av.y[] -= GRAVITY;   
}


#define circle(x, y, R) (sq(R) - sq(x) - sq(y))

event init (i = 0) {
//   #if TREE
//     refine ( (sq(x)+sq(y)<sq(1.1*R0)) && (sq(x)+sq(y)>sq(0.9*R0)) && level<MAX_LEVEL );
//   #endif
  fraction (f, circle(x, y, R0));
  foreach()
    vapor[] = f[]*vcs + (1. - f[])*cinf;
//   foreach_face()
//     uf.x[] = 0.;
//   boundary({vapor, uf});
 // CFL = 0.2;
}


event stability (i++) {
  
  foreach_face()
    uf_save.x[] = uf.x[];
  boundary((scalar*){uf_save});
  
  face vector ev[];
  
  vapor.D = D_V;
  vapor.peclet = vapor_peclet;
  vapor.inverse = true;
  
  phase_change_velocity (f, vapor, ev);
  boundary((scalar*){ev});
  
  foreach_face()
    uf.x[] += ev.x[];
  boundary((scalar*){uf});
  
}


event tracer_advection (i++) {
  foreach_face()
    uf.x[] = uf_save.x[];
  boundary((scalar*){uf});
}


event tracer_diffusion(i++) {

  #if TREE
    int max_level = MAX_LEVEL;
  #else
    int max_level = LEVEL;
  #endif
  
  vapor.D = D_V;
  vapor.tr_eq = vcs;
  vapor.inverse = true;
  dirichlet_diffusion (vapor, f, max_level, dt, dirichlet_time_factor);
}


#if TREE
event adapt (i++) {
  adapt_wavelet ({f, u}, (double[]){1e-3, 1e-2, 1e-2}, minlevel = MIN_LEVEL, maxlevel = MAX_LEVEL);
}
#endif


event interface (t = 3.*T_END/4.) {
  static FILE * fpshape = fopen("shape", "w");
  output_facets (f, fpshape);
  fflush(fpshape);
}


event log_console (i++) {
  if (i == 0)
    fprintf (fout, "t grid->tn dt mgp.i mgu.i\n");
  fprintf (fout, "%g %ld %g %d %d\n", t, grid->tn, dt, mgp.i, mgu.i);
  fflush (fout);
}


event slices (t=0; t += DELTA_T; t <= T_END) {
      
  // Vitor: You may want to build an animation further 
  // $convert -delay 10 -loop 0 *.ppm diffusion2D.gif
  
  char ppm_name_vof[80];
  sprintf(ppm_name_vof, "vof-%2.1f.ppm", t);
  output_ppm(f, file=ppm_name_vof, n=1<<MAX_LEVEL);
  
  char ppm_name_v[80];
  sprintf(ppm_name_v, "v-%2.1f.ppm", t);
  output_ppm(u.y, file=ppm_name_v, n=1<<MAX_LEVEL);
  
  char ppm_name_concentration[80];
  sprintf(ppm_name_concentration, "c-%2.1f.ppm", t);
  output_ppm(vapor, file=ppm_name_concentration, min=cinf, max=vcs, n=1<<MAX_LEVEL);
  
//   char ppm_name_temperature[80];
//   sprintf(ppm_name_temperature, "t-%2.1f.ppm", t);
//   output_ppm(temperature, file=ppm_name_temperature, min=ts, max=tinf, n=1<<MAX_LEVEL);
  
  char ppm_name_level[80]; 
  sprintf(ppm_name_level, "level-%2.1f.ppm", t);
  scalar ll[];
  foreach() {
    ll[] = level;
  }
  boundary({ll});
  output_ppm(ll, file=ppm_name_level, max=MAX_LEVEL, n=1<<MAX_LEVEL);
  
  #if _MPI
  scalar pid[];
  foreach() 
    pid[] = fmod(pid()*(npe() + 37), npe());    
  boundary ({pid});
  char ppm_name_pid[80]; 
  sprintf(ppm_name_pid, "pid-%2.1f.ppm", t);
  output_ppm(pid, file=ppm_name_pid, n=1<<MAX_LEVEL);
  #endif
  
}


event outputs (t = 0.; t += DELTA_T; t <= T_END) { 

  double effective_radius;
  effective_radius = pow(3.*statsf(f).sum, 1./3.);

  fprintf (stderr, "%.17g %.17g\n", t, effective_radius);
  fflush(stderr);
  //if (t <= T_END - 20.) {
    static FILE * fpvapor = fopen("vapor_profile", "w");
    static FILE * fpvaporresc = fopen("vapor_profile_resc", "w");
    for (double y = 0.; y <= 5.; y += 0.01)
      fprintf (fpvapor, "%g %g\n", y, interpolate (vapor, 0., y));
    for (double y = effective_radius; y <= 5.; y += 0.01)
      fprintf (fpvaporresc, "%g %g %g\n", effective_radius, y/effective_radius,
               interpolate (vapor, 0., y));
    fprintf (fpvapor, "\n");
    fprintf (fpvaporresc, "\n");
    fflush(fpvapor);
    fflush(fpvaporresc);
  //}
  
  scalar vapor_draw[];
  foreach() {
    f[] = clamp(f[], 0., 1.);
    vapor_draw[] = - vapor[];
  }
  boundary({f, vapor_draw});

  view (fov = 15, width = 640, height = 640, samples = 1, relative = false,
        tx = 0., ty = 0., bg = {BG, BG, BG});
  clear();
  draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
           fc = {BG, BG, BG});
  squares ("vapor_draw", min = - vcs, max = vcs, linear = false,
           map = cool_warm);
  mirror (n = {1., 0., 0.}, alpha = 0.) {
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
             fc = {BG, BG, BG});
    squares ("vapor_draw", min = - vcs, max = vcs, linear = false,
             map = cool_warm);
  }
  mirror (n = {0., 1., 0.}, alpha = 0.) { // vapor
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
             fc = {BG, BG, BG});
    squares ("vapor_draw", min = - vcs, max = vcs, linear = false,
             map = cool_warm);
    mirror (n = {1., 0., 0.}, alpha = 0.) {
        draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
                 fc = {BG, BG, BG});
        squares ("vapor_draw", min = - vcs, max = vcs, linear = false,
                 map = cool_warm);
    }
  }
  save ("video_static_drop.mp4");
}





