// http://basilisk.fr/sandbox/antko/evaporation/evaporation.c

/* Compile with:
 * 
 * 3D with MPI and graphics: CC99='mpicc -std=c99' qcc -Wall -O2 -grid=octree -D_MPI=1 evaporation.c -o evaporation -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
 *      
 * 2D with MPI and graphics: CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 droplet.c -o droplet -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
 * 
 * 2D without MPI: qcc -Wall -O2 droplet.c -o droplet -lm
 *  
 * Run with:
 * 
 * with MPI: mpirun -np 4 ./evaporation
 * without MPI: ./droplet
 * 
 * Optional:
 * 
 * qcc -source droplet.c
 * mpicc -O2 -Wall -std=c99 -D_MPI=1 _droplet.c -o droplet -lm
 * mpirun -np 2 ./droplet
 * 
 */

// Drop evaporation

// We investigate the evaporation of a single droplet in a still and relatively dry environment

// For this investigation we use – in an axisymmetric setting – the VOF description of the interface, the Bell-Collela-Glaz advection solver and the fully implicit reaction-diffusion solver

#if 0
#include "axi.h"
#endif

#include "advection.h"
#include "vof.h"
#include "diffusion.h"

// We allocate several scalars and vector fields to describe both the interface and the concentration field

scalar f[];
scalar * tracers = NULL, * interfaces = {f};
scalar concentration[], temperature[];

// We non-dimensionalise the problem with:
// the initial radius of the drop R, 
// the diffusion coefficient of liquid vapour in air D
// and the saturation concentration of water in air c0
// (this quantity having the dimension of a density M.L^-3)

// The diffusion equation for the vapour in the non-dimensional space reads:
// dc/dt = Delta c,
// appropriately completed with Dirichlet conditions at the surface of the drop and at infinity
// c=1 : at the drop surface
// c=coo/c0 : far from the drop

// As the interface recedes at a (dimensional) velocity
// ve = (D/rho).Delta c ~ (D/rho).(c0/R) = V
// a Peclet number Pe = VR/D = c0/rho enters in the problem description
// Typically, for the problem of a water droplet evaporating in dry air, this Peclet number is O(10^-5), meaning that the problem is really dominated by diffusion

// Here, we choose Pe=10^-3, and set the relative humidity of the surroundings (cinf) to 20%

const double Peclet_number = 1e-3;
const double cs = 1.0, ts = 0.;
const double cinf = 0.2, tinf = 1.;

// We set some numerical constants and placeholders for e.g. the effective radius of the droplet

const double tend = 1000.;

// Refinement level (level, n): (5, 32) (6, 64) (7, 128) (8, 256) (9, 512) (10, 1024)
int LEVEL = 5, MINLEVEL = 5, MAXLEVEL = 9;

double effective_radius = 1.;
double dt, dirichlet_time_constant;

mgstats mgconcentration, mgtemperature;

// Thanks to symmetry, we only solve a quarter of the domain, and requires the concentration to drop at its asymptotic value “at infinity” (that is, at the box boundary)

concentration[left] = dirichlet(cinf);
concentration[right] = dirichlet(cinf);
concentration[top] = dirichlet(cinf);
concentration[bottom] = dirichlet(cinf);
concentration[front] = dirichlet(cinf);
concentration[back] = dirichlet(cinf);

temperature[left] = dirichlet(tinf);
temperature[right] = dirichlet(tinf);
temperature[top] = dirichlet(tinf);
temperature[bottom] = dirichlet(tinf);
temperature[front] = dirichlet(tinf);
temperature[back] = dirichlet(tinf);

// The main function of the program (where we set the domain geometry to be ten times larger than the drop. Vitor: "drop radius")

int main() {
  
  size (20.);
  
  origin (-10., -10., -10.); // (0., 0.): left-bottom corner
  
  N = 1 << LEVEL;
  
  init_grid (N);
  
  run();
  
}

// The initial position of the interface is defined with this function

// r^2 = x^2 + y^2 + z^2      
#define circle(x,y,z) (sq(1.) - sq(x) - sq(y) -sq(z))

// Before the first step, we initialise the concentration field to c0 within the drop, and coo everywhere else. 
// The diffusion coefficient is set to one, but beware that in Basilisk the diffusion coefficient also incorporates information about the metrics, namely the face length fm.x()

event init (t = 0) {
  
  refine ( (sq(x)+sq(y)+sq(z)<sq(1.1*0.5*2)) && (sq(x)+sq(y)+sq(z)>sq(0.9*0.5*2)) && level < MAXLEVEL );
  
  fraction (f, circle(x,y,z));
  
  foreach() {
    concentration[] = f[]*cs + (1 - f[])*cinf;
    temperature[] = f[]*ts + (1 - f[])*tinf;
  }
  boundary({concentration, temperature});
  
  CFL = 0.2;
  
}

// A function to rescale normals so that they are unit vectors w.r.t. the 2-norm (by default, the 1-norm is adopted for efficiency purposes)

coord normal (Point point, scalar c) {
  
  coord n = mycs (point, c);
  
  double nn = 0.;
  
  foreach_dimension()
    nn += sq(n.x);
  nn = sqrt(nn);
  
  foreach_dimension()
    n.x /= nn;
  
  return n;
  
}


// Notes on metrics from http://basilisk.fr/sandbox/amansur/Rising/2D_axi
// In 2D axisymmetric, the cell metric cm changes, the volume/area of a cell is proportional to r (i.e. y). 
// The simulation is based on 2D parameters, to get the 2D axisymmetry each cell has to be multiplied by 2.pi.r. 
// cm[] becomes y, but not 2.pi.r!

// cm[] = y

// Also, the length scale factor or the face metric fm changes. This face metric is used when dealing with flux. Each flux should be 2.pi.r.Delta either on the x or y axis. 
// Here again the 2.pi.r is missing.

// fm.x[] = fm.y[] = y


// The core routine of the code. Several key operations are performed here:

// - Computation of the concentration field,
// - Derivation of the concentration gradient,
// - Setting of the advection velocity field for the interface according to Fick’s law.

// An important point is how the immersed Dirichlet boundary condition is handled here. To ensure that the concentration within the drop stays at c0
// and is not washed away, we introduce a source term in the diffusion equation. 
// This source term is only active within the drop and is as large as the difference between the actual concentration and c0 is. This modified equation therefore reads

// dc/dt = Delta c + (f/tal).(c0-c)

// Here 'f' is the liquid phase fraction, and 'tal' is a time constant, which has to be shorter than the typical diffusive timescale in order for the control to be effective. 
// In the following, we set 'tal' to a tenth of the smallest diffusive time computed at the cell level

event velocity (i++) {
  
  // The processes at play are slow for most of the simulation, and the solver is fully implicit so we can go full throttle and set huge timesteps like 10. 
  // However in the final step the velocity is expected to diverge, so we have to be careful and bound the timestep
  
  dt = dtnext (timestep (u, 10.));
  dirichlet_time_constant = 1e-1*sq(L0/(1 << MAXLEVEL));
  
  scalar volumic_metric[], dirichlet_source_term[], dirichlet_feedback_term[];
  
  face vector diffusion_coef[], concentration_flux[];
  face vector thermalDiff_coef[];
  
  foreach() {
    volumic_metric[] = cm[];
    dirichlet_source_term[] = cm[]*cs*f[]/dirichlet_time_constant;
    dirichlet_feedback_term[] = -cm[]*f[]/dirichlet_time_constant;
  }
  
  foreach_face() {
    diffusion_coef.x[] = fm.x[];
    // Vitor: Arbitrary thermal diffusion coefficient 0.001x concentration diffusion coefficient
    thermalDiff_coef.x[] = 0.001*fm.x[];
  }
  
  boundary({volumic_metric,dirichlet_source_term,dirichlet_feedback_term,diffusion_coef, thermalDiff_coef});
  
  mgconcentration = diffusion (concentration, dt, D=diffusion_coef, r=dirichlet_source_term, beta=dirichlet_feedback_term, theta=volumic_metric);
  mgtemperature = diffusion (temperature, dt, D=thermalDiff_coef);
  
  
  // Having the concentration we derive the concentration gradient, which is a face vector field:
  
  foreach_face()
    concentration_flux.x[] = (concentration[] - concentration[-1])/Delta;  
  boundary((scalar*){concentration_flux});
  
  // With the concentration gradient we can now advance the interface, following the lines drawn in meanflow.c. 
  // We define in particular the velocity of the interface as the product between the Peclet number, the local normal flux and the normal. 
  // Note that u.x is weighted by the metric (embodied in the diffusion coefficient).
  
  foreach_face() {
    
    u.x[] = 0.;
    
    if (f[] > 0. && f[] < 1.) {
      
      coord n = normal (point, f);
      
      u.x[] = (n.x > 0. ?
	       Peclet_number*diffusion_coef.x[]*concentration_flux.x[1] :
	       Peclet_number*diffusion_coef.x[]*concentration_flux.x[-1]);
      
    }
    else if (f[-1] > 0. && f[-1] < 1.) {
      
      coord np = normal (neighborp(-1), f);
      
      u.x[] = (np.x > 0. ?
	       Peclet_number*diffusion_coef.x[]*concentration_flux.x[1] :
	       Peclet_number*diffusion_coef.x[]*concentration_flux.x[-1]);
      
    }
  }
  boundary((scalar*){u});
  
}
  
// We now write several post-processing events, e.g. to track the effective radius of the drop:
  
event compute_radius (i++) {
  stats s = statsf (f);
  effective_radius = pow(3.*s.sum, 1./3.);
}

event log (i++) {
  if (i == 0)
    fprintf (stdout, "t grid->tn\n");
  fprintf (stdout, "%g %ld\n", t, grid->tn);
}

event logfraction (t += 1.) {
  
  static FILE * fr = fopen ("radius.csv", "w");
  if (i == 0)
    fprintf (fr, "t, r\n");  
  fprintf (fr, "%f %.12f\n", t, effective_radius);
  fflush (fr);
  
}

// During the lifetime of the drop, we also track the concentration and fraction profiles.

event logfile (t = 30; t += 10.; t < 580) {
  
  static FILE * fpinterf = fopen("interfaces.dat", "w");
  output_facets (f, fpinterf);
  fflush (fpinterf);

  static FILE * fp = fopen("cprof.dat", "w");
  static FILE * fpcresc = fopen("cprof_resc.dat", "w");
  static FILE * fpf = fopen("fprof.dat", "w");
  for (double y = 0.; y <= 5.; y += 0.01) {
    fprintf (fp, "%g %g\n", y,
      interpolate (concentration, 0., y));
    fprintf (fpf, "%g %g\n", y,
      interpolate (f, 0., y));
  }
  for (double y = effective_radius; y <= 5.; y += 0.01) {
    fprintf (fpcresc, "%g %g %g\n", effective_radius, y/effective_radius,
      interpolate (concentration, 0., y));
  }
  fprintf (fp, "\n");
  fprintf (fpcresc, "\n");
  fprintf (fpf, "\n");
  fflush (fp);
  fflush (fpcresc);
  fflush (fpf);
  
}


event output (t=0; t+=100; t<tend) {
      
  // Vitor: You may want to build an animation further 
  // $convert -delay 10 -loop 0 *.ppm diffusion2D.gif
  
  char ppm_name_vof[80];
  sprintf(ppm_name_vof, "vof-%2.1f.ppm", t);
  output_ppm(f, file=ppm_name_vof, n=1<<MAXLEVEL);
  
  char ppm_name_concentration[80];
  sprintf(ppm_name_concentration, "c-%2.1f.ppm", t);
  output_ppm(concentration, file=ppm_name_concentration, min=cinf, max=cs, n=1<<MAXLEVEL);
  
  char ppm_name_temperature[80];
  sprintf(ppm_name_temperature, "t-%2.1f.ppm", t);
  output_ppm(temperature, file=ppm_name_temperature, min=ts, max=tinf, n=1<<MAXLEVEL);
  
  char ppm_name_level[80]; 
  sprintf(ppm_name_level, "level-%2.1f.ppm", t);
  scalar ll[];
  foreach() {
    ll[] = level;
  }
  boundary({ll});
  output_ppm(ll, file=ppm_name_level, max=MAXLEVEL, n=1<<MAXLEVEL);
  
  #if _MPI
  scalar pid[];
  foreach() 
    pid[] = fmod(pid()*(npe() + 37), npe());    
  boundary ({pid});
  char ppm_name_pid[80]; 
  sprintf(ppm_name_pid, "pid-%2.1f.ppm", t);
  output_ppm(pid, file=ppm_name_pid, n=1<<MAXLEVEL);
  #endif
  
}




#if 1
event adapt (i++) {
  adapt_wavelet ({f, concentration, temperature}, (double[]){1e-3, 1e-3, 1e-3}, maxlevel = MAXLEVEL);
}
#endif
  
  
  
event snapshot3D (t=0; t+=100; t<tend) {
  
  // We dump a snapshot which can be used to restart the simulation and plot with OpenGL
  char name[80];
  sprintf (name, "dump-%2.1f", t);
  dump (file = name);  
  
}  
// Results...  
  
  

#if 0
event gfsview (i += 10) {
  foreach()
    f[] = clamp(f[],0,1);
  boundary ({f});
  static FILE * fpview = popen ("gfsview2D evaporation.gfv", "w");
  output_gfs (fpview);
}
#endif

#if 0
event movies (t += 15; t <= tend) {
  static FILE * fpmov =
    popen ("gfsview-batch2D evaporation.gfv | ppm2gif > evap.gif", "w");

  // we clamp volume fraction because gfsview does not like over/undershoots
  foreach()
    f[] = clamp(f[],0,1);
  boundary ({f});
  output_gfs (fpmov);
  
  fprintf (fpmov, "Save stdout { format = PPM width = 300 height = 300}\n");
}
#endif


