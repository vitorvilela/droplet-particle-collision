/*
 * This Basilisk setup is intended to simulate the droplet-particle impact cases described in:
 * 
 * Mitra S and Evans G (2018) Dynamic Surface Wetting and Heat Transfer in a Droplet-Particle System of Less Than Unity Size Ratio. 
 * Front. Chem. 6:259. doi: 10.3389/fchem.2018.00259
 *  
 */

/*
 * Summary:
 * 
 * Droplet-Particle collision with and without heat transfer (and evaporation) 
 * Transport formulation: conservative approach
 * Immersed boundary: imposed momentum
 * Contact-angle: no model
 * 
 * Warning: check if you need to create a ./results folder
 *
 * Compile with:
 * 
 * 3D with MPI and graphics: CC99='mpicc -std=c99' qcc -Wall -O2 -grid=octree -D_MPI=1 droplet.c -o droplet -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
 *      
 * 2D with MPI and graphics: CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 droplet.c -o droplet -lm -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11
 * 
 * 2D without MPI: qcc -Wall -O2 droplet.c -o droplet -lm
 *  
 * Run with:
 * 
 * with MPI: mpirun -np 2 ./droplet
 * without MPI: ./droplet
 * 
 * Optional:
 * 
 * qcc -source droplet.c
 * mpicc -O2 -Wall -std=c99 -D_MPI=1 _droplet.c -o droplet -lm
 * mpirun -np 2 ./droplet
 * 
 */

#include "momentum.h"
#include "tension.h"
//#include "view.h"

#if dimension == 3
#include "lambda2.h"
#endif

#define cf



// Domain

// The computational domain (12mm × 12mm × 9mm) for 3D simulations

#define LARGER 2 
#define DOMAIN_LENGTH LARGER*12.e-3

// Gravity acceleration
#define GRAVITY 9.81

// Ambient pressure at boundary
#define PRESSURE 101000


// Operating conditions used in the cold state interactions

// d{d,0} (mm)
// 2.9 ± 0.1

// Case     v{0} (m.s−1)     Re     We     Ca
// 1        0.32             932    4      0.004
// 2        0.54             1594   12     0.007
// 3        0.70             2039   19     0.009

#define CASE 1  

// Velocity of the droplet just before the impact
#define V0 0.32

// Reference length is the droplet diameter
#define D0 2.9e-3



// Fluid properties 
// Liquid: water, Gas: air, Solid: glass

// Cold state ambient condition (20 °C)
// Film boiling regime (350 °C)

// Cold state, physical properties of droplet used to report the dimensionless numbers are ρl = 998.2 kg. m−3; µ = 0.001 kg.m−1.s−1, σlg = 0.073N.m−1

// All the thermo-physical properties (density, viscosity, surface tension, heat capacity, and thermal conductivity) of the gas and liquid phase used in the simulations were set as temperature dependent polynomials

// At 20 °C and atmospheric pressure, water has @engineersedge.com:
// density of 998.21 kg/m³
// dynamic viscosity of 1.002e-3 kg/m.s 
// specific heat at constant pressure of 4181.8 J/kg.K
// thermal conductivity of 0.5984 W/m.k

// At 20 °C and atmospheric pressure, dry air has @theengineeringmindset.com:
// density of 1.2047 kg/m³
// dynamic viscosity of 1.8205e-5 kg/m.s 
// specific heat at constant pressure of 1006.1 J/kg.K
// thermal conductivity of 0.025596 W/m.k

#define RHO_L 998.2
#define RHO_G 1.2
#define RHO_S 500.

#define MU_L 1.0e-3
#define MU_G 1.8e-6

// Water–Air
// Temperature-Sigma (mN·m−1)
// 20 °C 72.86 ± 0.05
// 21.5 °C 72.75
// 25 °C 71.99
#define SIGMA 0.073



// Particle

// On a 10 mm solid brass sphere
// Before each droplet deposition, the sphere surface was carefully cleaned with acetone and allowed for sufficient time to dry
// Surface temperature was controlled by a PID controller connected with an embedded T-type thermocouple and a cartridge heater placed in a well- insulated billet

scalar particle[];

#define Dp 10.e-3
#define Xp 0.
#define Yp 0.
#define Up 0.
#define Vp 0.
#if dimension == 3
#define Zp 0.
#define Wp 0.
#endif



// Droplet

// Filtered water droplets of diameter ~2.9 +-0.1 mm
// A droplet delivery system with adjustable height, using a 21G hypodermic nozzle and a precision syringe pump, (∼10–150mm from the apex point of the sphere) was utilized

#define Dd D0
#define Xd 0.
const double Yd = Yp + 0.515*Dp + 0.5*Dd; // considering an additional gap of 1.5%Dp 
#define Ud 0.
#define Vd V0
#if dimension == 3
#define Zd 0.
#define Wd 0.
#endif

#if dimension == 2
const double initial_area = pi*Dd; // circunference (m)
#else
const double initial_area = pi*Dd*Dd; // superficial area (m²)
#endif



// Dimensionless parameters

#define Re RHO_L*V0*Dd/MU_L
#define We RHO_L*V0*V0*Dd/SIGMA 



// Spatial resolution

// The mesh comprising ~0.32 million hexahedral cells
// Size of cells were kept lower in the vicinity of spherical surface for better resolution of the three-phase contact line and gradually coarser away from the solid surface
// Total 10,592 cells were patched to resolve the droplet

#define SPATIALSCALER 1.e6

// Refinement level (level, n): (5, 32) (6, 64) (7, 128) (8, 256) (9, 512) (10, 1024) (11, 2048) (12, 4096)
int LEVEL = 5;
int LEVEL_MAX = 8;


// Temporal resolution

// Droplet-sphere interactions were captured using a Phantom v311 high speed camera at 2000 frames per second in shadowgraphymode using backlighting and a diffuser screen
// Total cell number was decided based on a trade- off between a reasonable agreement with the experimental data and computational time which took on average ∼2–4 days per workstation.

// All simulations were performed for a duration of 10–20ms depending on impact Weber number using a time step of 10−6 s

#define LAST_FRAME 20.e-3
#define TIMESCALER 1.e3
#define SMALL 1.e-7

// Checking times and simulation duration
const double checkat = 1./2000;
const double print3Dat = 1.e-3; 
const double print2Dat = 1.e-3;
const double duration = 1.01*LAST_FRAME;



// Functions for analyses

// Droplet boundary was marked to separate it from the background and area equivalent diameter was determined
// Centroid of the marked droplet was tracked prior to impact to estimate the impingement velocity
// Contact angles were determined on both left and right side of the marked interface by computing the inside angle between the two tangents—one drawn to the interface and the other on the spherical surface both passing through the three-phase contact line intersection points



int main() {
    
  L0 = DOMAIN_LENGTH;

  #if dimension == 2    
  origin(-0.5*DOMAIN_LENGTH, 0);
  #else
  origin(-0.5*DOMAIN_LENGTH, 0, -0.5*DOMAIN_LENGTH);
  #endif    
    
  rho1 = RHO_L;
  rho2 = RHO_G;    

  mu1 = MU_L;
  mu2 = MU_G;
  
  f.sigma = SIGMA;

  // Decreasing the tolerance on the Poisson solve improves the results 
  TOLERANCE = 1e-4;
        
  // Initial uniform mesh size 
  N = 1 << LEVEL;  
  init_grid(N);
  int M = 1 << LEVEL_MAX;
  
   
  printf ("\nCase: %d\tRe: %d\tWe: %d\nCoarser and finer grids (micrometers): %f\t%f\n", CASE, (int)(Re), (int)(We), (DOMAIN_LENGTH/N)*SPATIALSCALER, (DOMAIN_LENGTH/M)*SPATIALSCALER);
  
  //#if _OPENMP || _MPI
  //if (pid() == 0) {
  //#endif  
    //char msg[80];
    //sprintf (msg, "\nCase: %d\tRe: %d\tWe: %d\nCoarser and finer grids (micrometers): %f\t%f\n", CASE, (int)(Re), (int)(We), (DOMAIN_LENGTH/N)*SPATIALSCALER, (DOMAIN_LENGTH/M)*SPATIALSCALER);
    //printf ("%s", msg);
    //static FILE * fmsg = fopen ("stdout.txt", "w");
    //fprintf (fmsg, "%s", msg);  
    //fflush (fmsg);
  //#if _OPENMP || _MPI  
  //}
  //#endif   
  
  run(); 
  
}



// Boundary conditions

// x
q.n[right] = neumann(0);
p[right] = dirichlet(PRESSURE);

q.n[left] = neumann(0);
p[left] = dirichlet(PRESSURE);

// y
q.n[bottom] = dirichlet(0);
p[bottom] = neumann(0);

q.n[top] = neumann(0);
p[top] = dirichlet(PRESSURE);

// z
#if dimension == 3
q.n[front] = neumann(0);
p[front] = dirichlet(PRESSURE);

q.n[back] = neumann(0);
p[back] = dirichlet(PRESSURE);
#endif



event init(i = 0) {
  
  if (!restore (file = "dump")) {  
    
    #if dimension == 2
    coord vd = {Ud, -Vd};   
    refine ( ( sq(x-Xp) + sq(y-Yp)< sq(1.1*0.5*Dp) ||  sq(x-Xd) + sq(y-Yd) < sq(1.1*0.5*Dd) ) && level < LEVEL_MAX );   
    fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );  
    fraction ( f, sq(0.5*Dd) - sq(x-Xd) - sq(y-Yd) );
    #else
    coord vd = {Ud, -Vd, Wd};
    // Particle and Droplet Refinement
    refine ( ( sq(x-Xp) + sq(y-Yp) + sq(z-Zp) < sq(1.1*0.5*Dp) ||  sq(x-Xd) + sq(y-Yd) + sq(z-Zd) < sq(1.1*0.5*Dd) ) && level < LEVEL_MAX );
    // Particle Fraction
    fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) - sq(z-Zp) );
    // Droplet Fraction
    // in fractions.h
    // #define fraction(f,func)
    // func = positive inside and negative outside circle/sphere (droplet)
    // void fractions { c[] = ni ? plane_volume (n, α/ni) : s.x[]; }
    fraction ( f, sq(0.5*Dd) - sq(x-Xd) - sq(y-Yd) - sq(z-Zd) );
    #endif
    
    foreach() 
      foreach_dimension()      
        q.x[] = f[]*rho1*vd.x;
    boundary ((scalar *){q}); 
    
#ifndef cf
    foreach()
      cf[] = f[];
#endif
    
  }
  
}



event stationary_particle (i++) {
  
  #if dimension == 2
  coord vp = {Up, Vp};  
  fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );
  #else
  coord vp = {Up, Vp, Wp};  
  // Stationary particle
  fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) - sq(z-Zp) );
  #endif
  
  // We then use this (solid) volume fraction field to impose the corresponding velocity (momentum) field inside the solid (as a volume-weighted average)  
  foreach()
    foreach_dimension()
      q.x[] = particle[]*RHO_S*vp.x + (1. - particle[])*q.x[];
  boundary ((scalar *){q});
  
}



// Add the acceleration of gravity in the downward (-y) direction
event acceleration (i++) {
  
  face vector av = a;
  foreach_face(y)
    av.y[] -= GRAVITY;  
  
}




// Mesh adaptivity

event adapt (i++) {
  
  // Extract velocity from momentum
  vector u[];
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  boundary((scalar *){u});
  
#if dimension == 2
  adapt_wavelet({f, particle, u}, (double[]){1.e-3, 1.e-3, 1.e-2, 1.e-2}, LEVEL_MAX);
#else
  // Refinement threshold vector for f, particle and u vector
  adapt_wavelet({f, particle, u}, (double[]){1.e-3, 1.e-3, 1.e-2, 1.e-2, 1.e-2}, LEVEL_MAX);  
#endif  
  
  //unrefine (y<(1/2)*DOMAIN_LENGTH);
  
}



event logfile (i++) {
  
  static FILE * flf = fopen ("simulation.csv", "w");
    
  if (i == 0)
    fprintf (flf, "t, t*, dt, mgp.i, mgu.i, grid->tn, perf.t, perf.speed\n");
  
  fprintf (flf, "%g, %g, %g, %d, %d, %ld, %g, %g\n", 
                t, t*V0/D0, dt, mgp.i, mgu.i, grid->tn, perf.t, perf.speed);
  
  fflush (flf);
  
}




// Check droplet position at regular intervals and save .ppm images of the simulation

event snapshot (t=0; t+=print2Dat; t<duration) {
    
  char ppm_name_vof[80];
  sprintf(ppm_name_vof, "vof-%2.3f.ppm", t*TIMESCALER); // naming as miliseconds
  
  char ppm_name_u[80];
  sprintf(ppm_name_u, "u-%2.3f.ppm", t*TIMESCALER);   
  
  char ppm_name_lambda2[80];
  sprintf(ppm_name_lambda2, "lambda2-%2.3f.ppm", t*TIMESCALER); 
  
  char ppm_name_level[80]; 
  sprintf(ppm_name_level, "level-%2.3f.ppm", t*TIMESCALER);
    
  // VOF cut-off
  scalar ff[];
  foreach() {
    ff[] = f[] < 1e-4 ? 0 : f[] > 1. - 1e-4 ? 1. : f[];
  }
  boundary({ff});
    
  scalar ll[];
  foreach() {
    ll[] = level;
  }
  boundary({ll});
    
  // Cells for which m is negative will be black  
  scalar m[];
  foreach() {
    m[] = 0.5 - particle[];
  }
  boundary ({m});
  
  vector u[];  
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  boundary((scalar *){u});
  
  scalar l2[];
  lambda2 (u, l2);
      
  output_ppm(ff, file=ppm_name_vof, mask=m, n=1<<LEVEL_MAX);
  output_ppm(u.x, file=ppm_name_u, n=1<<LEVEL_MAX);
  output_ppm(l2, file=ppm_name_lambda2, n=1<<LEVEL_MAX);
  output_ppm(ll, file=ppm_name_level, max=LEVEL_MAX, n=1<<LEVEL_MAX);
    
}


event snapshot3D (t=0; t+=print3Dat; t<duration) {
  
  // We dump a snapshot which can be used to restart the simulation
  char name[80];
  sprintf (name, "dump-%2.3f", t*TIMESCALER); // naming as miliseconds
  dump (file = name);  
  
}
