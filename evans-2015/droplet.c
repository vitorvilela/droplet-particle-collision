/*
 * Droplet-Particle collision - Three-phase VoF / Non-conservative approach
 *
 * Compile with:
 * 
 *      3D with MPI: CC99='mpicc -std=c99' qcc -Wall -O2 -grid=octree -D_MPI=1 droplet.c -o droplet -lm
 *      2D with MPI: CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 droplet.c -o droplet -lm
 *      2D without MPI: qcc -Wall -O2 droplet.c -o droplet -lm
 * 
 * Run with:
 * 
 *      with MPI: mpirun -np 2 ./droplet >> ./results/out.dat
 *      without MPI: ./droplet
 * 
 * Optional:
 * 
 * qcc -source droplet.c
 * mpicc -O2 -Wall -std=c99 -D_MPI=1 _droplet.c -o droplet -lm
 * mpirun -np 2 ./droplet
 * 
 */

/*
 * This Basilisk setup is intended to simulate the droplet-particle collision cases described in:
 * 
 * Journal. 
 * Title.
 * Author, year.
 * 
 * Objective: 
 * Methods/Models:
 * Analyses: 
 * 
 */

#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "tension.h"

/*
 * Length scale
 */
#define D0 5.98e-3

/*
 * Fluid flow properties
 */
#define RHO_L 998.2
#define RHO_G 1.
#define RHO_S 2450.

#define MU_L 1.e-3
// Using a greater MU_G until droplet rests on the still particle
#define MU_G 1.e-6
//#define MU_G 1.e-6
#define MU_S 1.e-1 

// SIGMA_12=SIGMA_1+SIGMA_2 for Liquid-Gas
#define SIGMA_1 0.0728        
#define SIGMA_2 0.00001      
#define SIGMA_3 0.9          

#define PRESSURE 101000

#define LARGER 3
#define HIGHER 1.5

/* 
 * Domain lenght
 */
#define DOMAIN_LENGTH LARGER*10*D0

/*
 * Polypropylene still particle
 */
#define Dp D0
#define Xp 0.
#define Yp HIGHER*0.5*DOMAIN_LENGTH

/*
 * Droplet
 */
#define Dd 3.41e-3
#define Xd 0.
const double Yd = Yp + 0.5*Dp + 0.51*Dd;

/*
 * Glass moving particle
 */
#define Dg 1.13e-3
#define Xg 0.
const double Yg = (Yp + 0.5*Dp + 0.5*Dd) + 0.5*Dd + 2.5*Dg;




/*
 * Checking time and simulation duration
 */
#define nondimensional_checkat 50000
#define nondimensional_duration 1500000
const double checkat = nondimensional_checkat*Dd/(9.81/Dd); 
const double duration = nondimensional_duration*Dd/(9.81/Dd);
const double stillDroplet = 1000000*Dd/(9.81/Dd);

/*
 * Refinement level
 * (level, n): (4, 16) (5, 32) (6, 64) (7, 128) (8, 256) 
 */
int LEVEL = 5;
int LEVEL_MAX = 11;

double volref[3];

int main() {
    
    /*
     * Domain Length
     */
    L0 = DOMAIN_LENGTH;
        
    
    /* 
     * Left-bottom origin at
     */
    origin(-0.5*DOMAIN_LENGTH, 0);
      
    
    /* 
     * Density of phases 1-Liquid 2-Gas  3-Solid
     */
    rho1 = RHO_L;
    rho2 = RHO_G;    
    rho3 = RHO_S;
    
    /*
     * Dynamic viscosity of phases 1-Liquid 2-Gas 3-Solid
     */
    mu1 = MU_L;
    mu2 = MU_G; //mu2 = MU_G ? t < stillDroplet : MU_G*MU_G;
    mu3 = MU_S;
    
    
    /*
     * Surface tension 
     */
    f1.sigma = SIGMA_1;
    f2.sigma = SIGMA_2;
    f3.sigma = SIGMA_3;
    
    
    /*
     * Decreasing the tolerance on the Poisson solve improves the results
     */    
    TOLERANCE = 1e-4;
    
    
    /* 
     * Initial uniform mesh size 
     */
    N = 1 << LEVEL;
    init_grid(N);

    
    run(); 
  
}


/*
 * Boundary conditions
 */
u.n[bottom] = neumann(0);
p[bottom] = dirichlet(PRESSURE);

u.n[top] = neumann(0);
p[top] = dirichlet(PRESSURE);

u.n[right] = neumann(0);
p[right] = dirichlet(PRESSURE);

u.n[left] = neumann(0);
p[left] = dirichlet(PRESSURE);


scalar particle[];

event init(i = 0) {
    
  
  if (!restore (file = "dump")) {  
    
    refine ( ( fabs(x-Xp) < fabs(1.1*0.5*Dp) ) && level < LEVEL_MAX);
    
    // Polypropylene still particle
//     refine ( ( sq(x-Xp) + sq(y-Yp) < sq(1.1*0.5*Dp) ) && level < LEVEL_MAX);
    fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );
    
    // Droplet
//     refine ( ( sq(x-Xd) + sq(y-Yd) < sq(1.1*0.5*Dd) ) && level < LEVEL_MAX);
    fraction ( f1, sq(0.5*Dd) - sq(x-Xd) - sq(y-Yd) );
        
    // Glass moving particle
//     refine ( ( sq(x-Xg) + sq(y-Yg) < sq(1.1*0.5*Dg) ) && level < LEVEL_MAX);
    fraction ( f3, sq(0.5*Dg) - sq(x-Xg) - sq(y-Yg) );
        
    foreach() {     
        f2[] = 1.0 - f1[] - f3[];
    }
    boundary ({f2}); 
    
    volref[0] = statsf(f1).sum;
    volref[1] = statsf(f2).sum;
    volref[2] = statsf(f3).sum;
    
  }
}


event stationary_particle (i++) {
  
  // The velocity of the particle
  coord vp = {0., 0.}; 
  
  fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );
  
  foreach()
    foreach_dimension()
      u.x[] = particle[]*vp.x + (1. - particle[])*u.x[];
  boundary ((scalar *){u});
  
  if (t < stillDroplet) {
    fraction ( f3, sq(0.5*Dg) - sq(x-Xg) - sq(y-Yg) );
    
    foreach()
      foreach_dimension()
	u.x[] = f3[]*vp.x + (1. - f3[])*u.x[];
    boundary ((scalar *){u});
  }
  
}


/*
 * Add the acceleration of gravity in the downward (-y) direction
 */
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 9.81;  
}


/*
 * Mesh adaptivity
 */
event adapt (i++) {
  
  /*
   * Refinement threshold vector for f, u.x and u.y 
   */
  adapt_wavelet({f1, f3, particle, u}, (double[]){1.e-3, 1.e-3, 1.e-3, 1.e-2, 1.e-2}, LEVEL_MAX);
  
  unrefine (y<(1/2)*DOMAIN_LENGTH); 
  
}


// Every ten timesteps, we output the time, timestep, volume or area, position and velocity of the droplet
event logfile (i += 10) {
    
  // Summation of: centroid times area (2D) or volume (3D) of each liquid (f[]) cell 
  double xd = 0., yd = 0.;
  double xp = 0., yp = 0.;
  
  // Total droplet area or volume (before breakup) or total liquid parcels area or volume
  double sd = 0.;
  double sp = 0.;
  
  // Summation of: velocity times area (2D) or volume (3D) of each liquid (f[]) cell 
  double vdx = 0., vdy = 0.;
  double vpx = 0., vpy = 0.;
      
  foreach(reduction(+:xd) reduction(+:yd) 
	  reduction(+:vdx) reduction(+:vdy) 
	  reduction(+:sd)
	  reduction(+:xp) reduction(+:yp) 
	  reduction(+:vpx) reduction(+:vpy) 
	  reduction(+:sp)) {
    
      // Area (2D) or volume (3D) of each cell
      // In common.h: 
      //   #define dv() (sq(Delta)*cm[]) (2D)
      //   #define dv() (cube(Delta)*cm[]) (3D)
      double dvd = f1[]*dv();
      double dvp = f3[]*dv();
  
      xd += x*dvd;
      yd += y*dvd;    
      vdx += u.x[]*dvd;
      vdy += u.y[]*dvd;    
      sd += dvd;
      
      xp += x*dvp;
      yp += y*dvp;    
      vpx += u.x[]*dvp;
      vpy += u.y[]*dvp;    
      sp += dvp;
    
  }

  // xd/sd - body centroid position
  // vdx/sd - body velocity
  fprintf (fout,
	   "t: %1.1e  dt: %1.1e, darea: %1.1e, xd: %1.1e, yd: %1.1e, vdx: %1.1e, vdy: %1.1e, parea: %1.1e, xp: %1.1e, yp: %1.1e, vpx: %1.1e, vpy: %1.1e\n", 
	   t, dt, sd, xd/sd, yd/sd, vdx/sd, vdy/sd, sp, xp/sp, yp/sp, vpx/sp, vpy/sp);
  fflush (fout);
  
}


/*
 * Check droplet position at regular intervals 
 * and save .ppm images of the simulation 
 */
event snapshot (t=0; t+=checkat; t<=duration) {
    
    // We dump a snapshot which can be used to restart the simulation
    char dump_name[80];
    sprintf(dump_name, "./results/dump-%i.ppm", (int)(t*fabs(9.81/Dd)/Dd));
    dump (file = dump_name);
  
    char ppm_name_vof[80];
    sprintf(ppm_name_vof, "./results/vof-%i.ppm", (int)(t*fabs(9.81/Dd)/Dd));
    
    char ppm_name_level[80];
    sprintf(ppm_name_level, "./results/level-%i.ppm", (int)(t*fabs(9.81/Dd)/Dd));
    
    
    /*
    * VoF cut-off
    */
    scalar ff1[];
    scalar ff3[];
    foreach() {
      ff1[] = f1[] < 1e-4 ? 0 : f1[] > 1. - 1e-4 ? 1. : f1[];
      ff3[] = f3[] < 1e-4 ? 0 : f3[] > 1. - 1e-4 ? 1. : f3[];
    }
    boundary({ff1, ff3});
    
    scalar ff[];
    foreach() {
        ff[] = ff1[] + 50*ff3[];
    }
    
    scalar ll[];
    foreach() {
        ll[] = level;
    }
    boundary({ll});
    
    // Cells for which m is negative will be black in the movie.  
    scalar m[];
    foreach()
        m[] = 0.5 - particle[];
    boundary ({m});

    output_ppm(ff, file=ppm_name_vof, mask=m, n=1<<LEVEL_MAX);
//     output_ppm(ff, file=ppm_name_vof, n=1<<LEVEL_MAX);
    output_ppm(ll, file=ppm_name_level, max=LEVEL_MAX, mask=m, n=1<<LEVEL_MAX);
    //output_ppm (omega, fp, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m, min=-10, max=10, linear=true);
    
    
// #if _OPENMP || _MPI
//     char ppm_name_pid[80];
//     sprintf(ppm_name_pid, "./results/pid-%i.ppm", (int)(t*fabs(9.81/Dd)/Dd));
//     scalar pid[];
//     foreach()
//         pid[] = tid();
//     double tmax = npe() - 1;
//     output_ppm (pid, file=ppm_name_pid, max = tmax, n=1<<LEVEL_MAX);
// #endif
    
    
}


    