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
 *      with MPI: mpirun -np 2 ./droplet
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
 * Fuel Processing Technology. 
 * Numerical investigation of heavy fuel droplet-particle collisions in the injection zone of a 
 * Fluid Catalytic Cracking reactor, Part I: Numerical model and 2D simulations.
 * Manolis Gavaises, 2016.
 *  
 */

#include "navier-stokes/centered.h"
#include "three-phase.h"
#include "tension.h"

/* 
 * Non-dimensional numbers 
 * Case #1
 */

#define Re 2272
#define We 4266

/*
 * Length scale
 */
#define D0 75.e-6
#define U0 15. 

/*
 * Fluid flow properties
 */
#define RHO_L 801.16                    // from paper
#define RHO_G 0.872                     // from peacesoftware.de
#define RHO_S 400.                      // any intermediate

#define MU_L RHO_L*U0*D0/Re             // from paper f(Re) 3.97e-04       // 5.03e-5 from chemeo.com 
#define MU_G 3.74e-5                    // from peacesoftware.de
#define MU_S 1.e-1                      // any big

// SIGMA_12=SIGMA_1+SIGMA_2             // from paper f(We) 3.17e-03       // 0.0283 from api datashet
#define SIGMA_1 RHO_L*U0*U0*D0/We       // Liquid
#define SIGMA_2 0.00001                 // Gas
#define SIGMA_3 0.9                     // Solid

#define PRESSURE 202650                 // from paper





/*
 * Particle
 */
#define Dp D0
#define Xp 0.
#define Yp 5.6*D0


/*
 * Droplet
 */
#define Dd D0
#define Xd 0.
const double Yd = Yp + 0.5*Dp + Dd;
#define Ud 0
#define Vd 15.


/* 
 * Domain lenght
 */
#define DOMAIN_LENGTH 12*D0


/*
 * Checking time and simulation duration
 */
#define nondimensional_checkat 0.1
#define nondimensional_duration 2.2 // 1.7 + (0.5 of non-dimensional impact time)
const double checkat = nondimensional_checkat*Dd/Vd; 
const double duration = nondimensional_duration*Dd/Vd;


/*
 * Refinement level
 * (level, n): (4, 16) (5, 32) (6, 64) (7, 128) (8, 256) 
 */
int LEVEL = 5;
int LEVEL_MAX = 10;

double volref[3];

int main() {

    printf("\nMU_L: %1.2e\nSIGMA: %1.2e\n", MU_L, SIGMA_1);
    
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
    mu2 = MU_G;
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




event init(i = 0) {
    
  if (!restore (file = "dump")) {  

    // Particle
    refine ( ( sq(x-Xp) + sq(y-Yp) < sq(1.1*0.5*Dp) ) && level < LEVEL_MAX);
    fraction ( f3, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );
    
    // Droplet
    refine ( ( sq(x-Xd) + sq(y-Yd) < sq(1.1*0.5*Dd) ) && level < LEVEL_MAX);
    fraction ( f1, sq(0.5*Dd) - sq(x-Xd) - sq(y-Yd) );
        
    foreach() {     
        u.x[] = f1[]*Ud;
        u.y[] = f1[]*(-Vd);    
        f2[] = 1.0 - f1[] - f3[];
    }
    boundary ((scalar *){u}); 
    
    volref[0] = statsf(f1).sum;
    volref[1] = statsf(f2).sum;
    volref[2] = statsf(f3).sum;
    
  }
}


event stationary_particle (i++) {
  coord vp = {0., 0.}; // the velocity of the particle
  // Moving particle
  //fraction (particle, - (sq(x - vp.x*t) + sq(y - vp.y*t) - sq(0.5*Dp)));
  // Stationary particle
  fraction ( f3, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );
  
  // We then use this (solid) volume fraction field to impose the corresponding velocity (momentum) field inside the solid (as a volume-weighted average).
  
  foreach()
    foreach_dimension()
      u.x[] = f3[]*vp.x + (1. - f3[])*u.x[];
  boundary ((scalar *){u});
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
  adapt_wavelet({f1, f3, u}, (double[]){1.e-3, 1.e-3, 1.e-2, 1.e-2}, LEVEL_MAX);
  
//   unrefine (y<(1/2)*DOMAIN_LENGTH); // || y>(7/8)*DOMAIN_LENGTH); // || fabs(x)>(7/8)*DOMAIN_LENGTH);
  
}


/*
 * Every ten timesteps, we output the time, position, velocity of
 * the droplet, timestep etc... 
 */
event logfile (i += 10) {
  /*
  double xb = 0., yb = 0., sb = 0.;
  double vbx = 0., vby = 0.;
    
  vector u[];
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  boundary((scalar *){u});
  
  foreach(reduction(+:xb) reduction(+:yb) 
	  reduction(+:vbx) reduction(+:vby) 
	  reduction(+:sb)) {
    
      double dv = (1. - f[])*dv();
  
      xb += x*dv;
      yb += y*dv;    
      vbx += u.x[]*dv;
      vby += u.y[]*dv;    
      sb += dv;
    
  }
  
  droplet_track = (yb/sb)/D;
  droplet_u = vbx/fabs(sb);
  droplet_v = vby/fabs(sb);
  printf("Droplet-X = %f\n", xb/fabs(sb));
  printf("Droplet-Y / D = %f\n", droplet_track);  
  printf("At velocity u(%f) v(%f) \n", droplet_u, droplet_v);
  */
  fprintf (fout,
	   "t: %1.1e  dt: %1.1e \n", 
	   t, dt);
  fflush (fout);
  
}


/*
 * Check droplet position at regular intervals 
 * and save .ppm images of the simulation 
 */
event snapshot (t=0; t+=checkat; t<=1.1*duration) {
    
    // We dump a snapshot which can be used to restart the simulation
    dump (file = "./results/dump");
  
    char ppm_name_vof[80];
//     sprintf(ppm_name_vof, "vof-%1.f.ppm", t*1000);
    sprintf(ppm_name_vof, "./results/vof-%1.2f.ppm", t*Vd/Dd - 0.5); // tal = 0.5 at impact
        
    char ppm_name_level[80];
//     sprintf(ppm_name_level, "level-%1.f.ppm", t*1000);
    sprintf(ppm_name_level, "./results/level-%1.2f.ppm", t*Vd/Dd - 0.5);
    
    
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
        ff[] = ff1[] + 10*ff3[];
    }
    
    scalar ll[];
    foreach() {
        ll[] = level;
    }
    boundary({ll});
    
//     // Cells for which m is negative will be black in the movie.  
//     scalar m[];
//     foreach()
//         m[] = 0.5 - particle[];
//     boundary ({m});

//     output_ppm(ff, file=ppm_name_vof, mask=m, n=1<<LEVEL_MAX);
    output_ppm(ff, file=ppm_name_vof, n=1<<LEVEL_MAX);
    output_ppm(ll, file=ppm_name_level, max=LEVEL_MAX, n=1<<LEVEL_MAX);
    //output_ppm (omega, fp, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m, min=-10, max=10, linear=true);
    
    
// #if _OPENMP || _MPI
//     char ppm_name_pid[80];
// //     sprintf(ppm_name_pid, "pid-%1.f.ppm", t*1000);
//     sprintf(ppm_name_pid, "./results/pid-%1.2f.ppm", t*fabs(Vd)/Dd);
//     scalar pid[];
//     foreach()
//         pid[] = tid();
//     double tmax = npe() - 1;
//     output_ppm (pid, file=ppm_name_pid, max = tmax, n=1<<LEVEL_MAX);
// #endif
    
    
}


    