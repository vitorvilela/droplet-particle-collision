/*
 * Droplet-Particle collision - Immersed Boundary imposed Q - Conservative approach
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

#include "momentum.h"
#include "tension.h"
#define cf

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

#define SIGMA RHO_L*U0*U0*D0/We		// Liquid 3.17e-03

#define PRESSURE 202650  		// from paper

#define LARGER 1.5
#define HIGHER 1.2

/* 
 * Domain lenght
 */
#define DOMAIN_LENGTH LARGER*8*D0

/*
 * Particle
 */
#define Dp D0
#define Xp 0.
#define Yp HIGHER*5.6*D0

/*
 * Droplet
 */
#define Dd D0
#define Xd 0.
const double Yd = Yp + 0.5*Dp + Dd;
#define Ud 0
#define Vd 15.


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


int main() {

    printf("\nMU_L: %1.2e\nSIGMA: %1.2e\n", MU_L, SIGMA);
    
    /*
     * Domain Length
     */
    L0 = DOMAIN_LENGTH;
            
    /* 
     * Left-bottom origin at
     */
    origin(-0.5*DOMAIN_LENGTH, 0);
          
    /* 
     * Density of phases 1-Liquid 2-Gas 
     */
    rho1 = RHO_L;
    rho2 = RHO_G;    
        
    /*
     * Dynamic viscosity of phases 1-Liquid 2-Gas 
     */
    mu1 = MU_L;
    mu2 = MU_G;
        
    /*
     * Surface tension 
     */
    f.sigma = SIGMA;
        
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
q.n[bottom] = neumann(0);
p[bottom] = dirichlet(PRESSURE);

q.n[top] = neumann(0);
p[top] = dirichlet(PRESSURE);

q.n[right] = neumann(0);
p[right] = dirichlet(PRESSURE);

q.n[left] = neumann(0);
p[left] = dirichlet(PRESSURE);


scalar particle[];

event init(i = 0) {
  
  if (!restore (file = "dump")) {  

    // Particle
    refine ( ( sq(x-Xp) + sq(y-Yp) < sq(1.1*0.5*Dp) ) && level < LEVEL_MAX);
    fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );
    
    // Droplet
    refine ( ( sq(x-Xd) + sq(y-Yd) < sq(1.1*0.5*Dd) ) && level < LEVEL_MAX);
    fraction ( f, sq(0.5*Dd) - sq(x-Xd) - sq(y-Yd) );
	
    foreach() {      
      q.x[] = f[]*rho1*Ud;
      q.y[] = f[]*rho1*(-Vd);    
    }
    boundary ((scalar *){q}); 
    
#ifndef cf
    foreach()
      cf[] = f[];
#endif
    
  }
  
}

event stationary_particle (i++) {
  coord vp = {0., 0.}; // the velocity of the particle
  // Moving particle
  //fraction (particle, - (sq(x - vp.x*t) + sq(y - vp.y*t) - sq(0.5*Dp)));
  // Stationary particle
  fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );
  
  // We then use this (solid) volume fraction field to impose the corresponding velocity (momentum) field inside the solid (as a volume-weighted average).
  
  foreach()
    foreach_dimension()
      q.x[] = particle[]*RHO_S*vp.x + (1. - particle[])*q.x[];
  boundary ((scalar *){q});
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
   * Extract velocity from momentum
   */
  vector u[];
  foreach()
    foreach_dimension()
      u.x[] = q.x[]/rho[];
  boundary((scalar *){u});
  
  /*
   * Refinement threshold vector for f, u.x and u.y 
   */
  adapt_wavelet({f, u}, (double[]){1.e-3, 1.e-2, 1.e-2}, LEVEL_MAX);
  //adapt_wavelet({f}, (double[]){1.e-4}, LEVEL_MAX);

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
	   "t: %.8f  dt: %.8f \n", 
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
    
    // Cells for which m is negative will be black in the movie.  
    scalar m[];
    foreach()
        m[] = 0.5 - particle[];
    boundary ({m});
  
    
    vector u[];
    foreach()
      foreach_dimension()
        u.x[] = q.x[]/rho[];
    boundary((scalar *){u});

    output_ppm(ff, file=ppm_name_vof, mask=m, n=1<<LEVEL_MAX);
    output_ppm(ll, file=ppm_name_level, max=LEVEL_MAX, n=1<<LEVEL_MAX);
    //output_ppm (omega, fp, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m, min=-10, max=10, linear=true);
    
    
// #if _OPENMP || _MPI
//     char ppm_name_pid[80];
// //     sprintf(ppm_name_pid, "pid-%1.f.ppm", t*1000);
//     sprintf(ppm_name_pid, "./results/pid-%1.2f.ppm", t*Vd/Dd - 0.5); // tal = 0.5 at impact);
//     scalar pid[];
//     foreach()
//         pid[] = tid();
//     double tmax = npe() - 1;
//     output_ppm (pid, file=ppm_name_pid, max = tmax, n=1<<LEVEL_MAX);
// #endif
    
    
}


    