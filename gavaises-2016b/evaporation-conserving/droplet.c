/*
 * Droplet-Particle collision - Immersed Boundary imposed u - Non-conservative approach
 *
 * Compile with:
 * 
 *      3D with MPI: CC99='mpicc -std=c99' qcc -Wall -O2 -grid=octree -D_MPI=1 droplet.c -o droplet -lm
 *      2D with MPI: CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 droplet.c -o droplet -lm
 *      2D without MPI: qcc -Wall -O2 droplet.c -o droplet -lm
 * 
 * Run with:
 * 
 *      with MPI: mpirun -np 2 droplet
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

//#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tracer.h"

#include "elementary_body.h"

/* 
 * Non-dimensional numbers 
 * Case #1
 */

#define Re 1000 //2272
#define We 100 //4266

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

// We have two phases e.g. air and water. For large viscosity and density ratios, the harmonic mean for the viscosity tends to work better than the default arithmetic mean. We “overload” the default by defining the mu() macro before including the code for two phases.
//#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))


#define MU_L RHO_L*U0*D0/Re             // from paper f(Re) 3.97e-04       // 5.03e-5 from chemeo.com 
#define MU_G 3.74e-5                    // from peacesoftware.de

#define SIGMA RHO_L*U0*U0*D0/We		// Liquid 3.17e-03 // air-water: 0.0728

#define PRESSURE 202650  		// from paper

#define LARGER 1.5 // 1.5 //1.5
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
#define nondimensional_duration 3.0 // 1.7 + (0.9 of non-dimensional impact time)
const double checkat = nondimensional_checkat*Dd/Vd; 
const double duration = nondimensional_duration*Dd/Vd;



// Evaporation parameters
#define F_ERR 1e-10

#define vapor_peclet 1e-3
#define D_V 0.01 // 1
#define vcs 1.
#define cinf 0.2 //1. //0.2

#define temperature_peclet 1e-3
#define D_T 0.005
#define tcs 0.
#define tinf 0.2
#define tp 1.

#define dirichlet_time_factor 10.


// TODO D_V, D_T = f(x,y)

/*
 * Refinement level
 * (level, n): (4, 16) (5, 32) (6, 64) (7, 128) (8, 256) 
 */
int LEVEL = 6;
int LEVEL_MAX = 10;

scalar vapor[], temperature[];
scalar * tracers = {vapor, temperature};
face vector uf_save[];


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
       
    // TODO physical properties = f(temperature)
    
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
temperature[bottom] = dirichlet(tinf);
vapor[bottom] = neumann(0);
u.n[bottom] = neumann(0);
p[bottom] = dirichlet(PRESSURE);

temperature[top]   = dirichlet(tinf);
vapor[top] = neumann(0);
u.n[top] = neumann(0);
p[top] = dirichlet(PRESSURE);

temperature[right] = dirichlet(tinf);
vapor[right] = neumann(0);
u.n[right] = neumann(0);
p[right] = dirichlet(PRESSURE);

temperature[left] = dirichlet(tinf);
vapor[left] = neumann(0);
u.n[left] = neumann(0);
p[left] = dirichlet(PRESSURE);


scalar particle[];

event init(i = 0) {
  
  if (!restore (file = "dump")) {  

    // Single refinement statment to particle and droplet
    refine ( ( (sq(x-Xp) + sq(y-Yp)< sq(1.1*0.5*Dp)) ||  (sq(x-Xd) + sq(y-Yd) < sq(1.1*0.5*Dd)) ) && level < LEVEL_MAX );   
    
    // Particle
    //refine ( (sq(x-Xp)+sq(y-Yp)<sq(1.1*0.5*Dp)) && (sq(x-Xp)+sq(y-Yp)>sq(0.9*0.5*Dp)) && level<LEVEL_MAX);
//     refine ( (sq(x-Xp)+sq(y-Yp)<sq(1.1*0.5*Dp)) && level<LEVEL_MAX);
    fraction ( particle, sq(0.5*Dp)-sq(x-Xp)- sq(y-Yp) );
    
    // Droplet
    //refine ( (sq(x-Xd)+sq(y-Yd)<sq(1.1*0.5*Dd)) && (sq(x-Xd)+sq(y-Yd)>sq(0.9*0.5*Dd)) && level<LEVEL_MAX);
//     refine ( (sq(x-Xd)+sq(y-Yd)<sq(1.1*0.5*Dd)) && level<LEVEL_MAX);
    fraction ( f, sq(0.5*Dd) - sq(x-Xd) - sq(y-Yd) );
	
    foreach() {      
      u.x[] = f[]*Ud;
      u.y[] = f[]*(-Vd);  
      vapor[] = f[]*vcs + (1. - f[])*cinf;
      temperature[] = f[]*tcs + (1. - f[])*tinf;
      vapor[] = particle[]*cinf + (1. - particle[])*vapor[]; 
      temperature[] = particle[]*tp + (1. - particle[])*temperature[]; 
    }
    boundary ((scalar *){u}); 
    boundary ({vapor, temperature});
    
  }
  
}

// TODO evaporation = f(temperature)

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
    int max_level = LEVEL_MAX;
  #else
    int max_level = LEVEL;
  #endif
  
  vapor.D = D_V;
  vapor.tr_eq = vcs;
  vapor.inverse = true;
  dirichlet_diffusion (vapor, f, max_level, dt, dirichlet_time_factor);
      
  temperature.D = D_T;
  temperature.peclet = temperature_peclet;
  temperature.inverse = true;
  dirichlet_diffusion (temperature, f, max_level, dt, dirichlet_time_factor);
  
}

event stationary_particle (i++) {
  
  coord vp = {0., 0.}; // the velocity of the particle
  // Moving particle
  //fraction (particle, - (sq(x - vp.x*t) + sq(y - vp.y*t) - sq(0.5*Dp)));
  // Stationary particle
  fraction ( particle, sq(0.5*Dp) - sq(x-Xp) - sq(y-Yp) );
  
  // We then use this (solid) volume fraction field to impose the corresponding velocity (momentum) field inside the solid (as a volume-weighted average).
  
  foreach() {
    vapor[] = particle[]*cinf + (1. - particle[])*vapor[]; 
    temperature[] = particle[]*tp + (1. - particle[])*temperature[];
    foreach_dimension()
      u.x[] = particle[]*vp.x + (1. - particle[])*u.x[];    
  }
  boundary ((scalar *){u});
  boundary ({vapor, temperature});
  
}

// TODO flow field = f(temperature source term)

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
  adapt_wavelet({f, particle, u}, (double[]){1.e-3, 1.e-3, 1.e-2, 1.e-2}, LEVEL_MAX);
  //adapt_wavelet({f}, (double[]){1.e-4}, LEVEL_MAX);

  //   unrefine (y<(1/2)*DOMAIN_LENGTH); // || y>(7/8)*DOMAIN_LENGTH); // || fabs(x)>(7/8)*DOMAIN_LENGTH);  
}



event log_console (i++) {
  if (i == 0)
    fprintf (fout, "t dt mgp.i mgu.i grid->tn \n");
  fprintf (fout, "%g %g %d %d %ld\n", t, dt, mgp.i, mgu.i, grid->tn);
  fflush (fout);
}

/*
 * Every ten timesteps, we output the time, position, velocity of
 * the droplet, timestep etc... 
 */
// event logfile (i += 10) {
//   /*
//   double xb = 0., yb = 0., sb = 0.;
//   double vbx = 0., vby = 0.;
//     
//   vector u[];
//   foreach()
//     foreach_dimension()
//       u.x[] = q.x[]/rho[];
//   boundary((scalar *){u});
//   
//   foreach(reduction(+:xb) reduction(+:yb) 
// 	  reduction(+:vbx) reduction(+:vby) 
// 	  reduction(+:sb)) {
//     
//       double dv = (1. - f[])*dv();
//   
//       xb += x*dv;
//       yb += y*dv;    
//       vbx += u.x[]*dv;
//       vby += u.y[]*dv;    
//       sb += dv;
//     
//   }
//   
//   droplet_track = (yb/sb)/D;
//   droplet_u = vbx/fabs(sb);
//   droplet_v = vby/fabs(sb);
//   printf("Droplet-X = %f\n", xb/fabs(sb));
//   printf("Droplet-Y / D = %f\n", droplet_track);  
//   printf("At velocity u(%f) v(%f) \n", droplet_u, droplet_v);
//   */
//   fprintf (fout,
// 	   "t: %.8f  dt: %.8f \n", 
// 	   t, dt);
//   fflush (fout);
//   
// }






event snapshot (t=0; t+=checkat; t<=duration) {
  
    // We dump a snapshot which can be used to restart the simulation
  char dump_name[80];
//     sprintf(dump_name, "dump-%1.3f", t*1000);
  sprintf(dump_name, "dump-%1.1f", t*Vd/Dd); // - 0.9);
    dump (file = dump_name);
  
    char ppm_name_vof[80];
//     sprintf(ppm_name_vof, "vof-%1.3f.ppm", t*1000);
    sprintf(ppm_name_vof, "vof-%1.1f.ppm", t*Vd/Dd); // - 0.9); // tal = 0.9 at impact
        
    char ppm_name_level[80];
//     sprintf(ppm_name_level, "level-%1.3f.ppm", t*1000);
    sprintf(ppm_name_level, "level-%1.1f.ppm", t*Vd/Dd); // - 0.9);
    
    char ppm_name_vapor[80];
//       sprintf(ppm_name_vapor, "c-%1.3f.ppm", t*1000);
    sprintf(ppm_name_vapor, "c-%1.1f.ppm", t*Vd/Dd); // - 0.9);
    
    char ppm_name_temperature[80];
    sprintf(ppm_name_temperature, "temp-%1.1f.ppm", t*Vd/Dd);
 
    
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
    foreach() {
        m[] = 0.5 - particle[];       
    }
    boundary ({m});
    

    output_ppm(ff, file=ppm_name_vof, mask=m, n=1<<LEVEL_MAX);
    output_ppm(ll, file=ppm_name_level, max=LEVEL_MAX, n=1<<LEVEL_MAX);
    output_ppm(vapor, file=ppm_name_vapor, min=cinf, max=vcs, n=1<<LEVEL_MAX);
    output_ppm(temperature, file=ppm_name_temperature, min=tcs, max=tp, n=1<<LEVEL_MAX);
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


    