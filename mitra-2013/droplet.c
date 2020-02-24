/* 
 * This Basilisk setup is intended to simulate the droplet-particle collision cases described in:
 * 
 * Chemichal Engineering Science 100 (2013) 105-119
 * Droplet impact dynamics on a spherical particle 
 * Subhasish Mitra, Mayur J. Sathe, Elham Doroodchi, Ranjeet Utikar, Milin K. Shah, Vishnu Pareek, Jyeshtharaj B. Joshi c, Geoffrey M. Evans
 * 
 * Objective:
 *    - Subcooled droplet impact on a highly thermally conductive spherical surface was investigated both theoretically and experimentally.
 *    - The droplet spreading patterns in cold condition and film boiling regime were simulated using the 3D CFD models.
 *  
 * Analyses: 
 *    - The effect of Weber number on spreading of droplets of three different liquids namely water, isopropyl alcohol and acetone.
 *    - The droplet shape evolution and surface wetting upon droplet impact at surface temperatures ranging between 20 1C and 250 1C were investigated using a high speed camera. * 
 *    
 * Validation:
 *    - Maximum droplet spread was measured and compared with available correlations.
 *    - Despite a very small temperature drop in the film boiling regime indicating small fraction of vaporization, Schlieren imaging of acetone droplets showed qualitative vapour field around the rebounding droplets.
 *  
 * Outcomes: 
 *    - Generally wetting contact was observed at surface temperatures below or close to saturation temperature whilst a non-wetting contact was exhibited at surface temperatures significantly greater than the saturation temperature. The drop in surface temperature was found to be significantly lower in this non-wetting contact regime which led to significant reduction in heat transfer coefficient.
 * 
 */

/*
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
 * On cluster:
 * 
 * qcc -source droplet.c
 * mpicc -O2 -Wall -std=c99 -D_MPI=1 _droplet.c -o droplet -lm
 * mpirun -np 2 ./droplet
 * 
 */


// Models
// @paper: a volume of fluid (VOF) approach commercial solver ANSYS Fluent (version: 14)
// @paper: the vapour layer in the film boiling regime has not been modelled in the current study (see Section 4 for detail explanation). 
//         instead, the simplified approach proposed by Karl et al. (1996) was used to simulate film boiling regime cases
// @paper: second order upwind scheme was used for discretization of momentum and energy equation
//         volume fraction was discretized using Geo-Reconstruct scheme
//         pressure was discretized using Presto scheme while pressure–velocity coupling was obtained by SIMPLE scheme
// @paper: a first order implicit time stepping method was used in all the simulations


// @basilisk: We use the centered Navier–Stokes solver, two-phase flow and the momentum-conserving option. Note that the momentum-conserving option is crucial to obtain stable solutions for this air-water density ratio configuration.
#include "navier-stokes/centered.h"
#include "navier-stokes/conserving.h"
#include "two-phase.h"
#include "tension.h"
#include "tracer.h"

#include "elementary_body.h"


// Cases
// @paper: two different cases were simulated
//         impact of water and isopropyl alcohol droplet at 20 °C and at 250 °C
//         case 1 - softly deposited water droplet (height of fall <10 mm) on the surface of the particle was studied at a surface temperature of 100 °C,
//                  aiming at measuring evaporation time

#define CASE 1


// Temperature
// @paper: setup with a temperature controlled test piece, at temperatures ranging between 20 °C and 250 °C
// @paper: a 12 mm x 20 mm size hollow aluminium block was used as a heating billet with a 200 W cylindrical cartridge heater placed inside the block
// @paper: a steady state surface temperature was ensured during the experiments before the data were collected

// Experimental particle surface temperature (°C) x time (s)
#define Tp(t) ( 99.97734 - 5.627058e-1*pow(t,1) + 6.454108e-2*pow(t,2) - 7.905534e-3*pow(t,3) + 6.935052e-4*pow(t,4) - 3.706304e-5*pow(t,5) + 1.18456e-6*pow(t,6) - 2.270449e-8*pow(t,7) + 2.556086e-10*pow(t,8) - 1.557908e-12*pow(t,9) + 3.967671e-15*pow(t,10) )


// Dimensionless parameters
// @paper: examine the effect of Weber number on the droplet deformation and evaporation
// @paper: the Weber numbers in the experiments were varied by changing the height of the needle relative to the heated surface ranging from 10 mm to 100 mm
// @paper: Weber number range of 8–84, 14–136 and 13–130
//#define Re 1000
#define We 8



// Fluid properties
// @paper: distilled water, isopropyl alcohol and acetone

//define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)

#define RHO_L 801.16                   
#define RHO_G 0.872                    
#define RHO_S 400.                     

#define MU_L RHO_L*U0*D0/Re            
#define MU_G 3.74e-5                 

#define SIGMA RHO_L*U0*U0*D0/We






// Physical and computational domain
// @paper: domain dimensions were decided based on the experimental observations of maximum spread diameter of the impinging droplet
//         a domain size of 12mm x 9mm x 12mm was used for 3D simulations
#define LARGER 1.5
#define DOMAIN_LENGTH LARGER*8*D0



// Particle
// @paper: a 10e-3 m diameter brass particle was placed on top of the heating billet
// @paper: the polished brass particles used in the experiments had surface roughness of 0.12e-6 m
// @paper: before each experiment the brass particle was thoroughly cleaned with soft liquid soap followed by isopropyl alcohol and acetone rinsing 
//         to remove any dirt or greasy substances from the surface
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
// @paper: the height of the needle was varied to change the Weber number of the impinging droplet
// @paper: three different fluids—water, isopropyl alcohol and acetone
// @paper: a 15.6 ml water droplet with corresponding diameter of 3.1e-3
#define Dd 3.1e-3
#define DROPLET_PARTICLE_GAP 0.1*Dd
#define Xd 0.
const double Ydf = Yp + 0.5*Dp + 0.5*Dd + DROPLET_PARTICLE_GAP;
#define Yd (Yp + 0.5*Dp + 0.5*Dd + DROPLET_PARTICLE_GAP)
#define Ud 0.
#define Vd 15.
#if dimension == 3
#define Zd 0.
#define Wd 0.
#endif


// Spatial and temporal resolution
// @paper: a CMOS type high speed digital camera (Dantec IDT) was used to capture the images of droplet impact on the particle
// @paper: the typical frame rate and shutter speeds used in the experiments were 1000–1400 frames/s and 50–500 ms
// @paper: all simulations were performed using a time step of 10^-6 with 20–30 iterations per timestep
// @paper: a hexahedral mesh
//         total 321,538 hexahedral cells were used, which was found to calculate the droplet volume with less than 1% deviation
//         cells were dense near the spherical surface and were coarse gradually towards the boundaries for better resolving of the interface near the solid surface
//         total 10,592 cells were patched to model the droplet in 3D mesh. For 3D simulation of isopropyl alcohol droplet, 4,412 cells were patched
#define REFERENCE_FRAMERATE 1/1400 // 7.14e-4 s



// Outputs
// @paper: the camera was placed parallel to the impact surface at a zero degree angle to the setup
// @paper: three different imaging methodologies were carried out to study the droplet surface interaction:
//         - the illumination of the light was made uniform in the images by placing white paper diffusers behind the impacting object when 
//           (droplet interaction dynamics on the surface) was of interest
//         - the shadowgraphy technique was applied by placing a plastic diffuser screen between the test setup and the light source when high contrast images were required for
//           (measuring maximum spread of the droplets on the surface)
//         - Schlieren imaging was carried out in order 
//           (to visualize the vapour trail around the impinging acetone droplets on the hot surface)



// 

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



// Boundary condition
// @paper: - no-slip boundary condition was applied at the particle wall in the cold simulation cases while 
//         a free slip condition was applied on the surface in the film boiling regime cases
// @paper: - a contact angle was applied as a wall boundary condition
//         - for simulating isopropyl alcohol droplet impact, a static contact angle of 0 °C was used
//         in case of water, both static contact angle (90°) and an experimentally measured dynamic contact angle profile was used
//         - using static contact angle, the simulations cannot fully capture the dynamics of the system
//         CFD modelling using a constant contact angle delays in capturing this recoiling phase
// @paper: - pressure inlet boundary condition was applied on all the surrounding faces
//         - a free-slip condition was applied on the wall along with a constant contact angle of 180° to model non-wetting behaviour of the vapour layer


#define PRESSURE 101000

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



// Outcome - maximum spread diameter, peak and dip time (Fig. 4, Fig. 5)
// @paper: the free surface of the droplet deforms continuously due to conversion of impact kinetic energy into surface energy and viscous dissipation
// @paper: the maximum spread diameter which occurs in between 4.3 ms and 7.1 ms
//         the droplet takes shape of a doughnut at this phase
//         the peak of droplet visible at 4.3 ms takes a dip at 7.1 ms initiating the recoil phase
//         the model with static contact angle boundary conditions predicts the maximum spread at a later time than that of the model with the dynamic contact angle boundary condition


// Outcome - oscilattion period
// @paper: oscillation period can be estimated from below expression proposed by Schiaffino and Sonin (1997)
//         t_osc = sqrt( Dd^3.rho_l/sigma )
// @paper: using above expression, oscillation time was found to be 20.6 ms which was found to be 30% higher than 
//         the experimentally observed oscillation period of first spreading and recoiling phase



// Outcome - average evaporation time (Fig. 6)
// @paper: softly deposited water droplet (height of fall <10 mm)
// @paper: surface temperature of 100 °C
// @paper: the average evaporation time for the water droplet at the surface temperature of 100 °C was measured to be 30.27 s


scalar vapor[], temperature[];
scalar * tracers = {vapor, temperature};
face vector uf_save[];

int main() {

    printf("\nYpf: %1.6f\nYp: %1.6f\n", Ypf, Yp);
    
    printf("\nMU_L: %1.2e\nSIGMA: %1.2e\n", MU_L, SIGMA);
    
    
    // Domain length    
    L0 = DOMAIN_LENGTH;
            
     
    // Left-bottom origin at
    origin(-0.5*DOMAIN_LENGTH, 0);
       
    // TODO physical properties = f(temperature)
 
    // Density of phases 1-Liquid 2-Gas 

    rho1 = RHO_L;
    rho2 = RHO_G;    
        
    
    // Dynamic viscosity of phases 1-Liquid 2-Gas 
    mu1 = MU_L;
    mu2 = MU_G;
       
    
    // Surface tension coefficient
    f.sigma = SIGMA;
        
    
    // @basilisk.fr: decreasing the tolerance on the Poisson solve improves the results
    // @paper: a residual of 10^-4 was set for convergence of continuity, momentum and volume fraction equations while
    //         a residual of 10^-6 was used for energy equation
    TOLERANCE = 1e-4;
        
    
    
    // Initial uniform mesh
    N = 1 << LEVEL;
    init_grid(N);
    
    run(); 
  
}





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


    