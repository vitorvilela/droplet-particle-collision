/* 
 * This Basilisk setup is intended to simulate the droplet evaporation case described in:
 * 
 * The Journal of Physical Chemistry, 2017 
 * Marangoni Flow Induced Evaporation Enhancement on Binary Sessile Drops
 * Pin Chen, Souad Harmand, Safouene Ouenzerfi, Jesse Schiffler
 * 
 * Objective:
 *    - The evaporation processes of pure water, pure 1-butanol, and 5% 1-butanol aqueous solution drops on heated hydrophobic substrates are investigated 
 *      to determine the effect of temperature on the drop evaporation behavior.
 *    - Another objective of this study is to determine the relation between the substrate temperature and the Marangoni effect as well as their influence on the evaporation rate.
 *  
 * Analyses: 
 *    - The evolution of the parameters (contact angle, diameter, and volume) during evaporation.
 *    - Infrared thermal mapping of the drop surface. 
 *    
 * Validation:
 *    - A series of empirical equations for predicting the evaporation rates.
 *  
 * Outcomes: 
 *    - The pure 1-butanol drop does not show any thermal instability at different substrate temperatures, 
 *      while the convection cells created by the thermal Marangoni effect appear on the surface of the pure water drop from 50°C.
 * 
 */

/* Introduction:
 * 
 *    - Applications: coatings [1], ink printing [2,3], and cooling systems [4]
 *    - Related experimental [5-9], analytical [6,10-12] and numerical [13-15] works
 *    - Regimes: constant contact angle regime, constant contact radius regime, and a combination of both regimes
 * 
 */

/* Phenomenology:
 * 
 *    - When a pure solution drop is deposited on a heated substrate, the temperature at the bottom of the drop is close to that of the substrate 
 *      while the temperature at the top of the drop is similar to the ambient temperature.
 *    - When the temperature difference as well as the surface tension difference between the top and the periphery of the drop becomes suitably large, 
 *      the liquid will flow from the low surface tension place to the high surface tension place on the drop surface (thermal Marangoni flow).
 *    - The surface tension difference can be also created by non-uniformity in composition and the generated flow in this case is termed as the solutal Marangoni flow.
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
 * 
 * Evaporation references at Basilisk:
 * 
 * http://basilisk.fr/sandbox/qmagdelaine/phase_change/1_elementary_body/static_drop.c#physical-parameters
 * http://basilisk.fr/sandbox/qmagdelaine/phase_change/1_elementary_body/static_drop.c#physical-parameters
 * 
 */


// Models / Methods
// @paper: 


// @basilisk: We use the centered Navier–Stokes solver, two-phase flow and the momentum-conserving option. Note that the momentum-conserving option is crucial to obtain stable solutions for this air-water density ratio configuration.

#include "navier-stokes/centered.h"
#include "navier-stokes/conserving.h"
#include "two-phase.h"
#include "tension.h"
#include "tracer.h"

#include "elementary_body.h"


// Cases
// @paper: Ts = 22°C, 35°C, 50°C, 60°C, 70°C, and 80°C
// @paper: The diameter of tested drop is always smaller than the capillary length (2.7 mm at standard temperature and pressure)
// @paper: pure water, pure 1-butanol, and 5% 1-butanol aqueous solution

// Case 0: 
// Ts = Ta = 22°C + 273.15 [K], relative humidity 50%
// Fluid: pure water
// Droplet initial diameter: Dd = 1.3e-3
// Gravity: 0.
#define CASE0


// Temperature
#define Ta (22+273.15)
#define Ts (22+273.15)
#define RH 0.5
#define GRAVITY 0.


// Fluid properties
// @paper: pure water, pure 1-butanol, and 5% 1-butanol aqueous solution drops
// @paper: ambient temperature and relative humidity can be controlled
// @web.colby.edu/ch217public/files/2012/04/density-of-water.xlsx: rho [kg/m³], T [°C] .:. [K]-->[°C]
// @http://ddbonline.ddbst.de/VogelCalculation/VogelCalculationCGI.exe?component=Water: mu [Pa.s], T [K]
// @https://en.wikipedia.org/wiki/Density_of_air#Humid_air: rho [kg/m³] for humid air, T [°C] .:. [K]-->[°C] and [K], when indicated
// @http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/fprops/propsoffluids/node5.html: mu [kg/m.s]
// @white, 2006: sigma [J/m²], T [°C] .:. [K]-->[°C]

// Partial pressure of water vapour [hPa = 100 Pa]
#define psat (0.01 * 6.1078*pow(10.;7.5*(T-273.15)/(T-273.15)+237.3))

// Partial pressure of water vapour [Pa]
#define pv RH*psat

// Absolute pressure
#define PRESSURE 101000

// Partial pressure of dry air [Pa], where p stands for the absolute pressure
#define pd (PRESSURE-pv)

#define rho(f, T) (clamp(f,0,1)*(rho1(T) - rho2(T)) + rho2(T))
#define rho1(T) (1000 * (1 - ((T-273.15)+288.9414)/(508929.2*((T-273.15)+68.12963))*((T-273.15)-3.9863)^2))
#define rho2(T) (pd/(287.05*T) + pv/(461.495*T))

#define mu(f, T) (clamp(f,0,1)*(mu1(T) - mu2(T)) + mu2(T))
#define mu1(T) (0.001 * exp(-3.7188 + 578.919/(-137.546+T) ))
#define mu2(T) (1.458e-6*pow(T;3/2)/(T+110.4))                

#define SIGMA(T) (0.076 - 0.00017*(T-273.15))



// Physical and computational domain
// @paper: The substrates are placed in a vapor chamber (14 cm × 12.4 cm × 7.5 cm)

#define DOMAIN_LENGTH 14.e-3


// Wall / Substrate
// @paper: The drops are deposited on the hydrophobic silicon substrate, below which a heater connected to an electrical controller is attached
// @paper: six different substrate temperatures from Ts = 22 to 80°C



// Droplet
// @paper: "approximated by 1.3 mm from Figure 2. 

#define Dd 1.3e-3
#define Xd 0.
#define Yd 0.
#define Ud 0.
#define Vd 0.
#if dimension == 3
#define Zd 0.
#define Wd 0.
#endif


// Spatial and temporal resolution
// @paper: The infrared camera 640 × 512 pixels, 200Hz, 15 µm detector pitch
// @paper: Camera (Allied Vision Technologies, 15 fps, 780 × 580 pixels)

#define REFERENCE_FRAMERATE 1/15
#define FRAMERATE REFERENCE_FRAMERATE/2

// Refinement level (level, n): (4, 16) (5, 32) (6, 64) (7, 128) (8, 256) (9, 512) (10, 1024)
int LEVEL = 6;
int LEVEL_MAX = 10;






/*
 * Checking time and simulation duration
 */
#define nondimensional_checkat 0.1
#define nondimensional_duration 3.0 // 1.7 + (0.9 of non-dimensional impact time)
const double checkat = nondimensional_checkat*Dd/Vd; 
const double duration = nondimensional_duration*Dd/Vd;





// Evaporation parameters
// @paper: the ambient temperature in the vapor chamber was stabilized at Ta = 22°C, and the relative humidity was maintained at 50% throughout the experiments
// @basilisk: Typically, for the problem of a water droplet evaporating in dry air, this Peclet number is O(10−5), meaning that the problem is really dominated by diffusion
// Bolz and Tuve (1976) referred at docs/chapter91.pdf: D [m²/s], T [K]

#define F_ERR 1e-10

// Peclet number Pe = cs/rho2
#define peclet cs/rho2(T)

// Mass diffusion coefficient of vapour in air
#define D_V (-2.775e-6 + 4.479e-8*T + 1.656e-10*T*T)
// Saturation concentration of water in air [kg/m³]
// #define cs 0.5
// Vapour concentration in the droplet
#define vcs 1.
// Vapour concentration in air
#define cinf RH

// #define temperature_peclet 1e-3
// #define D_T 0.005
#define tcs Ta
#define tinf Ta
// #define tp 1.

#define dirichlet_time_factor 10.

// TODO D_V, D_T = f(x,y)



// Boundary condition

temperature[bottom] = dirichlet(Ts);
vapor[bottom] = neumann(0.);
u.n[bottom] = dirichlet(0.);
p[bottom] = neumann(0.);

temperature[top]   = dirichlet(Ta);
vapor[top] = neumann(0);
u.n[top] = dirichlet(0.);
p[top] = neumann(0.);

temperature[right] = dirichlet(Ta);
vapor[right] = neumann(0.);
u.n[right] = dirichlet(0.);
p[right] = neumann(PRESSURE);

temperature[left] = dirichlet(Ta);
vapor[left] = neumann(0.);
u.n[left] = dirichlet(0.);
p[left] = neumann(0.);



// Outcome - describe
// @paper: 



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
        
    
    // @basilisk: decreasing the tolerance on the Poisson solve improves the results
    
    TOLERANCE = 1e-4;
            
    
    // Initial uniform mesh
    N = 1 << LEVEL;
    init_grid(N);
    
    run(); 
  
}



// Initial conditions
// @paper: without air current
// @paper:

event init(i = 0) {
  
  if (!restore (file = "dump")) {  

    // Single refinement statment to particle and droplet
    refine ( (sq(x-Xd) + sq(y-Yd) < sq(1.1*0.5*Dd)) && level < LEVEL_MAX );   
    
    
    // Droplet
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



// TODO flow field = f(temperature source term)

/*
 * Add the acceleration of gravity in the downward (-y) direction
 */
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 0.;  
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



// Output
// @paper: The infrared camera (FLIR X6580SC, 640 × 512 pixels, 200Hz, 15 µm detector pitch) is installed on top for infrared thermal mapping and 
//         visualization of thermal instabilities on the surface of the droplets.

// Output
// @paper: CCD (charge-coupled device) camera (Allied Vision Technologies, 15 fps, 780 × 580 pixels) is used to record the evaporation process of the droplets for profile analysis 
//         to measure the contact angle, volume, diameter, and height of the sessile droplets during evaporation 



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


/* References:
 * 
 * (1) Kimura, M.; Misner, M. J.; Xu, T.; Kim, S. H.; Russell, T. P. Long-range ordering of diblock copolymers induced by droplet pinning. Langmuir 2003, 19(23), 9910-9913.
 * (2) De Gans, B. J.; Duineveld, P. C.; Schubert, U. S. Inkjet printing of polymers: state of the art and future developments. Adv. Mater. 2004, 16, 203–213.
 * (3) Park, J.; Moon, J. Control of colloidal particle deposit patterns within picoliter droplets ejected by ink-jet printing. Langmuir 2006, 22(8), 3506-3513.
 * (4) Su, X.; Zhang, M.; Han, W.; Guo, X. Enhancement of heat transport in oscillating heat pipe with ternary fluid. Int. J. Heat Mass Transf. 2015, 87, 258-264.
 * (5) Erbil, H. Y.; McHale, G.; Newton, M. I. Drop evaporation on solid surfaces: constant contact angle mode. Langmuir 2002, 18(7), 2636-2641.
 * (6) Nguyen, T. A. ; Nguyen, A. V. ; Hampton, M. A. ; Xu, Z. P.; Huang, L. ; Rudolph, V. Theoretical and experimental analysis of droplet evaporation on solid surfaces. Chem. Eng. Sci. 2012, 69(1), 522-529.
 * (7) Kelly-Zion, P. L.; Pursell, C. J.; Vaidya, S.; Batra, J. Evaporation of sessile drops under combined diffusion and natural convection. Colloids Surf., A. 2011, 381(1), 31-36.
 * (8) Bennacer, R.; Sefiane, K. Vortices, dissipation and flow transition in volatile binary drops. J. Fluid Mech. 2014, 749, 649-665.
 * (9) Shi, L.; Shen, P.; Zhang, D.; Lin, Q.; Jiang, Q. Wetting and evaporation behaviors of water–ethanol sessile drops on PTFE surfaces. Surf. Interface Anal. 2009, 41(12‐ 13), 951-955.
 * (10) Masoud, H.; Felske, J. D. Analytical solution for inviscid flow inside an evaporating sessile drop. Phys. Rev. E: Stat. Phys., Plasmas, Fluids,. 2009, 79(1), 016301.
 * (11) Sefiane, K.; Bennacer, R. An expression for droplet evaporation incorporating thermal effects. J. Fluid Mech. 2011, 667, 260-271.
 * (12) Bormashenko, E. Young, Boruvka–Neumann, Wenzel and Cassie–Baxter equations as the transversality conditions for the variational problem of wetting. Colloids Surf., A. 2009, 345(1), 163-165.
 * (13) Korlie, M. S. Three-dimensional computer simulation of liquid drop evaporation. Comput. Math. Appl. 2000, 39(12), 43-52.
 * (14) Widjaja, E.; Harris, M. T. Numerical study of vapor phase-diffusion driven sessile drop evaporation. Comput. Chem. Eng. 2008, 32(10), 2169-2178.
 * (15) Barash, L. Y.; Bigioni, T. P.; Vinokur, V. M.; Shchur, L. N. Evaporation and fluid dynamics of a sessile drop of capillary size. Phys. Rev. E. 2009, 79(4), 046301.
 * (16) Picknett, R. G.; Bexon, R. The evaporation of sessile or pendant drops in still air. J. Colloid Interface Sci. 1977, 61(2), 336-350.
 * (17) Rowan, S. M.; Newton, M. I.; McHale, G. Evaporation of microdroplets and the wetting of solid surfaces. J. Phys. Chem. 1995, 99(35), 13268-13271.
 * (18) Shanahan, M. E. R.; Bourges, C. Effects of evaporation on contact angles on polymer surfaces. Int. J. Adhes. Adhes. 1994, 14(3), 201-205.
 * (19) Bourges-Monnier, C.; Shanahan, M. E. R. Influence of evaporation on contact angle. Langmuir 1995, 11(7), 2820-2829.
 * 
 */


    