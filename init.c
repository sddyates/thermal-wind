/*================================================================================*/
/*
	This init file initialises an isothermal parker wind around an sun like
	star with a jutier like exoplanet with an orbital radius of 0.047 AU.
	Due to the close proximity of the planet to its hoast star the planets
	surface is being iradiated by intence UV radiation which induces a
	planetary wind. The planet, like the star, has a dipolar magnetosphere
	and hence the planetary wind is initilised in the same way as the stellar
	wind.

	The model is described in detail in the paper:
	Matsakos T., Uribe1 A., Königl A., 2015, A&A, 578, A6
*/
/*================================================================================*/
#include "pluto.h"
#include <math.h>

double ParkerVelocity(double, double, double, double, double);

void Init (double *v, double x1, double x2, double x3)
{

  double velocity;

#if COOLING != NO
  g_minCoolingTemp = g_inputParam[surface_temperature];
  g_maxCoolingRate = 0.1;
#endif

#if ROTATING_FRAME == YES
  double period = (g_inputParam[Omega]*86400/UNIT_TIME);
  g_OmegaZ = 2.0*CONST_PI/period;
#endif

  g_smallPressure = 1.0e-12/UNIT_PRESSURE; // Small value for pressure fix.
  g_gamma = 5.0/3.0;

  double mmw = MeanMolecularWeight(v);

  /* --- Stellar properties --- */

  double mass = g_inputParam[M_star]*CONST_Msun/UNIT_MASS;
  double radius = g_inputParam[R_star]*CONST_Rsun/UNIT_LENGTH;
  double magnetic_field_polar = g_inputParam[B_star]/UNIT_B;
  double surface_density = g_inputParam[star_density]/UNIT_DENSITY;

  /* --- Isothermal --- */

  double temperature = g_inputParam[T_star];
  double escape_velocity = sqrt(2.0*UNIT_G*mass);
  double sound_speed = sqrt(CONST_kB*temperature/(mmw*CONST_amu))/UNIT_VELOCITY;
  double lambda = 0.5*pow(escape_velocity/sound_speed,2);
  double critical_radius = UNIT_G*mass/(2.0*sound_speed*sound_speed);

  velocity = ParkerVelocity(
    sound_speed, escape_velocity, x1, critical_radius, radius);

  double surface_pressure = sound_speed*sound_speed*surface_density;

  double pressure = surface_pressure*exp(lambda*(radius/x1 - 1.0) - 0.5*pow(velocity/sound_speed,2));
  double density = (surface_density/surface_pressure)*pressure;

  /* --- Adiabatic. --- */

  //double surface_pressure = surface_density*(g_gamma - 1.0);

  // Hydrostatic.
  //double density = surface_density*(g_gamma - 1.0)/g_gamma*UNIT_G*mass/radius*pow(radius/x1, g_gamma/(g_gamma - 1.0));
  //double pressure = surface_density*pow(radius/x1, 1.0/(g_gamma - 1.0));

  //double pressure = CONST_kB*temperature*density/(mmw*CONST_amu);

  v[RHO] = density;
  EXPAND(v[VX1] = velocity;,
         v[VX2] = 0.0;,
         v[VX3] = 0.0;);
  v[PRS] = pressure;

#if PHYSICS == MHD
  #if BACKGROUND_FIELD == NO
    EXPAND(v[BX1] = magnetic_field_polar*pow(x1,-3)*cos(x2);,
           v[BX2] = (magnetic_field_polar/2.)*pow(x1,-3)*sin(x2);,
           v[BX3] = 0.0;)
  #endif
  #if BACKGROUND_FIELD == YES
    EXPAND(v[BX1] = 0.0;,
           v[BX2] = 0.0;,
           v[BX3] = 0.0;)
  #endif
#endif

}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

void Analysis (const Data *d, Grid *grid){}

#if BACKGROUND_FIELD == YES
void BackgroundField (double x1, double x2, double x3, double *B0)
{
  double Bq = g_inputParam[B_star]/UNIT_B;
  EXPAND(B0[0] = Bq*pow(x1,-3)*cos(x2);,
         B0[1] = (Bq/2.0)*pow(x1,-3)*sin(x2);,
         B0[2] = 0.0;)
}
#endif

void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
{

  int i, j, k;
  double velocity;

  double vradial, vtheta, vphi;
  double ***vx1 = d->Vc[VX1];
  double ***vx2 = d->Vc[VX2];
  double ***vx3 = d->Vc[VX3];
  double ***rho = d->Vc[RHO];
  double ***prs = d->Vc[PRS];

#if PHYSICS == MHD
  double ***bx1 = d->Vc[BX1];
  double ***bx2 = d->Vc[BX2];
  double ***bx3 = d->Vc[BX3];
#endif
  double *x1 = grid->xgc[IDIR];
  double *x2 = grid->xgc[JDIR];
  double *x3 = grid->xgc[KDIR];

  x1 = grid->xgc[IDIR];
  x2 = grid->xgc[JDIR];
  x3 = grid->xgc[KDIR];

  // - quantity | value | units.

  // general quatities.
  double vector[1];
  double mmw = MeanMolecularWeight(vector);

  // Planet & Star properties.
  double T_parker = g_inputParam[T_star];
  double T_inflow = g_inputParam[surface_temperature];
  double mass = g_inputParam[M_star]*CONST_Msun/UNIT_MASS;
  double radius = g_inputParam[R_star]*CONST_Rsun/UNIT_LENGTH;
  double magnetic_field_polar = g_inputParam[B_star]/UNIT_B;

  //double M_dotcgs = g_inputParam[M_dot_astro]*CONST_Msun/tyear;
  //double mass_loss = M_dotcgs*UNIT_TIME/UNIT_MASS;

  double escape_velocity = sqrt(2.0*UNIT_G*mass/radius);

  double sound_speed_parker = sqrt(CONST_kB*T_parker/(mmw*CONST_amu))
                              /UNIT_VELOCITY;
  double sound_speed_inflow = sqrt(CONST_kB*T_inflow/(mmw*CONST_amu))
                              /UNIT_VELOCITY;

  double critical_radius = UNIT_G*mass/(2.0*sound_speed_parker*sound_speed_parker);

  double omega = g_inputParam[Omega]*2.67e-6*UNIT_TIME;

  //RH = 1.0e-15/UNIT_DENSITY;//mass_loss/(4.0*CONST_PI*parker[0]);
  //double velocity = ParkerVelocity(sound_speed_parker, escape_velocity,
  //  1.0, critical_radius, radius);

  //double surface_density = mass_loss/(4.0*CONST_PI*velocity);

  double surface_density = g_inputParam[star_density]/UNIT_DENSITY;

  //double density = surface_density*(g_gamma - 1.0)/g_gamma*UNIT_G*mass/radius*pow(radius/1.0, g_gamma/(g_gamma - 1.0));

  //double pressure = surface_density*pow(radius/1.0, 1.0/(g_gamma - 1.0));

  velocity = ParkerVelocity(
    sound_speed_parker, escape_velocity, 1.0, critical_radius, radius);

  double surface_pressure = surface_density*T_inflow/(KELVIN*mmw);

  //double surface_pressure = sound_speed_inflow*sound_speed_inflow*surface_density;
  //double surface_pressure = surface_density*(g_gamma - 1.0);

  double pressure = surface_pressure;

  double density = surface_density;

  //parker[0] = sound_speed*pow(escape_velocity/(2.0*sound_speed), 2)
  //    *exp(-escape_velocity/(2.0*pow(sound_speed, 2)) + 0.666666)/UNIT_VELOCITY;
  //RHo = mass_loss/(4.0*CONST_PI*parker[0]);
  //Pres = sound_speed*sound_speed*RH;// should be *1.0/g_gamma;

  int ghost = 2;//(NX1_TOT - NX1)/2;

  if(side == X1_BEG){
    if(box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){

        /*
        if (vx1[k][j][ghost] > sound_speed){
          EXPAND(vradial = sound_speed;,
                 vtheta = 0.0;,
                 vphi = 0.0;)
        } else if (vx1[k][j][ghost] < -sound_speed){
          EXPAND(vradial = -sound_speed;,
                 vtheta = 0.0;,
                 vphi = 0.0;)
        } else if (vx1[k][j][ghost] < fabs(sound_speed)) {
        */

        EXPAND(vx1[k][j][i] = velocity;,
               vx2[k][j][i] = 0.0;,
               vx3[k][j][i] = 0.0;);


        rho[k][j][i] = density;
        //EXPAND(vx1[k][j][i] = 2.0*vx1[k][j][ghost] - vx1[k][j][ghost+1];,
        //       vx2[k][j][i] = 2.0*vx2[k][j][ghost] - vx2[k][j][ghost+1];,
        //       vx3[k][j][i] = 0.0;)
        prs[k][j][i] = pressure;

        #if PHYSICS == MHD
          #if BACKGROUND_FIELD == NO
            EXPAND(bx1[k][j][i] =
                     magnetic_field_polar*pow(x1[i], -3)*cos(x2[j]);,
                   bx2[k][j][i] =
                     (magnetic_field_polar/2.0)*pow(x1[i], -3)*sin(x2[j]);,
                   bx3[k][j][i] = 0.0;)
          #endif
          #if BACKGROUND_FIELD == YES
            EXPAND(bx1[k][j][i] = 0.0;,
                   bx2[k][j][i] = 0.0;,
                   bx3[k][j][i] = 0.0;)
          #endif
        #endif
      }
    }
  }
}


#if BODY_FORCE != NO
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
{
  double mass = g_inputParam[M_star]*CONST_Msun/UNIT_MASS;
  double gravity = -UNIT_G*mass/x1/x1;
  if (x1 > 1.0){
    g[IDIR] = gravity;
    g[JDIR] = 0.0;
    g[KDIR] = 0.0;
  }
}
#endif


#if USER_SOURCE == YES
void UserSource(double *v, double *S, double x1, double x2, double x3)
{
    double radius = 1.0;
    double Hs = g_inputParam[scale_hight]*CONST_Rsun/UNIT_LENGTH;

    //double ph_flux = g_inputParam[energy_flux]/(UNIT_MASS/pow(UNIT_TIME, 3));  // [erg cm^-2 s^1]
    double H_0 = g_inputParam[energy_flux]/(UNIT_MASS/(UNIT_LENGTH*pow(UNIT_TIME, 3)));

    //double Qh = ph_flux/Hs*pow(radius/x1, 2)*exp(-(x1 - radius)/Hs);
    double Qh = H_0*exp(-(x1 - radius)/Hs);

    /*
    double radius = 1.0;
    double unit_energy = (UNIT_MASS*pow(UNIT_LENGTH, 2)/pow(UNIT_VELOCITY, 2));
    double h0 = 4.9128e-7/(unit_energy/(pow(UNIT_LENGTH, 3)*UNIT_TIME);
    double l0 = 0.7;
    double Qh = h0*exp(-(x1 - radius)/l0);
    double Alf_pressure = 8.36e-2/UNIT_PRESSURE;

    if (x1 < radius) {

      double v2 = EXPAND(v[VX1]*v[MX1], + v[MX2]*v[MX2], + v[MX3]*v[MX3]);

      double m2 = EXPAND(v[MX1]*v[MX1], + v[MX2]*v[MX2], + v[MX3]*v[MX3]);

      double kinE = 0.5*m2/v[RHO];

      double pressure = (g_gamma - 1.0)*(v[ENG] - kinE);

      pressure += Alf_pressure;

      double energy = 0.5*v[RHO]*v2 + pressure*(1.0/(g_gamma - 1.0));

      // This is wrong!!!!!!!!

    }
    */

    S[RHO] = 0.0;
    EXPAND(S[MX1] = 0.0;,
           S[MX2] = 0.0;,
           S[MX3] = 0.0;)
#if PHYSICS == MHD
    EXPAND(S[BX1] = 0.0;,
           S[BX2] = 0.0;,
           S[BX3] = 0.0;)
#endif
    S[ENG] = Qh;

}
#endif


/*
  In order to obtain the velocity as a function of radial
  distance from either the star or planet it is nessacery
  to solve parkers equation and find its roots. There are
  two roots one for the subsonic and one for the supersonic
  regions. the transition between the two rigines is given
  by the critical point (see above expressions). The
  roots are found using the Newton raphson tschnique.
*/
double ParkerVelocity(double sound_speed, double escape_velocity, double x1, double critical_radius, double radius)
{
  int iteration = 0;
  double psi, step, dfn, fn;
  double lambda = 0.5*pow(escape_velocity/sound_speed,2);
  double eta = x1/radius;
  double b = -3.0 - 4.0*log(lambda/2.0) + 4.0*log(eta) + 2.0*(lambda/eta);

  if (x1 <= critical_radius) {
    psi = 2.e-8;
  } else {
    psi = 2.5;
  }
  do {
    fn = -psi + log(psi) + b;
    dfn = -1.0 + (1.0/psi);
    step = fn/dfn;
    psi = psi - step;
    iteration++;
  } while (fabs(step) > 1.0e-8 && iteration < 100000);
  if (isnan(psi)){
    psi = 0.0;
  }

  //parker[0] = sound_speed*sqrt(psi);
  //parker[1] = Pres*exp(lambda*(radius/x1-1.0)-0.5*pow(parker[0]/sound_speed,2));
  //parker[2] = (RHo/Pres)*parker[1];

  return sound_speed*sqrt(psi);

}
