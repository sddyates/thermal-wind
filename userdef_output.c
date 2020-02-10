#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;

  double a11, a12, a13, a21, a22, a23, a31, a32, a33;

  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  double ***rho = d->Vc[RHO];
  double ***prs = d->Vc[PRS];

  double ***vx1 = d->Vc[VX1];
  double ***vx2 = d->Vc[VX2];
  double ***vx3 = d->Vc[VX3];

  #if PHYSICS == MHD
    double ***bx1 = d->Vc[BX1];
    double ***bx2 = d->Vc[BX2];
    double ***bx3 = d->Vc[BX3];
  #endif

  double ***Bx1_tot = GetUserVar("Bx1_tot");
  double ***Bx2_tot = GetUserVar("Bx2_tot");
  double ***Bx3_tot = GetUserVar("Bx3_tot");

  double ***vx1_tot = GetUserVar("vx1_tot");
  double ***vx2_tot = GetUserVar("vx2_tot");
  double ***vx3_tot = GetUserVar("vx3_tot");

  double ***Qh = GetUserVar("Qh");
  double ***Temp = GetUserVar("Temp");

#if USER_SOURCE == YES
  double ph_flux = g_inputParam[energy_flux]*pow(UNIT_TIME, 3)/UNIT_MASS;
  double Hs = g_inputParam[scale_hight];
#endif

  DOM_LOOP(k,j,i){

#if USER_SOURCE == YES
    double Qh_source = ph_flux/Hs*pow(1.0/x1[i], 2)*exp(-(x1[i] - 1.0)/Hs);
    Qh[k][j][i] = Qh_source;
#else
    Qh[k][j][i] = 0.0;
#endif

    double v[1];
    double mmw = MeanMolecularWeight(v);
    Temp[k][j][i] = prs[k][j][i]/rho[k][j][i]*KELVIN*mmw;
    //Temp[k][j][i] = prs[k][j][i]*UNIT_VELOCITY*UNIT_VELOCITY*mmw*CONST_mH/(CONST_kB*rho[k][j][i]);

#if PHYSICS == MHD
    double B1 = bx1[k][j][i];
    double B2 = bx2[k][j][i];
    double B3 = bx3[k][j][i];

#if BACKGROUND_FIELD == YES
    double Bq = g_inputParam[B_star]/UNIT_B;
    B1 += Bq*pow(x1[i],-3)*cos(x2[j]);
    B2 += (Bq/2.0)*pow(x1[i],-3)*sin(x2[j]);
    B3 += 0.0;
#endif


    a11 = sin(x2[j])*cos(x3[k]); a12 = cos(x2[j])*cos(x3[k]); a13 = -sin(x3[k]);
    a21 = sin(x2[j])*sin(x3[k]); a22 = cos(x2[j])*sin(x3[k]); a23 = cos(x3[k]);
    a31 = cos(x2[j]);            a32 = -sin(x2[j]);           a33 = 0.0;

    double B1C = B1*a11 + B2*a12 + B3*a13;
    double B2C = B1*a21 + B2*a22 + B3*a23;
    double B3C = B1*a31 + B2*a32 + B3*a33;


    EXPAND(Bx1_tot[k][j][i] = B1C*UNIT_B;,
           Bx2_tot[k][j][i] = B2C*UNIT_B;,
           Bx3_tot[k][j][i] = B3C*UNIT_B;)
#endif

#if PHYSICS == HD
    EXPAND(Bx1_tot[k][j][i] = 0.0;,
           Bx2_tot[k][j][i] = 0.0;,
           Bx3_tot[k][j][i] = 0.0;)
#endif


    double v1C = vx1[k][j][i]*a11 + vx2[k][j][i]*a12 + vx3[k][j][i]*a13;
    double v2C = vx1[k][j][i]*a21 + vx2[k][j][i]*a22 + vx3[k][j][i]*a23;
    double v3C = vx1[k][j][i]*a31 + vx2[k][j][i]*a32 + vx3[k][j][i]*a33;

    EXPAND(vx1_tot[k][j][i] = v1C*UNIT_VELOCITY;,
           vx2_tot[k][j][i] = v2C*UNIT_VELOCITY;,
           vx3_tot[k][j][i] = v3C*UNIT_VELOCITY;)

  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/*
 *
 *
 *************************************************************** */
{
  Image *image;

#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}
