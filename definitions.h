#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     VECTOR
#define  FORCED_TURB                    NO
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            13

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   DIV_CLEANING
#define  BACKGROUND_FIELD               YES
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             SUPER_TIME_STEPPING
#define  VISCOSITY                      SUPER_TIME_STEPPING
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  M_star                         0
#define  T_star                         1
#define  R_star                         2
#define  B_star                         3
#define  Omega                          4
#define  star_density                   5
#define  energy_flux                    6
#define  scale_hight                    7
#define  surface_temperature            8
#define  kappa                          9
#define  transition_temperature         10
#define  resistivity_eta                11
#define  viscosity_nu                   12

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  INTERNAL_BOUNDARY              NO
#define  LIMITER                        MINMOD_LIM
#define  UNIT_DENSITY                   (1.0e-15)
#define  UNIT_LENGTH                    CONST_Rsun
#define  UNIT_VELOCITY                  (1.0e5)
#define  UNIT_PRESSURE                  (UNIT_VELOCITY*UNIT_VELOCITY*UNIT_DENSITY)
#define  UNIT_MASS                      (UNIT_DENSITY*pow(UNIT_LENGTH,3))
#define  UNIT_TIME                      (UNIT_LENGTH/UNIT_VELOCITY)
#define  UNIT_MASS                      (UNIT_DENSITY*pow(UNIT_LENGTH,3))
#define  UNIT_G                         (CONST_G*((1/pow(UNIT_LENGTH,3))*UNIT_MASS*pow(UNIT_TIME,2)))
#define  UNIT_kB                        ((CONST_kB*pow(UNIT_TIME,2))/(UNIT_DENSITY*pow(UNIT_LENGTH,5)))
#define  UNIT_B                         (sqrt(4.0*CONST_PI*UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY))
#define  tyear                          3.15569e+7
#define  tday                           8.64e+4
#define  VTK_TIME_INFO                  YES
#define  VTK_VECTOR_DUMP                NO
#define  GLM_EXTENDED                   NO
#define  CHOMBO_REF_VAR                 RHO
#define  PPM_ORDER                      4
#define  SHOW_TIME_STEPS                NO
#define  CAK                            NO
#define  USER_SOURCE                    YES

/* [End] user-defined constants (do not change this line) */
