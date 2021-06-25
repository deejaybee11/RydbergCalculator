#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include <stdlib.h>

#pragma once
//Rubidium quantum defect coefficients
const double coeffS12[2] = {3.13118078, 0.1787};
const double coeffP12[2] = {2.6548849, 0.29006};
const double coeffP32[2] = {2.6416737, 0.29507};
const double coeffD32[2] = {1.348094811, -0.60544};
const double coeffD52[2] = {1.346462211, -0.59404};
const double coeffF52[2] = {0.0165192, -0.085};
const double coeffF72[2] = {0.0165437, -0.086};
const double coeffG72[2] = {0.004, 0.0};
const double coeffG92[2] = {0.004, 0.0};
const double coeffH92[2] = {0.000, 0.0};
const double R_c = 10973731.568160;
//Rubidium dispersion coefficients for the screened Coulomb potential
//index [l,i] where l is the current angular momentum number l and i picks from alpha0,1,2,3 and rc
//from PRA 49, 982
const double disp[4][5] = {{3.69628474,1.64915255,-9.86069196,0.19579987,1.66242117},{4.44088978,1.92828831,-16.79597770,-0.81633314,1.50195124},{3.78717363,1.5727864,-11.65588970,0.52942835,4.86851938},{2.39848933,1.76810544,-12.07106780,0.77256589,4.79831327}};
const double c_alpha_c = 9.0760;//318.1; //in a.u.
//Physical constants
const double c_c = 2.99792458e8;
const double c_h = 6.626e-34;
const double c_hbar = 1.054e-34;
const double c_pi = 3.14159265359;
const double c_a0 = 1.6605390666e-27;
const double c_Ry = 10973731.56816;
const double c_me = 9.1093837015e-31;
const double c_amu = 1.67262192369e-27;
const double c_e = 1.602176634e-19;
const double c_k = 1.380649e-23;
const double c_e0 = 8.8541878128e-12;
const double c_alpha = 0.0072973525693;





#endif // __CONSTANTS_H
