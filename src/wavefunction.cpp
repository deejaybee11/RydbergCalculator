#include "wavefunction.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>

#include "mkl.h"
#include "mathimf.h"
#include "atom.hpp"
#include "constants.hpp"

double Wavefunction::effective_charge(Atom &myAtom, double r) {
	
	double alpha1 = disp[myAtom.l][0];
	double alpha2 = disp[myAtom.l][1];
	double alpha3 = disp[myAtom.l][2];
	double alpha4 = disp[myAtom.l][3];
	double Qe = 0;

	Qe = (1 + (myAtom.Z-1) * exp(-alpha1 * r) - r * (alpha3 + alpha4 * r) * exp(-alpha2 * r));
	return Qe;
}


//core potential includes the core polarizability term
double Wavefunction::core_potential(Atom &myAtom, double r) {

	double rc = disp[myAtom.l][4];	
	double Vc = -effective_charge(myAtom, r) / r - c_alpha_c/(2*pow(r, 4.0)) * (1 - exp(-pow(r/rc, 6.0)));
	return Vc;
}

double Wavefunction::spin_orbit_interaction(Atom &myAtom, double r) {
	
	double Vso = 0.5 * pow(c_alpha, 2.0) / (2*pow(r, 3.0)) * (myAtom.j*(myAtom.j + 1) - myAtom.l*(myAtom.l + 1) - myAtom.s*(myAtom.s + 1));
  	return Vso;	
}


double Wavefunction::model_potential(Atom &myAtom, double r) {

	double V;
	if (myAtom.l < 4) {
		V = spin_orbit_interaction(myAtom, r) + core_potential(myAtom, r);
	}
	else {
		V = (-1.0/r) + spin_orbit_interaction(myAtom, r); 
	}

	return V;
}

double Wavefunction::g(Atom &myAtom, double x) {

	double r = x * x;
	double mass = myAtom.M*c_amu;
	double mu = (mass - c_me) / mass;
	double g = -3/(4.0 * r) + 4*r*(2* mu * (get_state_energy(myAtom)/27.211 - model_potential(myAtom, r)) - myAtom.l*(myAtom.l + 1)/pow(r, 2.0));
	return g;
}


void Wavefunction::numerov_integration(Atom &myAtom, double r_inner, double r_outer, double step) {
	
	double *res = 0;
	int npts = (int)((sqrt(r_outer)-sqrt(r_inner))/step);
    myAtom.wflength = npts;
	printf("Time to start the integration with npts=%d\n", npts);
	printf("The state energy is E=%f\n", get_state_energy(myAtom));
	printf("The quantum defect is delta=%f\n", get_defect(myAtom));
	res = (double*)mkl_malloc(2*npts*sizeof(double), 64);
	myAtom.radialWavefunction = (double*)mkl_malloc(npts*sizeof(double), 64);
	myAtom.r_vec = (double*)mkl_malloc(npts*sizeof(double), 64);

    r_inner = std::max(4.0*step, r_inner);
	
	//perform a couple of independent calculations to setup a couple of points to allow back integration
	int index = npts-1;
	double x = sqrt(r_inner) + step*((double)npts - 1);
	double step2 = step*step;
	double init1 = 0.01;
	double init2 = 0.01;

	res[index] = (2 * (1 - 5.0/12.0*step2*g(myAtom, x))*init1 - (1 + 1/12.0 * step2 * g(myAtom, x+step))*init2)/(1+1/12.0*step2*g(myAtom, x-step));
    std::cout << "The g function thing = " << g(myAtom, x) << std::endl;
	res[index+npts] = x;
	x = x - step;
	index = index - 1;

	res[index] = (2 * (1 - 5.0/12.0*step2*g(myAtom, x))*res[index+1] - (1 + 1/12.0 * step2 * g(myAtom, x+step))*init1)/(1+1/12.0*step2*g(myAtom, x-step));
	res[index+npts] = x;

	double maxVal = 0;
	int check = 0;
	int fromLastMax = 0;

	while (index > check) {
		
		index = index - 1;
		x = x-step;
		res[index] = (2 * (1 - 5.0/12.0*step2*g(myAtom, x))*res[index+1] - (1 + 1/12.0 * step2 * g(myAtom, x+step))*res[index+2])/(1+1/12.0*step2*g(myAtom, x-step));
		res[index + npts] = x;
		if (fabs(res[index] * sqrt(x)) > maxVal) {
				maxVal = fabs(res[index]*sqrt(x));
		}
		else {
			fromLastMax += 1;
			if (fromLastMax > 50) {
				check = index;
			}
		}
	}

    int divPoint = 0;
    while ((index > 0) && (divPoint == 0)) {
        index = index - 1;
        x = x - step;
        res[index] = (2 * (1 - 5.0/12.0 * step2 * g(myAtom, x)) * res[index+1] - (1 + 1/12.0 * step2 * g(myAtom, x+step)) * res[index + 2]) / (1 + 1/12.0 * step2 * g(myAtom, x-step));
        res[index + npts] = x;

        if ((divPoint == 0) && (fabs(res[index] * sqrt(x)) > maxVal)) {
            divPoint = index;
            while ((fabs(res[divPoint]) > fabs(res[divPoint + 1])) && (divPoint < check)) {
                divPoint += 1;
            }
            if (divPoint > check) {
                printf("WHOLE THING BROKE\n");
            }
        }
    }
	
    double sum = 0;

    double dr = res[10+npts]*res[10+npts] - res[9+npts]*res[9+npts];
	for (int i=0; i<npts; i++) {
		myAtom.radialWavefunction[i] = res[i]*sqrt(res[i+npts]);
		myAtom.r_vec[i] = res[i+npts]*res[i+npts];
        sum += pow(myAtom.radialWavefunction[i], 2.0)*dr;
	}
    std::cout << "The sum is = " << sum << std::endl;
    //Normalize the wavefunction
    for (int i=0; i < npts; i++) {
        myAtom.radialWavefunction[i] = myAtom.radialWavefunction[i] / sqrt(sum);
    }

    mkl_free(res);
}

double Wavefunction::get_defect(Atom &myAtom) {

	switch(myAtom.l) {
		case 0: return coeffS12[0] + coeffS12[1]/pow(myAtom.n - coeffS12[0], 2.0);
		case 1: return myAtom.j==0.5 ? coeffP12[0] + coeffP12[1]/pow(myAtom.n - coeffP12[0], 2.0) : coeffP32[0] + coeffP32[1]/pow(myAtom.n - coeffP32[0], 2.0);
		case 2: return myAtom.j==1.5 ? coeffD32[0] + coeffD32[1]/pow(myAtom.n - coeffD32[0], 2.0) : coeffD52[0] + coeffD52[1]/pow(myAtom.n - coeffD52[0], 2.0);
		case 3: return myAtom.j==2.5 ? coeffF52[0] + coeffF52[1]/pow(myAtom.n - coeffF52[0], 2.0) : coeffF72[0] + coeffF72[1]/pow(myAtom.n - coeffF72[0], 2.0);
		case 4: return myAtom.j==3.5 ? coeffG72[0] + coeffG72[1]/pow(myAtom.n - coeffG72[0], 2.0) : coeffG92[0] + coeffG92[1]/pow(myAtom.n - coeffG92[0], 2.0);
		case 5: return coeffH92[0] + coeffH92[1]/pow(myAtom.n - coeffH92[0], 2.0);
	}
	return 0;
}
double Wavefunction::get_state_energy(Atom &myAtom) {

	double mass = myAtom.M*c_amu;
	double mu = (mass - c_me)/mass;
    double defect = get_defect(myAtom); 
    double scaledRydberg = R_c * mu * 1.239841984e-6;
	double stateEnergy = -scaledRydberg / (pow(myAtom.n - defect, 2.0));
	return stateEnergy;
}
