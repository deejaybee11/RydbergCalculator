#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "mathimf.h"
#include "atom.hpp"
#include "wavefunction.hpp"
#include "constants.hpp"
#include "additional_functions.hpp"

int main() {


	double	step = 0.001;
	Atom atom1(15,0,0.5,0.5,37,87);
	double r_inner = pow(c_alpha_c, 1.0/3.0);
	double r_outer = 2.0*atom1.n*(atom1.n+15);
	printf("Created and atom with n=%d, l=%d, j=%f, s=%f, Z=%d\n", atom1.n, atom1.l, atom1.j, atom1.s, atom1.Z);
	printf("Attempt to solve ground state wavefunction, expect it to break\n");
	Wavefunction::numerov_integration(atom1, r_inner, r_outer, step);
	printf("solved????\n");
	printf("1st point of the wavefunction = %.20f\n", atom1.radialWavefunction[0]);
    save_file(atom1.r_vec, atom1.radialWavefunction, atom1.wflength, "testwf.txt");
	return 0;
}
