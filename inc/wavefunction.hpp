#ifndef __WAVEFUNCTION_H
#define __WAVEFUNCTION_H

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "atom.hpp"


namespace Wavefunction {

	double get_state_energy(Atom &myAtom);
	double core_potential(Atom &myAtom, double r);
	double effective_charge(Atom &myAtom, double r);
	double model_potential(Atom &myAtom, double r);
	double spin_orbit_interaction(Atom &myAtom, double r);
	void numerov_integration(Atom &myAtom, double r_inner, double r_outer, double step);
	double g(Atom &myAtom, double x);
	double get_defect(Atom &myAtom);

};



#endif // __WAVEFUNCTION_H	
