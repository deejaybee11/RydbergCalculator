#ifndef __ATOM_H
#define __ATOM_H

#include <stdlib.h>
#include <stdio.h>

class Atom {

public:

		Atom(int n, int l, double j, double s, int Z, int M);
		~Atom();

		int n;
		int l;
		double j;
		int Z;
        int M;
		double s;
		double* radialWavefunction = 0;
		double* r_vec = 0;
        double wflength = 0;

};
			
#endif // __ATOM_H
