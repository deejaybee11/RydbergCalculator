#include "atom.hpp"

#include <stdio.h>
#include <stdlib.h>

#include "mkl.h"

Atom::Atom(int n, int l, double j, double s, int Z, int M) {
	
	this->n = n;
	this->l = l;
	this->j = j;
	this->s = s;
	this->Z = Z;
    this->M = M;
}


Atom::~Atom() {
    mkl_free(r_vec);
    mkl_free(radialWavefunction);
}
