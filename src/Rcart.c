#include "cart.h"

#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Memory.h>
#include <R_ext/Utils.h>

#ifdef CHECK_INTERRUPTS
#define INTERRUPT R_CheckUserInterrupts();
#else
#define INTERRUPT 
#endif

/*
  r_pop - the grid/2D-array of population values
            either a numeric matrix or a list of numeric vectors.

  gridx, gridy - vectors of x and y locations on the grid at which we
            want the transformed results.

  r_dims - integer vector of length 2 giving the dimensions of r_pop
           (recall that the columns correspond to the values of y at
           different x's. So r_dims = (ncol, nrow).)
  
  blur - a non-negative value giving a blur to the difussion.


 This modifies the value of gridx and gridy in place.
*/

SEXP
R_makecartogram(SEXP r_pop, SEXP gridx, SEXP gridy, SEXP r_dims, SEXP blur)
{
    int *dims;
    double **pop;
    int i, n;
    SEXP ans;

    dims = INTEGER(r_dims);    

    /* Note that if we get an interrupt, this will not be freed! 
       We may change cart.c to use R_alloc().
     */
    cart_makews(dims[0], dims[1]);

    
    pop = (double **) R_alloc(sizeof(double *), dims[0]);
    if(TYPEOF(r_pop) == VECSXP) {
	/* This leads to non-contiguous values. Not certain if this is
	 * allowed in the cart.c code. */
	for(i = 0; i < dims[0]; i++) {
	    pop[i] = REAL(VECTOR_ELT(r_pop, i));
	}
    } else {
	/* Given a numeric matrix, so just pull out the pointers to
	    the beginning of each column. Since we are using the
            columns for the Y values, this gives us the form of a 2D
            array in C with pop[i] giving a vector/array of the y's. 
          */
	double *data = REAL(r_pop);
	for(i = 0; i < dims[0]; i++) {
	    pop[i] = data + i * dims[1];
	}
    }

    INTERRUPT
    cart_transform(pop, dims[0], dims[1]);

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, gridx = Rf_duplicate(gridx));
    SET_VECTOR_ELT(ans, 1, gridy = Rf_duplicate(gridy));

    INTERRUPT
    cart_makecart(REAL(gridx), REAL(gridy), Rf_length(gridx), dims[0], dims[1], REAL(blur)[0]);

    INTERRUPT
    cart_freews(dims[0], dims[1]); 

    UNPROTECT(1);
    return(ans);
}


#define IJ(arr, i, j) \
   arr[ dims[0] * (i) + (j) ]

#include <stdio.h>

SEXP
R_predict(SEXP obj, SEXP r_x, SEXP r_y, SEXP r_ans, SEXP r_dims)
{
    int n, i;
    double *x_ans_ptr, *y_ans_ptr;
    double *x, *y;
    int ix, iy;
    double dx, dy;
    int *dims;
    double *gridx, *gridy;

    n = Rf_length(r_x);
    x_ans_ptr = REAL(VECTOR_ELT(r_ans, 0));
    y_ans_ptr = REAL(VECTOR_ELT(r_ans, 1));
    dims = INTEGER(r_dims);

    gridx = REAL(VECTOR_ELT(obj, 0));
    gridy = REAL(VECTOR_ELT(obj, 1));

    x = REAL(r_x);
    y = REAL(r_y);

    for(i = 0; i < n; i++) {
	ix = (int) x[i] - 1;
	iy = (int) y[i] - 1;
	dx = x[i] - ix - 1.;
	dy = y[i] - iy - 1.;

	if(ix < 0 || ix >= dims[0] || iy < 0 || iy >= dims[1]) {
	    continue;
	}

	x_ans_ptr[i] = (1.-dx) * (1.-dy)*IJ(gridx, ix, iy) + dx*(1.-dy)*IJ(gridx, ix+1, iy)
             + (1.-dx)*dy*IJ(gridx, ix, iy+1) + dx*dy*IJ(gridx, ix+1, iy+1);
	y_ans_ptr[i] = (1. - dx)*(1. - dy)*IJ(gridy, ix, iy) + dx*(1. - dy)*IJ(gridy, ix+1, iy)
             + (1. - dx)*dy*IJ(gridy, ix, iy+1) + dx*dy*IJ(gridy, ix+1, iy+1);
    }

    return(R_NilValue);
}
