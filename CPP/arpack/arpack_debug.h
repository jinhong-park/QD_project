#ifndef _ARPACK_DEBUG_

#define _ARPACK_DEBUG_

/*
  This function sets all ARPACK FORTRAN debug variables to zero.
*/
void TraceOff();

/*
  This function sets the values of ARPACK FORTRAN debug variables
*/
void TraceOn(char which, int aup2, int aitr, int gets, int eupd);

#endif

