#include "arpack_debug.h"

#define F77NAME(x) x ## _

typedef int integer;

extern "C"
{

struct {
    integer logfil, ndigit, mgetv0;
    integer msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
    integer mnaupd, mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd;
    integer mcaupd, mcaup2, mcaitr, mceigt, mcapps, mcgets, mceupd;
  } F77NAME(debug);
}

void TraceOff()
/*
  This function sets all ARPACK FORTRAN debug variables to zero.
*/

{

  F77NAME(debug).logfil = 6;
  F77NAME(debug).ndigit = 0;
  F77NAME(debug).mgetv0 = 0;
  F77NAME(debug).msaupd = 0;
  F77NAME(debug).msaup2 = 0;
  F77NAME(debug).msaitr = 0;
  F77NAME(debug).mseigt = 0;
  F77NAME(debug).msapps = 0;
  F77NAME(debug).msgets = 0;
  F77NAME(debug).mseupd = 0;
  F77NAME(debug).mnaupd = 0;
  F77NAME(debug).mnaup2 = 0;
  F77NAME(debug).mnaitr = 0;
  F77NAME(debug).mneigt = 0;
  F77NAME(debug).mnapps = 0;
  F77NAME(debug).mngets = 0;
  F77NAME(debug).mneupd = 0;
  F77NAME(debug).mcaupd = 0;
  F77NAME(debug).mcaup2 = 0;
  F77NAME(debug).mcaitr = 0;
  F77NAME(debug).mceigt = 0;
  F77NAME(debug).mcapps = 0;
  F77NAME(debug).mcgets = 0;
  F77NAME(debug).mceupd = 0;

}


void TraceOn(char which, int aup2, int aitr, int gets, int eupd)
{
if (which=='s') {
  F77NAME(debug).logfil = 6;
  F77NAME(debug).ndigit = -3;
  F77NAME(debug).msaupd = 1;
  F77NAME(debug).msaup2 = aup2;
  F77NAME(debug).msaitr = aitr;
  F77NAME(debug).msgets = gets;
  F77NAME(debug).mseupd = eupd;
} 
else if (which=='n') {
  F77NAME(debug).logfil = 6;
  F77NAME(debug).ndigit = -3;
  F77NAME(debug).mnaupd = 1;
  F77NAME(debug).mnaup2 = aup2;
  F77NAME(debug).mnaitr = aitr;
  F77NAME(debug).mngets = gets;
  F77NAME(debug).mneupd = eupd;

}
else if (which=='z') {
  F77NAME(debug).logfil = 6;
  F77NAME(debug).ndigit = -3;
  F77NAME(debug).mcaupd = 1;
  F77NAME(debug).mcaup2 = aup2;
  F77NAME(debug).mcaitr = aitr;
  F77NAME(debug).mcgets = gets;
  F77NAME(debug).mceupd = eupd;

}
}
