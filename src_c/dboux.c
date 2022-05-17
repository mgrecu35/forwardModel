//#include "builtin.hpp"
//#include "math.hpp"
//#include "emission.hpp"

//namespace __emission__ {

//using __math__::exp;
#include "math.h"
int __;
#define FAST_FOR(i, l, u) \
  for(__=l;__<u;__++) { i=__; //printf("%i \n",i);

#define END_FOR }

float dboux_(float *k, float *T, float *_, float *_0_, float *_1_, float *_2_, int* _3_) 
{
  float dz=*_;
  float eps=*_0_;
  float umu=*_1_;
  float Ts=*_2_;
  int n=*_3_;
  float Tb, att, dz2;
  int i;

  Tb = 0.0;
  att = 0.0;
  dz2 = ((dz/2.0)/umu);
  //printf("%i \n",n);
  // for(i=0;i<n;i++)
  FAST_FOR(i, 0, n)
    att = (att+k[n-1-i]*dz2);
    Tb = (Tb+k[n-1-i]*dz2*2*exp(-att)*T[n-1-i]);
    att = (att+k[n-1-i]*dz2);
  END_FOR
  

  Tb = Tb+eps*Ts*exp(-att);
  
  for(i=0;i<n;i++)
    {
      att = (att+k[i]*dz2);
      Tb = Tb+(1-eps)*k[i]*dz2*2*exp(-att)*T[i];
      att = att+k[i]*dz2;
    }

  return Tb;
}


void dboux2_(float *tb, float *tb0, float *emtb,
	     int *ifreq, int *ipolG, float *f)
{
  int i;
  for(i=0;i<9;i++)
    {
      f[i]=0.;
      if(tb[i]>tb0[i] && emtb[i]>tb0[i])
	{
	  if(tb[i]>emtb[i])
	    {
	      tb0[i]=emtb[i];
	      f[i]=1;
	    }
	  else
	    {
	      f[i]=(tb[i]-tb0[i])/(emtb[i]-tb0[i]);
	      tb0[i]=tb[i];

	    }
	}
    }
}
