//  SFM  04/06/2013  Added math.h
#include <math.h>


//[144][91][12]
//91,144,12
//col,row,mon
//0.1*geoData%sstdata(col+1,row+1,mm+1)
void pix_screen_ocean_( float *tb10h, float *tb10v,
			float *tb19h, float *tb19v, float *tb21v,
			float *tb37h, float *tb37v,
			float *tb85h, float *lat, float *lon,
			short int *sstdata,
			int *mm, int *rainflag, float *sst, float *sfc_wind,
			float *tpw) 
{
  float        T19H, T19V, T22V, T37H, T37V, T85H, T10H, T10V;
  float        d85_37, lwp, lwp_max;
  float        WARM85_OC, T22_COLD, T37_COLD;
  float        AMB_HI, AMB_LO, AMB_INT;
  float        DT_OC, PCP_INT, sstv, hfreez;
  int          i, j, row, col;

  float cw[6]={1.70156,-0.942322,-2.65166,1.85971,-0.227124,136.170};
  float cs[6]={3.98973,-1.87483,-2.11608,1.31708,-0.101934,-214.284};
  float ct[6]={1.26687,-2.01279,2.06152,1.30843,-0.980663,-344.024};
  float sstnew;
  WARM85_OC = 262.;
  T22_COLD  = 192. + 0.; //     ! This may be a bug is GSCAT
  T37_COLD  = 185. + 1.7;
  AMB_HI    = 262. + 7.1;
  AMB_LO    = 255. + 6.9;
  AMB_INT   = 158. + 5.3;
  DT_OC     =   0. + 0.5;
  PCP_INT   =  38. - 0.1;



 
  *rainflag=1;
  T10H=*tb10h;
  T10V=*tb10v;
  T19H=*tb19h;
  T19V=*tb19v;
  T22V=*tb21v;
  T37H=*tb37h;
  T37V=*tb37v;
  T85H=*tb85h;
  d85_37 = T85H-T37H ;
  if ((T22V >= 285. ) || (T37V >= 285.))
    {
      *rainflag = 0;
      goto endloop;
    }
  lwp = 0.399635*log(285. - T22V) - 1.40692*log(285. - T37V) + 4.299;
  row = (int)(*lon + 180.0)/2.5;
  col = (int)(*lat+ 90.0)/2.0; 
  sstv = 0.1*sstdata[144*91*(*mm)+(row)*91+col];
  hfreez=(sstv-273)/6.5;
  lwp_max = 0.2*hfreez/4.;
  /*	    printf("%g %g %g %g %g\n",lwp, lwp_max, hfreez, T22V,
    T37V);
  */
  if ((lwp < lwp_max) && ( T19H < 200.)) 
    {
      *rainflag = 0;
      goto endloop;
    }
  if ((T85H > WARM85_OC) && (T22V < AMB_HI))
    *rainflag = 1;
  else if (T22V <= T22_COLD)
    *rainflag= 0;
  else if ((d85_37 > DT_OC) && (T37H <= T37_COLD))
    *rainflag= 0;
  else if (T22V > (PCP_INT + 0.88*T19V))
    *rainflag= 1;
  else if ((T22V < AMB_LO) &&
	   (T22V < (AMB_INT + 0.49*T85H)))
    *rainflag= 1;
 
 endloop:
 if(*rainflag==1)
    {
      *sfc_wind=-99;
      *tpw=-99;
    }
  else
    {
      /*  *sfc_wind=147.9+1.0969*T19V-0.4555*T22V-1.7600*T37V+
	  0.786*T37H;
	  *tpw=15.*(23.82-4.059*log(280-T22V)+0.02451*(log(280-T22V)-T37V));
	  printf("%g %g \n",*sfc_wind,*tpw);
      */
    }
 /* *sfc_wind=cw[0]*T10V+cw[1]*T10H+cw[2]*T19V+cw[3]*T19H+
    cw[4]*T22V+cw[5];
 */
 if(*sfc_wind<0.1) *sfc_wind=0.1;
 /*
   sstnew=cs[0]*T10V+cs[1]*T10H+cs[2]*T19V+cs[3]*T19H+
   cs[4]*T22V+cs[5];
   printf("%g %g \n",sstnew,*sst);
   *tpw=ct[0]*T10V+ct[1]*T10H+ct[2]*T19V+ct[3]*T19H+
   ct[4]*T22V+ct[5];
   if(*tpw<0.1) *tpw=0.1;
 */
}
