#include <flann/flann.h>

#include <stdio.h>
#include <stdlib.h>

/*allocate(hFreqPRg(dPRData%n1c21,49,4))
  call flannint(geoloc,  prgeoloc, hFreqTbs, hFreqPRg, (iEnd+1-iStart)*81, 2, &
  49*dPRData%n1c21)
*/

void flannint_(float *dataset, float *cset, float *gmiT, float *prT, 
	       int *rows, int *cols, int *tcount)
{

  int nn, i, k;
  int* result;
  float* dists;
  struct FLANNParameters p;
  float speedup;
  flann_index_t index_id;
  printf("%i %i %i \n",*cols,*rows,*tcount);

  //  *cols = 2;

  nn = 4;
  result = (int*) malloc(*tcount*nn*sizeof(int));
  dists = (float*) malloc(*tcount*nn*sizeof(float));
    
  p = DEFAULT_FLANN_PARAMETERS;
  p.algorithm = FLANN_INDEX_KDTREE;
  p.trees = 8;
  p.log_level = FLANN_LOG_INFO;
  p.checks = 64;
  p.random_seed=1976;
  printf("Computing index.\n");
  index_id = flann_build_index(dataset, *rows, *cols, &speedup, &p);
  
  flann_find_nearest_neighbors_index(index_id, cset, 
				     *tcount, result, dists,
				     nn, &p);
  int itb=0, ii;
  float tb=0;
  for(i=0;i<*tcount;i++)
    {
      /*printf("%4i %4i %4i  \n",
	result[i*4],result[i*4+1],result[i*4+2]);
	printf("%6.2f %6.2f %6.2f  \n",
	dists[i*4],dists[i*4+1],dists[i*4+2]);
      */
      tb=0;
      float sum=0;
      for(ii=0;ii<4;ii++)
	{
	  if(dists[i*4+ii]<0.1 && gmiT[result[i*4+ii]]>50)
	    {
	      tb+=1./(dists[i*4+ii]+1e-3)*gmiT[result[i*4+ii]];
	      sum+=1./(dists[i*4+ii]+1e-3);
	    }
	}
      if(sum>0.0001)
	tb/=sum;
      else
	tb=-99;
      prT[itb]=tb;
      //printf("%i %i %g \n",itb,i,tb);
      itb++;
    }
  flann_free_index(index_id, &p);
  //free(dataset);
  // free(testset);
  free(result);
  free(dists);
  
  //return 0;

}
