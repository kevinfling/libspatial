#include <libspatial/packed_bvh.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static uint32_t rng=7;
static uint32_t rnd(void){rng^=rng<<13;rng^=rng>>17;rng^=rng<<5;return rng;}
int main(void){
 for(uint32_t dims=3;dims<=256;dims*=2){
  uint32_t n=1024;float *points=malloc((size_t)n*dims*sizeof(float));uint64_t *keys=malloc(n*sizeof(*keys));float *q=malloc(dims*sizeof(float));
  for(size_t i=0;i<(size_t)n*dims;++i)points[i]=((int)(rnd()%20001)-10000)/10000.f;
  for(uint32_t i=0;i<n;++i)keys[i]=((uint64_t)rnd()<<32)|rnd();
  spatial_packed_bvh b={0};
  if(spatial_packed_bvh_build(&b,points,keys,n,dims,32))return 1;
  spatial_packed_bvh_queue *heap=malloc(b.node_count*sizeof(*heap));
  spatial_packed_bvh_hit got[10],expected[10];
  for(int qi=0;qi<100;++qi){
   for(uint32_t d=0;d<dims;++d)q[d]=((int)(rnd()%20001)-10000)/10000.f;
   for(int j=0;j<10;++j){expected[j].id=UINT32_MAX;expected[j].distance2=INFINITY;}
   for(uint32_t i=0;i<n;++i){spatial_packed_bvh_hit h={i,spatial_packed_dist2(q,points+(size_t)i*dims,dims)};int pos=0;
    while(pos<10&&!spatial_packed_hit_better(h,expected[pos]))++pos;
    if(pos<10){for(int j=9;j>pos;--j)expected[j]=expected[j-1];expected[pos]=h;}}
   size_t k=spatial_packed_bvh_knn(&b,q,10,got,(spatial_packed_bvh_scratch){heap,b.node_count});
   if(k!=10)return 2;
   for(int j=0;j<10;++j)if(got[j].id!=expected[j].id || fabs(got[j].distance2-expected[j].distance2)>1e-12){fprintf(stderr,"fail dims %u query %d pos %d\n",dims,qi,j);return 3;}
  }
  free(heap);spatial_packed_bvh_destroy(&b);free(points);free(keys);free(q);
  printf("packed bvh dims=%u pass\n",dims);
 }
 return 0;
}
