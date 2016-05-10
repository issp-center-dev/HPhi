#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

unsigned long int snoob(unsigned long int x);

int main(int argc, char* argv[]){
 unsigned long int i,j,max_i;
 int cnt;
 double org,ratio;

 max_i  = (unsigned long int)pow(2,48);
 //i      = (unsigned long int)pow(2,32);
 i      = 1+2+4+8+16;
 org    = pow(2,64);
 cnt    = 0;
 printf(" %d %lu \n",cnt,i);
 while(j < max_i){
   j     = snoob(i); 
   ratio = 1.0*j/org;
   printf(" %d %lu : ratio=%lf\n",cnt,j,ratio);
   i     = j;
   cnt  +=1;
 }

 return 0;
}

unsigned long int snoob(unsigned long int x){
  unsigned long int smallest, ripple, ones;
  smallest = x &(-x);
  ripple   = x+ smallest;
  ones     = x ^ ripple;
  ones     = (ones>>2)/smallest;
  return   ripple|ones;
}
