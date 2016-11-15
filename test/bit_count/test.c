#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int pop(unsigned  int x);

int main(int argc, char* argv[]){
 unsigned long int i,j,max_i;
 unsigned int u_i,l_i,i_32;
 long int cnt,u_cnt,l_cnt;
 double org,ratio;

 i_32   = (unsigned long int)pow(2,32)-1;
 max_i  = (unsigned long int)pow(2,61);
 max_i  = max_i-1;
 //i      = (unsigned long int)pow(2,32);
 i      = 1+2+4+8+16;
 cnt = pop(i);
 printf("%lu %ld\n",i,cnt);
 cnt = pop(max_i-1);
 printf("%lu %ld\n",max_i,cnt);

// make 32bit
 l_i = max_i & i_32;
 u_i = max_i >> 32;
 u_cnt = pop(u_i);
 l_cnt = pop(l_i);
 printf("%lu: %u %u: %d %d\n",max_i,u_i,l_i,u_cnt,l_cnt);

 return 0;
}

// for 32 bit
int pop(unsigned int x){
  x = x - ((x>>1) & 0x55555555);
  x = (x & 0x33333333)+ ((x>>2)& 0x33333333);
  x = (x+(x>>4)) & 0x0F0F0F0F;
  x = x+ (x>>8);
  x = x+ (x>>16);
  return  x & 0x0000003F;
}

/*
long int pop(unsigned long int x){
  unsigned long int y;
  y = x * 0x0002000400080010ULL;
  y = y & 0x1111111111111111ULL;
  y = y * 0x1111111111111111ULL;
  y = y >> 60;
  return  y;
}
*/
