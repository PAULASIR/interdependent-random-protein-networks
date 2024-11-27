#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<limits.h>


void main()
 {
   float k, i, sf, N, a;
   FILE *fp; 
   fp= fopen("sp2.dat","w");
   
   N = 200;
   a =2;
   
   for(k = 10; k <= 40000; k+=10){
      sf = N*8/ (float) k;
      fprintf(fp, "%f  %f\n", k, sf);
     }
   } 


