/*random network wirth active and inactive nodes*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<limits.h>


int nn=10,n=3;
double omega,eps,eps1,lambda;
double y[110][220], yy[110][220];
double alpha, beta;
double meanx, meany;
double mn;


void RK4(int,int,double,double,double[110][220], double[110][220],
                   void (*DGL)(double,double[110][220],double[110][220]),
                   void (*DGL1)(double,double[110][220],double[110][220]));
void DGL(double, double[110][220],double[110][220]);
void DGL1(double, double[110][220],double[110][220]);



/****************************function to generate random coupling matrices***********/

void assign(int row, int col, int arr[row][col])
{
    srand(time(0));

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++) {
            arr[i][j] = rand()%2;
        }
    }
}
/****************************************************/

void main()
{
//nn= no 0f oscillators, n=dimension of the model//
double t,h,pi;
double z[105000];
double amplitude,vmin,vmax,c_max,c_min;
double x_max,x_min;
double sum, amp;
double sum1, sum2;

int i,j,k;
srand(time(0));

FILE *fp1, *fp2, *fp3, *fp4;
fp1=fopen("1.dat","w");
fp2=fopen("2.dat","w");
fp3=fopen("3.dat","w");
fp4=fopen("4.dat","w");

/*initial conditions*/
for(j=1;j<=nn;j++)
 {
  y[j][1]=(float) rand()/(double)RAND_MAX*-2.0+1.0; 
  y[j][2]=(float) rand()/(double)RAND_MAX*-2.0+1.0;
  y[j][3]=0.05;
 }


/*initial conditions*/
for(j=1;j<=nn;j++)
 {
  yy[j][1]=(float) rand()/(double)RAND_MAX*-2.0+1.0; 
  yy[j][2]=(float) rand()/(double)RAND_MAX*-2.0+1.0;
  yy[j][3]=0.05;
 }


   beta = 0.06;
   omega = 2.0; 
   alpha = 0.06;

   scanf("%lf",&eps);
   scanf("%lf",&eps1);
   scanf("%lf",&lambda);
//***time step***//
 h=0.01; t=0.0;

x_max= -INT_MAX; x_min=INT_MAX;
c_max=c_min=0;
vmax= vmin=0;


//------------------------------------------time loop-----------------------------
 for(k=1;k<=60000;k++)
   {               
      t=h*(double)(k);
      RK4(nn,n,h,t,y,yy,DGL,DGL1); 

for(j=1;j<=nn;j++)
    {          
/*    leaving transients*/               
	if(k>=1000) 
	  {          
/*---------------------------(avg)amplitude calculation---------------------*/
       //  fprintf(fp2,"%d   %f %f \n", j, t, y[j][1]);
      sum1=0.0;
       fprintf(fp1,"%f \t",t); 
          for(i=1;i<=nn;i++)
             fprintf(fp1,"%f  %f\t",y[i][1], sum1+=y[i][1]);    
       fprintf(fp1," \n");
      
     
      sum2=0.0;
       fprintf(fp2,"%f \t",t); 
          for(i=1;i<=nn;i++)
             fprintf(fp2,"%f  %f\t",yy[i][1], sum2+=yy[i][1]);    
      fprintf(fp2," \n");
     
                    
       fprintf(fp3,"%f  %f  %f   %f\n", t, sum1/nn, sum2/nn, (fabs(sum1/nn) - fabs(sum2/nn)));                            
                     
         
         mn = sum1/nn;
                      
        z[k] = mn;
	    if( z[k] > x_max)
             x_max = z[k];
            if (z[k] < x_min)
              x_min = z[k];   
  
          }
	}    
}

   amp = fabs(x_max -x_min)/2;
   
   printf("%f \n", amp);


printf("process over!!\n");
}
/***********************RK4 SUBROUTINE*********************************/
void RK4(int nn,int n,double h,double t,double y[110][220], double yy[110][220],
	   void (*DGL)(double,double[110][220],double[110][220]),
	    void (*DGL1)(double,double[110][220],double[110][220]))
   
      {   
             int i,j,k;
          
          
	   double k1[110][220],k2[110][220],k3[110][220],k4[110][220];
	   double yaux[110][220];

	   double kk1[110][220],kk2[110][220],kk3[110][220],kk4[110][220];
	   double yaux1[110][220];
	  
	  
	   DGL(t,y,k1);
	   for(j=1;j<=nn;j++)
	   {
          for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k1[j][i]/2.0;
	   }
	   
	   
	   DGL1(t,yy,kk1);
	   for(j=1;j<=nn;j++)
	   {
           for(i=1;i<=n;i++)
	     yaux1[j][i]=yy[j][i]+h*kk1[j][i]/2.0;
	   }
	   
	   
	   DGL(t+h/2.0,yaux,k2);
	   for(j=1;j<=nn;j++)
	   {
          for(i=1;i<=n;i++)
	    yaux[j][i]=y[j][i]+h*k2[j][i]/2.0;
	   }
	   
	  DGL1(t+h/2.0,yaux1,kk2);
	   for(j=1;j<=nn;j++)
	   {
          for(i=1;i<=n;i++)
	    yaux1[j][i]=yy[j][i]+h*kk2[j][i]/2.0;
	   }
	   
	   
	  DGL(t+h/2.0,yaux,kk3);
	   for(j=1;j<=nn;j++)
	    {
            for(i=1;i<=n;i++)
	      yaux[j][i]=y[j][i]+h*k3[j][i];
	    }   

	 
	  DGL1(t+h/2.0,yaux1,k3);
	   for(j=1;j<=nn;j++)
	   {
            for(i=1;i<=n;i++)
	      yaux1[j][i]=yy[j][i]+h*kk3[j][i];
	   }
	   
	   
	   DGL(t+h,yaux,k4);
	   for(j=1;j<=nn;j++)
	   {
             for(i=1;i<=n;i++)
	       y[j][i]=y[j][i]+h*((k1[j][i]+2*k2[j][i]+2*k3[j][i]+k4[j][i])/6.0);
	   }
	   
	    
	   DGL1(t+h,yaux1,kk4);
	   for(j=1;j<=nn;j++)
	    {
             for(i=1;i<=n;i++)
	       yy[j][i]=yy[j][i]+h*((kk1[j][i]+2*kk2[j][i]+2*kk3[j][i]+kk4[j][i])/6.0);
	    }	   
	   
}

//*********************FUNCTION SUBROUTINE********************************//
void DGL(double t,double y[110][220],double F[110][220])
 {

 int row = nn;
 int col = nn;
 int arr[row][col];
 double meanx,meany;
 double eta, a;
 
 int i,j;

 assign(row, col, arr);
 
  
  meanx = 0.0;
  for(j=1;j<=nn;j++)
      meanx += y[j][1];
      
      meanx = meanx /nn;
      
      
       
  meany = 0.0;
  for(j=1;j<=nn;j++)
      meany += yy[j][1];
      
      meany = meany /nn;
  
  
  for(j = 1; j <= nn; j++) {  
     double a = (float) rand()/(double)RAND_MAX*-0.1+0.05;    
       if(j >= nn/2)
          eta = 0;
           else
          eta =1;
                     
 for(i = 1; i <= nn; i++) {  
 F[j][1]=  y[j][2] + eps*arr[j][i]*(y[i][1] - y[j][1]) + eps1*(meany - meanx) + lambda*y[j][3] +eta*a;                                                                    
     F[j][2]=  -alpha*y[j][2] -omega*y[j][1] - beta*y[j][1]*y[j][1]*y[j][1]; 
     F[j][3] = -meanx - meany - y[j][3];
      }
   }
}

/******************************************************/
//*********************FUNCTION SUBROUTINE********************************//
void DGL1(double t,double yy[110][220],double F[110][220])
 {

 int row = nn;
 int col = nn;
 int arr[row][col];
 double meanx,meany;
 double a,eta;
 
 int i,j;

 assign(row, col, arr);
 
  
  meanx = 0.0;
  for(j=1;j<=nn;j++)
      meanx += y[j][1];
      
      meanx = meanx /nn;
      
      
       
  meany = 0.0;
  for(j=1;j<=nn;j++)
      meany += yy[j][1];
      
      meany = meany /nn;
  
  
  for(j = 1; j <= nn; j++) {
    double a = (float) rand()/(double)RAND_MAX*-0.1+0.05;    
       if(j >= nn/2)
          eta = 0;
           else
          eta =1;
                         
     for(i = 1; i <= nn; i++) {   
     F[j][1]=  yy[j][2] + eps*arr[j][i]*(yy[i][1] - yy[j][1]) + eps1*(meany - meanx) + lambda*yy[j][3] +eta*a;                                                                     
     F[j][2]=  -alpha*yy[j][2] -omega*yy[j][1] - beta*yy[j][1]*yy[j][1]*yy[j][1]; 
     F[j][3] = -meanx - meany - yy[j][3];
      }
   }
}







