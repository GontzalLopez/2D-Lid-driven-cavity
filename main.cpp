/*
2D Lid-driven cavity problem with FVM
SIMPLE algorithm, staggered grid, Hybrid scheme, SOR algorithm for pressure calculation
Theory: H.Versteeg W.Malalasekra "An introduction to CFD"

Gontzal Lopez Ruiz / Sept-2021
*/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//Grid size 
#define  nx 100
#define  ny 100

//*****Aux. functions*****//
//************************//

void setOne(double A[][nx+2])
{
	for(int i=0;i<(nx+1);i++)
	{
		for(int j=0;j<(ny+1);j++)
		{
			A[i][j]=1.0;
		}
	}
}

double max3(double a, double b, double c)
{
	return (a>b) ? ((a>c) ? a : c) : ((b>c) ? b : c);
}

//************************//
//************************//


int main()
{
//Cavity dimensions
double Lx=1; 
double Ly=1;
double dx=Lx/nx;
double dy=Ly/ny;

double mu=0.01;
double rho=1.0;

//relaxation parameters
double alphap = 0.4;
double alpha = 0.8;
int itP=200;
int uvIt=12;

double contError;
double Fw_x,Fe_x,Fs_x,Fn_x,Dw_x,De_x,Ds_x,Dn_x;
double Fw_y,Fe_y,Fs_y,Fn_y,Dw_y,De_y,Ds_y,Dn_y;
double ae[ny+2][nx+2]={0.0};
double aw[ny+2][nx+2]={0.0};
double an[ny+2][nx+2]={0.0};
double as[ny+2][nx+2]={0.0};
double aPu[ny+2][nx+2];
double aPv[ny+2][nx+2];
double aPp[ny+2][nx+2];
setOne(aPu);
setOne(aPv);
setOne(aPp);

int un=1.0; //north inlet BC [m/s]
double u[ny+2][nx+1]={0.0};
double uStar[ny+2][nx+1]={0.0};
double uNew[ny+2][nx+1]={0.0};
double uFinal[ny+1][nx+1]={0.0};

double v[ny+1][nx+2]={0.0};
double vStar[ny+1][nx+2]={0.0};
double vNew[ny+1][nx+2]={0.0};
double vFinal[ny+1][nx+1]={0.0};

double du[ny+2][nx+1]={0.0};
double dv[nx+1][nx+2]={0.0};


double p[ny+2][nx+2];
double pc[ny+2][nx+2]={0.0};
setOne(p);
double pFinal[nx+1][ny+1]={0.0};
double b[nx+2][ny+2]={0.0};
	
//coordinates&directions:

//i=0 top row
//i=max bottom row
//j=0 left column 
//j=max rigth column 

//u velocity inlet BC
for(int j = 0; j<(nx+1); j++)
{
	u[0][j]=1;     // north inlet
	uStar[0][j]=1;     // north inlet		
}

// Start loop

cout<<"Starting time loop..."<<endl;

double error=1;
double error_req=1e-07;
int itnum=1;

while(error>error_req)
{
	cout<<"Iteration Number : "<<itnum<<endl;

		//x momentum
		for(int i=1;i<ny+1;i++)
		{
			for(int j=1;j<nx;j++)
			{
				Fw_x=rho*0.5*dy*(u[i][j]+u[i][j-1]);
				Fe_x=rho*0.5*dy*(u[i][j+1]+u[i][j]);
				Fs_x=rho*0.5*dx*(v[i][j]+v[i][j+1]);
				Fn_x=rho*0.5*dx*(v[i-1][j]+v[i-1][j+1]);
				Dw_x=(mu*dy)/dx;
				De_x=(mu*dy)/dx;
				Ds_x=(mu*dx)/dy;
				Dn_x=(mu*dx)/dy;

				//Hybrid scheme
				ae[i][j]=max3(0.0,-Fe_x,(De_x-(0.5*Fe_x)));
				aw[i][j]=max3(0.0,Fw_x,(Dw_x+(0.5*Fw_x)));
				an[i][j]=max3(0.0,-Fn_x,(Dn_x-(0.5*Fn_x)));
				as[i][j]=max3(0.0,Fs_x,(Ds_x+(0.5*Fs_x)));
				aPu[i][j]=aw[i][j]+ae[i][j]+an[i][j]+as[i][j]+(Fe_x-Fw_x)+(Fn_x-Fs_x);
				aPu[i][j]=aPu[i][j]/alpha;
			}
		}
	
		for(int k=1;k<uvIt;k++)
		{
			for(int i=1;i<ny+1;i++)
			{
				for(int j=1;j<nx;j++)
				{
					uStar[i][j] =(1-alpha)*u[i][j]+(1.0/aPu[i][j])*(ae[i][j]*uStar[i][j+1]+aw[i][j]*uStar[i][j-1]+an[i][j]*uStar[i-1][j]+as[i][j]*uStar[i+1][j]+dy*(p[i][j+1]-p[i][j]));
				}
			}
		}
	
		//y momentum
		for(int i=1;i<ny;i++)
		{
			for(int j=1;j<nx+1;j++)
			{
				
				Fw_y=rho*0.5*dx*(u[i][j-1]+u[i+1][j-1]);
				Fe_y=rho*0.5*dx*(u[i][j]+u[i+1][j]);
				Fs_y=rho*0.5*dx*(v[i][j]+v[i+1][j]);
				Fn_y=rho*0.5*dx*(v[i-1][j]+v[i][j]);
				Dw_y=(mu*dx)/dy;
				De_y=(mu*dx)/dy;
				Ds_y=(mu*dy)/dx;
				Dn_y=(mu*dy)/dx;

				//Hybrid scheme
				ae[i][j]=max3(0.0,-Fe_y,(De_y-(0.5*Fe_y)));
				aw[i][j]=max3(0.0,Fw_y,(Dw_y+(0.5*Fw_y)));
				an[i][j]=max3(0.0,-Fn_y,(Dn_y-(0.5*Fn_y)));
				as[i][j]=max3(0.0,Fs_y,(Ds_y+(0.5*Fs_y)));
				
				aPv[i][j]=an[i][j]+as[i][j]+aw[i][j]+ae[i][j]+(Fe_y-Fw_y)+(Fn_y-Fs_y);
				aPv[i][j]=aPv[i][j]/alpha;
			}
		}
	
		for(int k=1;k<uvIt;k++)
		{		
			for(int i=1;i<ny;i++)
			{
				for(int j=1;j<nx+1;j++)
				{
					vStar[i][j] =(1-alpha)*v[i][j]+(1.0/aPv[i][j])*(ae[i][j]*vStar[i][j+1]+aw[i][j]*vStar[i][j-1]+an[i][j]*vStar[i-1][j]+as[i][j]*vStar[i+1][j]+dx*(p[i][j]-p[i+1][j]));					
				}
			}
		}
       
		//Mass imbalance
		error=0;
			for(int i=1;i<ny+1;i++)
			{
				for(int j=1;j<nx+1;j++)
				{
					b[i][j]=(uStar[i][j]-uStar[i][j-1])*dy+(vStar[i-1][j]-vStar[i][j])*dx;			
                    			error=error+b[i][j]*b[i][j];					
				}
			}

		//continuity residual
		error=sqrt(error);

		cout<<"Continuity error: "<<error<<endl;
		
			for(int i=1;i<ny+1;i++)
			{
				for(int j=1;j<nx+1;j++)
				{
					ae[i][j]=(dx*dy)/aPu[i][j];
					aw[i][j]=(dx*dy)/aPu[i][j-1];
					an[i][j]=(dy*dx)/aPv[i-1][j];
					as[i][j]=(dy*dx)/aPv[i][j];
					aPp[i][j]=ae[i][j]+aw[i][j]+an[i][j]+as[i][j];
				}
			}
					
			
		//SOR algorithm with w=1.7 for pressure
       		 for(int k=1;k<=itP;k++)
		{
			for(int i=1;i<ny+1;i++)
			{
				for(int j=1;j<nx+1;j++)
				{		
					pc[i][j]=pc[i][j]+(1.7/aPp[i][j])*(ae[i][j]*pc[i][j+1]+aw[i][j]*pc[i][j-1]+an[i][j]*pc[i-1][j]+as[i][j]*pc[i+1][j]+b[i][j]-pc[i][j]*aPp[i][j]);
				}
			}
		}

		//Pressure correction
       	        for(int i=1;i<ny+1;i++)
		{
           		 for(int j=1;j<nx+1;j++)
			{
				p[i][j]=p[i][j]+alphap*pc[i][j];
			}
		}

		//u velocity correction
		for(int i=1;i<ny+1;i++)
		{
			for(int j=1;j<nx;j++)
			{
				uStar[i][j]=uStar[i][j]+(dy/aPu[i][j])*(pc[i][j+1]-pc[i][j]);
			}
		}

		//v velocity correction
		for(int i=1;i<ny;i++)
		{
			for(int j=1;j<nx+1;j++)
			{
				vStar[i][j]=vStar[i][j]+(dx/aPv[i][j])*(pc[i][j]-pc[i+1][j]);
			}
		}

		//
		//update velocity values
		for(int i=0;i<ny+2;i++)
		{
			for(int j=0;j<nx+1;j++)
			{
				u[i][j]=uStar[i][j];
			}
		}

		for(int i=0;i<ny+1;i++)
		{
			for(int j=0;j<nx+2;j++)
			{
				v[i][j]=vStar[i][j];
			}
		}
itnum++;
}//end iteration

cout<<"End calculation"<<endl;

//************Post-process results************//
//********************************************//
//********************************************//

//Velocity interpolation

for(int i=0;i<ny+1;i++)
{
	for(int j=0;j<nx+1;j++)
	{
	 uFinal[i][j]=0.5*(u[i][j]+u[i+1][j]);
	 vFinal[i][j]=0.5*(v[i][j]+v[i][j+1]);
	 pFinal[i][j]=0.25*(p[i][j]+p[i][j+1]+p[i+1][j]+p[i+1][j+1]);
	}
}

//Output files

FILE *fdata;
fdata = fopen("./post/dataCavity.dat","w+t");

if ( fdata == NULL )
{
	printf("\nERROR when opening file\n");
	fclose( fdata );
}
else
{
//	fprintf( fdata, "Xcoord   Ycoord  Uvelocity\n");

	for (int j = 0;j<(nx+1);j++)
	{
	double xCoord=j*dx;

		for (int i = ny;i>=0;i--)
		{
			double yCoord=(ny-i)*dy;
			fprintf( fdata, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xCoord,yCoord,uFinal[i][j],vFinal[i][j],pFinal[i][j]);
		}
		fprintf(fdata,"\n");
	}	
}

fclose(fdata);

    // OUTPUT DATA
//FILE *fdata;
fdata = fopen("./post/dataCavity_prof.dat","w+t");

if ( fdata == NULL )
{
	printf("\nERROR when opening file\n");
	fclose( fdata );
}
else
{
//	fprintf( fdata, "Xcoord   Ycoord  Uvelocity\n");

	for (int j = 0;j<(nx+1);j++)
	{
	double xCoord=j*dx;
		fprintf( fdata, "%5.8lf\t%5.8lf\t%5.8lf\n", xCoord,uFinal[50][j],vFinal[50][j]);
	}	
}

fclose(fdata);


//FILE *fdata;
fdata = fopen("./post/dataCavity_uExp.dat","w+t");

if ( fdata == NULL )
{
	printf("\nERROR when opening file\n");
	fclose( fdata );
}
else
{
//	fprintf( fdata, "Xcoord   Ycoord  Uvelocity\n");

	for (int i = ny;i>=0;i--)
	{
		double yCoord=(ny-i)*dy;
		fprintf( fdata, "%5.8lf\t%5.8lf\n", yCoord,uFinal[i][50]);
	}	
}

fclose(fdata);
}

