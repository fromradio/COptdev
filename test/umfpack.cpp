// umfpack test

#include <stdio.h>
#include "umfpack.h"
#include <Header>

#define h long

h n = 5;
// COPT::COPTlong Ap[]={0,2,5,9,10,12};
COPT::COPTlong *Ap;
COPT::COPTlong Ai[]={0,1,0,2,4,1,2,3,4,2,1,4};
double Ax[]={2.,3.,3.,-1.,4.,4.,-3.,1.,2.,2.,6.,1.};
double* Axx = Ax;
double b[]={8.,45.,-3.,3.,19.};
double x[5];
double Info[UMFPACK_INFO],Control[UMFPACK_CONTROL];

int main(int argc,char*argv[])
{
	// double *null = (double*)NULL;
	Ap = new COPT::COPTlong[6];
	Ap[0]=0;
	Ap[1]=2;
	Ap[2]=5;
	Ap[3]=9;
	Ap[4]=10;
	Ap[5]=12;
	int i;
	void *Symbolic,*Numeric;
	COPT::umfpack_symbolic(n,n,Ap,Ai,Ax,&Symbolic,Control,Info);
	COPT::umfpack_numeric(Ap,Ai,Axx,Symbolic,&Numeric,Control,Info);
	COPT::umfpack_free_symbolic(&Symbolic,double());
	COPT::umfpack_solve(UMFPACK_A,Ap,Ai,Axx,x,b,Numeric,Control,Info);
	COPT::umfpack_free_numeric(&Numeric,double());
	for(i =0 ; i < n ; ++ i ) printf("x[%d] = %g\n",i,x[i]);
	umfpack_di_report_info(Control,Info);

	return 0;
}