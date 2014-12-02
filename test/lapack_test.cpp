#include <stdio.h>
#include <clapack.h>
#define size 3

int main(int argc,char* argv[])
{
	int i,j,c1,c2,pivot[size],ok;
	float A[size][size],b[size],AT[size*size];

	A[0][0] = 3.1; A[0][1]=1.3; A[0][2] = -5.7;
	A[1][0] = 1.0; A[1][1] = -6.9; A[1][2] = 5.8;
	A[2][0] = 3.4; A[2][1] = 7.2; A[2][2] = -8.8;

	b[0] = -1.3;
	b[1] = -0.1;
	b[2] = 1.8;

	for ( i = 0 ; i < size ; ++ i )
	{
		for ( j = 0 ; j < size ; ++ j )
			AT[j+size*i]=A[j][i];
	}

	c1 = size;
	c2 = 1;

	sgesv_(&c1,&c2,AT,&c1,pivot,b,&c1,&ok);

	for ( j = 0 ; j < size ; ++ j ) printf("%e\n",b[j]);
}