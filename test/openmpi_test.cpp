//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#include <mpi.h>
#include <Header>

typedef COPT::KernelTrait<double> 		kernel;
typedef kernel::Vector 					Vector;
typedef kernel::scalar					scalar;
typedef kernel::index 					Index;

const int sizevec = 4;
void summation( const int id , const int p );
/*		Compute the summation of two vectors
 *		The first processor is the main processor
 *
 */
int main( int argc , char* argv[])
{

	int id, p;
	

	MPI::Init(argc,argv);
	MPI::Status status;

	id = MPI::COMM_WORLD.Get_rank();
	p = MPI::COMM_WORLD.Get_size();
	// std::cout<<"size is "<<sizeof(scalar)<<' '<<sizeof(MPI::)<<' '<<sizeof(double)<<std::endl;
	// if(id==0)
	// {
	// 	scalar* s = new scalar[2];
	// 	s[0]=1.0;
	// 	s[1]=2.0;
	// 	MPI::COMM_WORLD.Send(&s[0],2,MPI::DOUBLE,1,1);
	// }
	// else if(id==1)
	// {
	// 	scalar* ss = new scalar[2];
	// 	MPI::COMM_WORLD.Recv(ss,2,MPI::DOUBLE,0,1,status);
	// 	std::cout<<ss[0]<<std::endl;
	// 	std::cout<<ss[1]<<std::endl;
	// }
	// std::cout<<id<<" "<<p<<std::endl;
	summation(id,p);
	MPI::Finalize();
}

void summation( const int id , const int p )
{
	scalar *s1 = NULL,*s2 = NULL;
	scalar *result = NULL;
	Index l;
	MPI::Status status;
	if ( p > 1 )
	{
		if (id < p-1 )
			l = sizevec/(p-1);
		else
			l = sizevec-sizevec/(p-1)*(p-2);
	}
	if ( p == 1 )
		// only one processor exists
	{
		Vector v1 = Vector::random(sizevec);
		Vector v2 = Vector::random(sizevec);
		std::cout<<"the first vector is "<<v1<<std::endl;
		std::cout<<"the second vector is "<<v2<<std::endl;
		std::cout<<"the summation is "<<v1+v2<<std::endl;
	}
	else
	{
		if( id == 0 )
		{
			Vector v1 = Vector::random(sizevec);
			Vector v2 = Vector::random(sizevec);
			std::cout<<"the first vector is "<<v1<<std::endl;
			std::cout<<"the second vector is "<<v2<<std::endl;
			s1 = v1.dataPtr();
			s2 = v2.dataPtr();
			int t1 = 1;
			int t2 = 2;
			for ( int i = 1 ; i < p-1 ; ++ i )
			{
				MPI::COMM_WORLD.Send(&s1[sizevec*(i-1)],sizevec/(p-1),MPI::DOUBLE,i,t1);
				MPI::COMM_WORLD.Send(&s2[sizevec*(i-1)],sizevec/(p-1),MPI::DOUBLE,i,t2);
			}
			MPI::COMM_WORLD.Send(&s1[sizevec*(p-2)],sizevec-sizevec/(p-1)*(p-2),MPI::DOUBLE,p-1,t1);
				MPI::COMM_WORLD.Send(&s2[sizevec*(p-2)],sizevec-sizevec/(p-1)*(p-2),MPI::DOUBLE,p-1,t2);
		}
		else
		{
			// std::cout<<"processor "<<id<<" with length "<<l<<std::endl;
			int t1 = 1,t2 = 2;
			s1 = new scalar[l];
			s2 = new scalar[l];
			MPI::COMM_WORLD.Recv(s1,l,MPI::DOUBLE,0,t1,status);
			MPI::COMM_WORLD.Recv(s2,l,MPI::DOUBLE,0,t2,status);
			std::cout<<id<<" receive success!"<<std::endl;
			result = new scalar[l];
			for ( int i = 0 ; i < l ; ++ i )
			{
				result[i] = s1[(id-1)*l+i]+s2[(id-1)*l+i];
			}
			int tag = 3;
			MPI::COMM_WORLD.Send(result,l,MPI::DOUBLE,0,tag);
		}

		if ( id == 0 )
		{
			Vector result(sizevec);
			int length=sizevec/(p-1);
			int tag = 3;
			for ( int i = 1 ; i < p - 1 ; ++ i )
			{
				MPI::COMM_WORLD.Recv(&result.dataPtr()[(i-1)*length],length,MPI::DOUBLE,i,tag,status);
			}
			MPI::COMM_WORLD.Recv(&result.dataPtr()[(p-2)*length],sizevec-(p-2)*length,MPI::DOUBLE,p-1,tag,status);
			std::cout<<"the summation is "<<result<<std::endl;
		}
	}
}

