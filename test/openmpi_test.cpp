//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#include <mpi.h>
#include <Header>

typedef COPT::KernelTrait<double> 		kernel;
typedef kernel::Vector 					Vector;
typedef kernel::scalar					scalar;
typedef kernel::index 					Index;

const int sizevec = 10000000;
void summation( const int id , const int p );
/*		Compute the summation of two vectors
 *		The first processor is the main processor
 *
 */
int main( int argc , char* argv[])
{

	int id, p;
	double wtime;

	MPI::Init(argc,argv);
	MPI::Status status;

	id = MPI::COMM_WORLD.Get_rank();
	p = MPI::COMM_WORLD.Get_size();

	if (id==0)
		wtime = MPI::Wtime();

	summation(id,p);

	if (id==0)
		wtime = MPI::Wtime()-wtime;

	if (id==0)
		std::cout<<"time elapsed: "<<wtime<<"s."<<std::endl;
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
		result = new scalar[sizevec];
		for ( int i = 0 ; i < sizevec ; ++ i )
			result[i] = v1(i)+v2(i);
		Vector r(sizevec,result);
	}
	else
	{
		if( id == 0 )
		{
			Vector v1 = Vector::random(sizevec);
			Vector v2 = Vector::random(sizevec);
			// std::cout<<"the first vector is "<<v1<<std::endl;
			// std::cout<<"the second vector is "<<v2<<std::endl;
			s1 = v1.dataPtr();
			s2 = v2.dataPtr();
			int t1 = 1;
			int t2 = 2;
			std::cout<<"send begins"<<std::endl;
			for ( int i = 1 ; i < p-1 ; ++ i )
			{
				// std::cout<<s1[sizevec/(p-1)*(i-1)]<<std::endl;
				MPI::COMM_WORLD.Send(&s1[sizevec/(p-1)*(i-1)],sizevec/(p-1),MPI::DOUBLE,i,t1);
				MPI::COMM_WORLD.Send(&s2[sizevec/(p-1)*(i-1)],sizevec/(p-1),MPI::DOUBLE,i,t2);
			}
			MPI::COMM_WORLD.Send(&s1[sizevec/(p-1)*(p-2)],sizevec-sizevec/(p-1)*(p-2),MPI::DOUBLE,p-1,t1);
			MPI::COMM_WORLD.Send(&s2[sizevec/(p-1)*(p-2)],sizevec-sizevec/(p-1)*(p-2),MPI::DOUBLE,p-1,t2);
			std::cout<<"send over!"<<std::endl;
			result = new scalar[sizevec];
		}
		MPI::COMM_WORLD.Bcast(result,sizevec,MPI::DOUBLE,0);
		if( id > 0 )
		{
			int t1 = 1,t2 = 2;
			s1 = new scalar[l];
			s2 = new scalar[l];
			MPI::COMM_WORLD.Recv(s1,l,MPI::DOUBLE,0,t1,status);
			MPI::COMM_WORLD.Recv(s2,l,MPI::DOUBLE,0,t2,status);
			// std::cout<<id<<" receive success!"<<std::endl;
			result = new scalar[l];
			for ( int i = 0 ; i < l ; ++ i )
			{
				std::cout<<id<<" "<<i<<std::endl;
				result[(id-1)*(sizevec/(p-1))+i] = s1[i]+s2[i];
				// std::cout<<result[i]<<std::endl;
			}
			// int tag = 3;
			// MPI::COMM_WORLD.Send(result,l,MPI::DOUBLE,0,tag);
			delete[] s1;
			delete[] s2;
		}
	}
}

