#include <Header>


typedef double		 						FT;
typedef COPT::KernelTrait<FT,int>			kernel;
typedef kernel::Array 						Array;
typedef kernel::Vector 						Vector;
typedef kernel::Matrix 						Matrix;
typedef COPT::TripletBase<FT,int>		 	Triplet;
typedef COPT::SpMatrixBase<FT,int>			SpMatrix;

int main(int argc,char* argv[])
{
	std::cout<<"!!!!"<<std::endl; 
	int rows = 10;
	int cols = 10;
	int elesize = 10;
	int* rowind = new int[elesize];
	int* colptr = new int[cols+1];
	FT*	vals = new FT[elesize];
	for ( int i = 0 ; i < elesize ; ++ i )
	{
		rowind[i]=i;
		vals[i]=2.0;
	}
	for ( int i = 0 ; i <= cols ; ++ i )
		colptr[i] = i;
	SpMatrix mat(rows,cols,elesize,colptr,rowind,vals);
	delete[]rowind;
	delete[]colptr;
	delete[]vals;
	// std::cout<<mat(0,1)<<std::endl;
	// std::cout<<mat(1,1)<<std::endl;  
	// std::cout<<mat(9,9)<<std::endl;
	Vector vec(10);
	for ( int i = 0 ; i < 10 ; ++ i )
		vec[i] = i;
	std::cout<<mat*vec<<std::endl;

	std::vector<Triplet> tris;
	tris.push_back(Triplet(0,0,1));
	tris.push_back(Triplet(1,0,2));
	tris.push_back(Triplet(0,1,1)); 
	tris.push_back(Triplet(0,0,2));
	tris.push_back(Triplet(0,1,2));
	tris.push_back(Triplet(1,1,3));
	tris.push_back(Triplet(1,0,4));
	tris.push_back(Triplet(0,0,5));
	tris.push_back(Triplet(0,1,7));
	// for ( int i = 0 ; i < tris.size() ; ++ i )
	// {
	// 	std::cout<<tris[i].rowIndex()<<' '<<tris[i].columnIndex()<<std::endl;
	// }
	// std::sort(tris.begin(),tris.end(),COPT::columnComparison<Triplet>());
	// std::cout<<"after sorting"<<std::endl;
	// for ( int i = 0 ; i < tris.size() ; ++ i )
	// {
	// 	std::cout<<tris[i].rowIndex()<<' '<<tris[i].columnIndex()<<std::endl;
	// }

	SpMatrix m;
	m.setFromTriplets(10,10,tris.begin(),tris.end());
	Vector v(10);
	v(0) = 1.0;
	std::cout<<m*v<<std::endl;
	v(0) = 0.0;
	v(1) = 1.0;
	std::cout<<m*v<<std::endl;
	std::cout<<m.toDenseMatrix()<<std::endl;
	mat = 10*mat ;
	std::cout<<mat.toDenseMatrix()<<std::endl;

	// std::cout<<(m*mat).toDenseMatrix()<<std::endl;
	std::cout<<(m+mat).toDenseMatrix()<<std::endl;
	std::cout<<(m-mat).toDenseMatrix()<<std::endl; 
	std::cout<<(mat-m).toDenseMatrix()<<std::endl;
	std::cout<<((mat-m)*mat).toDenseMatrix()<<std::endl; 

	SpMatrix mm(m);
	std::cout<<mm.toDenseMatrix()<<std::endl;
}