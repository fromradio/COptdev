#include <Header>


typedef double		 					FT;
typedef COPT::Array<FT> 				Array;
typedef COPT::VectorBase<FT>			Vector;
typedef COPT::MatrixBase<FT>	 		Matrix;
typedef COPT::KernelTrait<FT>			kernel;
typedef COPT::TripletBase<FT>		 	Triplet;
typedef COPT::SpMatrixBase<FT>			SpMatrix;

int main(int argc,char* argv[])
{
	std::cout<<"!!!!"<<std::endl; 
	size_t rows = 10;
	size_t cols = 10;
	size_t elesize = 10;
	size_t* rowind = new size_t[elesize];
	size_t* colptr = new size_t[cols+1];
	FT*	vals = new FT[elesize];
	for ( int i = 0 ; i < elesize ; ++ i )
	{
		rowind[i]=i;
		vals[i]=2.0;
	}
	for ( int i = 0 ; i <= cols ; ++ i )
		colptr[i] = i;
	SpMatrix mat(rows,cols,elesize,rowind,colptr,vals);
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
	m.setFromTriplets(10,10,tris);
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