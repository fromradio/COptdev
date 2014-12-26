// #include <umfpack.h>


// #include <Eigen/Dense>

// #include <Header>
// #include <IO>

// #include <Eigen/Sparse>
// #include <Eigen/src/UmfPackSupport/UmfPackSupport.h>
// #include <string>
// #include <vector>
// #include <fstream>
// #include <iostream>


// typedef Eigen::SparseMatrix<double>		SpMat;
// typedef Eigen::Triplet<double>			T;
// typedef Eigen::VectorXd 				Vector;
// typedef std::vector<T>					Triplets;


// typedef COPT::SpMatrixBase<double,int>	CSpMat;
// typedef COPT::VectorBase<double,int>	CVec;


// inline void readMtxFile( const std::string& filename , SpMat& mat )
// {
// 	if(filename.substr(filename.find_last_of(".")+1)=="mtx")
// 	{
// 		// deal the matrix
// 		std::ifstream fin(filename);
// 		if(!fin)
// 		{
// 			std::cerr<<"file does not exist!"<<std::endl;
// 			return;
// 		}
// 		std::string temp;
// 		int rows,cols,nnz;
// 		bool first = true;
// 		Triplets tris;
// 		while(fin)
// 		{
// 			std::getline(fin,temp);
// 			if(temp.c_str()[0] == '%')
// 				continue;
// 			else if(temp.empty())
// 				continue;
// 			else if(temp.c_str()[0] == ' ')
// 				continue;
// 			else
// 			{
// 				std::istringstream iss(temp);
// 				std::vector<std::string> tokens;
// 				std::copy(std::istream_iterator<std::string>(iss),
// 					std::istream_iterator<std::string>(),
// 					std::back_inserter(tokens));
// 				if( first )
// 				{
// 					std::cout<<"temp is "<<std::endl;
// 					if(tokens.size() == 3)
// 					{
// 						// matrix
// 						// if(sizeof(typename SpMatrix::Size)==4)
// 						// {
// 							rows = atoi(tokens[0].c_str());
// 							cols = atoi(tokens[1].c_str());
// 							nnz = atoi(tokens[2].c_str());
// 						// }
// 						// else
// 						// {
							
// 						// 	rows = atol(tokens[0].c_str());
// 						// 	cols = atol(tokens[1].c_str());
// 						// 	nnz = atol(tokens[2].c_str());
// 						// 	std::cout<<"rows is "<<rows<<std::endl;
// 						// }
// 						tris.reserve(nnz);
// 					}
// 					else
// 					{
// 						// throw COException("Error, the input should be a matrix not a vector!");
// 					}
// 					first = false;
// 				}
// 				else
// 				{
// 					// if(sizeof(typename SpMatrix::Size)==4)
// 					// {

// 						tris.push_back(T(atoi(tokens[0].c_str())-1,atoi(tokens[1].c_str())-1,atof(tokens[2].c_str())));
// 					// }
// 					// else
// 					// {
// 					// 	tris.push_back(T(atol(tokens[0].c_str())-1,atol(tokens[1].c_str())-1,atof(tokens[2].c_str())));
// 					// }
// 				}
// 			}
// 		}
// 		mat.resize(rows,cols);
// 		mat.setFromTriplets(tris.begin(),tris.end());
// 	}
// 	else
// 		return;
// }

int main(int argc ,char*argv[])
{
	// SpMat m;
	// Vector vec;
	// readMtxFile("data/ecology1.mtx",m);

	// // CSpMat mm;
	// // COPT::readMtxFile("data/bcsstk17.mtx",mm);
	// // std::cout<<"m's element size is "<<m.elementSize()<<std::endl;
	// Vector v;
	// // v.resize(m.rows());
	// v = Vector::Random(m.rows());
	// // initrand();
	// // for ( int i = 0 ; i < m.rows() ; ++ i )
	// // 	v[i] = randdouble();
	// vec = m*v;
	// // std::cout<<"v is "<<v<<std::endl;
	// // Eigen::SparseLU<SpMat> solver;
	// // solver.compute(m);
	// // Vector result = solver.solve(vec);
	// // Vector result = m.SimplicialLDLT().solve(vec);
	// // std::cout<<"result is "<<result<<std::endl;
	// // std::cout<<"ldlterror is "<<std::sqrt((result-v).squaredNorm())<<std::endl;
	// Eigen::UmfPackLU<SpMat> solver2;
	// solver2.compute(m);
	// Vector result = solver2.solve(vec);
	// // Vector result = m.SimplicialLDLT().solve(vec);
	// // std::cout<<"result is "<<result<<std::endl;
	// std::cout<<"umfpack error is "<<std::sqrt((result-v).squaredNorm())<<std::endl;

	// void* symbolic,*numeric;

	// // const double *control = new double[UMFPACK_CONTROL];
	// double *info = new double[UMFPACK_INFO];
	// int s1 = COPT::umfpack_symbolic(mm.rows(),mm.cols(),mm.columnPointer(),mm.rowIndex(),mm.values(),&symbolic,0,info);

	// std::cout<<s1<<std::endl;
	// int s2 = COPT::umfpack_numeric(mm.columnPointer(),mm.rowIndex(),mm.values(),symbolic,&numeric,0,info);
	// std::cout<<s2<<std::endl;
	// Vector vvv(vec.size());

	// CVec vv(v),vvec(vec),rr(vec.size());
	// int ec = COPT::umfpack_solve(UMFPACK_A,mm.columnPointer(),mm.rowIndex(),mm.values(),rr.dataPtr(),vvec.dataPtr(),numeric,0,info);

	// // std::cout<<"error is "<<std::sqrt((vvv-v).squaredNorm())<<std::endl;
	// std::cout<<"lala is "<<std::sqrt((vv-rr).squaredNorm())<<std::endl;

}