//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef MTX_FILE_HPP__
#define MTX_FILE_HPP__


namespace COPT
{

/************************Mtx File reader*********************/

template<class T>
inline void readMtxFile(const std::string&filename , T& t)
{
	readMtxFile(filename,t,typename T::ObjectCategory());
}


template<class Matrix>
inline void readMtxFile( const std::string& filename , Matrix& mat , const matrix_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)=="mtx")
	{
		std::ifstream fin(filename);
		if(!fin)
		{
			std::cerr<<"file does not exist!"<<std::endl;
			return;
		}
		std::string temp;
		typename Matrix::index rows,cols,i=0;
		bool first = true;
		while (fin)
		{
			std::getline(fin,temp);
			if(temp.c_str()[0] == '%' )
				continue;
			else if(temp.empty())
				continue;
			else if (temp.c_str()[0] == ' ')
				continue;
			else
			{
				std::istringstream iss(temp);
				std::vector<std::string> tokens;
				std::copy(std::istream_iterator<std::string>(iss),
					std::istream_iterator<std::string>(),
					std::back_inserter(tokens));
				if(first)
				{
					if(tokens.size()!=2)
					{
						std::cerr<<"the format of matrix market file is not right !"<<std::endl;
						return;
					}
					if(sizeof(typename Matrix::index)==4)
					{
						rows = atoi(tokens[0].c_str());
						cols = atoi(tokens[1].c_str());
					}
					else
					{
						rows = atol(tokens[0].c_str());
						cols = atol(tokens[1].c_str());
					}
					mat.resize(rows,cols);
					first = false;
				}
				else
				{
					mat.dataPtr()[i++]=atof(tokens[0].c_str());
				}
			}
		}
	}
	else
	{	
		std::cerr<<"Mtx file reader warning: File extension is not mtx!"<<std::endl;
		return;
	}
}

template<class SpMatrix>
inline void readMtxFile( const std::string& filename , SpMatrix& mat , const sp_matrix_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)=="mtx")
	{
		// deal the matrix
		std::ifstream fin(filename);
		if(!fin)
		{
			std::cerr<<"file does not exist!"<<std::endl;
			return;
		}
		std::string temp;
		typename SpMatrix::index rows,cols,nnz;
		bool first = true;
		std::vector<typename SpMatrix::Triplet> tris;
		while(fin)
		{
			std::getline(fin,temp);
			if(temp.c_str()[0] == '%')
				continue;
			else if(temp.empty())
				continue;
			else if(temp.c_str()[0] == ' ')
				continue;
			else
			{
				std::istringstream iss(temp);
				std::vector<std::string> tokens;
				std::copy(std::istream_iterator<std::string>(iss),
					std::istream_iterator<std::string>(),
					std::back_inserter(tokens));
				if( first )
				{
					if(tokens.size() == 3)
					{
						// matrix
						if(sizeof(typename SpMatrix::index)==4)
						{
							rows = atoi(tokens[0].c_str());
							cols = atoi(tokens[1].c_str());
							nnz = atoi(tokens[2].c_str());
						}
						else
						{
							
							rows = atol(tokens[0].c_str());
							cols = atol(tokens[1].c_str());
							nnz = atol(tokens[2].c_str());
						}
						tris.reserve(nnz);
					}
					else
					{
						throw COException("Error, the input should be a matrix not a vector!");
					}
					first = false;
				}
				else
				{
					if(sizeof(typename SpMatrix::index)==4)
					{

						tris.push_back(typename SpMatrix::Triplet(atoi(tokens[0].c_str())-1,atoi(tokens[1].c_str())-1,atof(tokens[2].c_str())));
					}
					else
					{
						tris.push_back(typename SpMatrix::Triplet(atol(tokens[0].c_str())-1,atol(tokens[1].c_str())-1,atof(tokens[2].c_str())));
					}
				}
			}
		}
		mat.fastSetFromTriplets(rows,cols,tris.begin(),tris.end());
		fin.close();
	}
	else
		return;
}

template<class Vector>
inline void readMtxFile( const std::string& filename , Vector& vec , const vector_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)=="mtx")
	{
		// deal the matrix
		std::ifstream fin(filename);
		if(!fin)
		{
			std::cerr<<"file does not exist!"<<std::endl;
			return;
		}
		std::string temp;
		typename Vector::index size,i = 0;
		bool first = true;
		while(fin)
		{
			std::getline(fin,temp);
			if(temp.c_str()[0] == '%')
				continue;
			else if(temp.empty())
				continue;
			else if(temp.c_str()[0] == ' ')
				continue;
			else
			{
				std::istringstream iss(temp);
				std::vector<std::string> tokens;
				std::copy(std::istream_iterator<std::string>(iss),
					std::istream_iterator<std::string>(),
					std::back_inserter(tokens));
				if( first )
				{
					if(tokens.size() == 2)
					{
						// matrix
						if(sizeof(typename Vector::index)==4)
						{
							if(atoi(tokens[1].c_str()) != 1 )
								throw COException("Error, the input is a matrix! Please check it out!");
							size = atoi(tokens[0].c_str());
						}
						else
						{
							if(atol(tokens[1].c_str()) != 1 )
								throw COException("Error, the input is a matrix! Please check it out!");
							size = atol(tokens[0].c_str());
						}
						vec.resize(size);
					}
					else
					{
						throw COException("Error, the input should be a vector not a sparse matrix!");
					}
					first = false;
				}
				else
				{
					vec[i] = atof(tokens[0].c_str());
					++ i;
				}
			}
		}
		fin.close();
	}
	else
		return;
}

/****************************Mtx file writer*********************************/
template<class T>
void writeMtxFile( const std::string& filename , const T& t )
{
	writeMtxFile(filename,t,typename T::ObjectCategory());
}

template<class Vector>
void writeMtxFile( const std::string& filename , const Vector& vec , const vector_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)!="mtx")
	{
		std::cerr<<"Mtx file writing warning: File extension is not mtx!"<<std::endl;
	}
	std::ofstream fout(filename);
	fout<<"%% MatrixMarket vector data"<<std::endl;
	fout<<"%% Generated by open source library COPT"<<std::endl;
	fout<<vec.size()<<' '<<1<<std::endl;
	for ( int i = 0 ; i < vec.size() ; ++ i )
	{
		fout<<vec(i)<<std::endl;
	}
}

template<class Matrix>
void writeMtxFile( const std::string& filename , const Matrix& mat , const matrix_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)!="mtx")
	{
		std::cerr<<"Mtx file writing warning: File extension is not mtx!"<<std::endl;
	}
	std::ofstream fout(filename);
	fout<<"%% MatrixMarket matrix data"<<std::endl;
	fout<<"%% Generated by open source library COPT"<<std::endl;
	fout<<mat.rows()<<' '<<mat.cols()<<std::endl;
	for ( int i = 0 ; i < mat.rows()*mat.cols() ; ++ i )
	{
		fout<<mat.dataPtr()[i]<<std::endl;
	}
}

template<class SpMatrix>
void writeMtxFile( const std::string& filename , const SpMatrix& spmat , const sp_matrix_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)!="mtx")
	{
		std::cerr<<"Mtx file writing warning: File extension is not mtx!"<<std::endl;
	}
	std::ofstream fout(filename);
	fout<<"%% MatrixMarket sparse matrix data"<<std::endl;
	fout<<"%% Generated by open source library COPT"<<std::endl;
	fout<<spmat.rows()<<' '<<spmat.cols()<<' '<<spmat.elementSize()<<std::endl;

}

}

#endif