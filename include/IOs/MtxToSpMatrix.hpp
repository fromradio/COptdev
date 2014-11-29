//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef MTX_TO_SP_MATRIX_HPP
#define MTX_TO_SP_MATRIX_HPP

namespace COPT
{
template<class SpMatrix>
inline void readMtxFile( const std::string& filename , SpMatrix& mat , const matrix_tag& )
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
		typename SpMatrix::Size rows,cols,nnz;
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
					std::cout<<"temp is "<<std::endl;
					if(tokens.size() == 3)
					{
						// matrix
						if(sizeof(typename SpMatrix::Size)==4)
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
							std::cout<<"rows is "<<rows<<std::endl;
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
					if(sizeof(typename SpMatrix::Size)==4)
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
		mat.setFromTriplets(rows,cols,tris.begin(),tris.end());
	}
	else
		return;
}

template<class Vector>
inline void readMtxFile( const std::string& filename , Vector& vec , const vector_tag& )
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
		typename Vector::Size size,i = 0;
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
						if(sizeof(typename Vector::Size)==4)
						{
							size = atoi(tokens[0].c_str());
						}
						else
						{
							size = atol(tokens[0].c_str());
						}
						vec.resize(size);
					}
					else
					{
						throw COException("Error, the input should be a vector not a matrix!");
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
	}
	else
		return;
}

template<class T>
inline void readMtxFile(const std::string&filename , T& t)
{
	readMtxFile(filename,t,typename T::Category());
}

} // End of namespace COPT

#endif