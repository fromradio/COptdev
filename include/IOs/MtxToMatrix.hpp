//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef MTX_TO_MATRIX_HPP__
#define MTX_TO_MATRIX_HPP__

namespace COPT
{
template<class Matrix>
inline void readMtxFile( const std::string& filename , Matrix& mat , const matrix_tag& )
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
		return;
}

template<class Vector>
inline void writeMtxFile(const std::string& filename )
{}


}
#endif