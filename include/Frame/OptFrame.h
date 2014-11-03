// Optimization Framework
// Powered by 'MathU'
// 	Leaded by Zhouwang Yang
//	Implemented by ...

#ifndef OPTFRAME_H
#define OPTFRAME_H

#include "BasicMath.h"
#include <iostream>
#include <vector>

namespace COpt
{

class VectorOptimization
{
	// a general class of vector optimization

public:
	// The type of terminal
	// 	Error: iteration terminals if error is less than threshold
	//	IterMax: iteration terminals if iteration number is bigger than threshold
	//	Both: the iteration terminals no matter which condition is satisfied
	enum Terminal{
		Error,
		IterMax,
		Both
	};

	class Parameter
	{
	public:
		// Constructor
		Parameter();

		// getter and setter
		int iterationThreshold();
		void setIterationThreshold( int  iter);
		double errorThreshold();
		void setErrorThreshold(double e);

	private:
		// threshhold of iteration number
		int __iterthresh;
		// threshhold of error
		double __errorthresh;
	};

	class OptOutput
	{
	public:
		// Constructor
		OptOutput();
		// append one error
		void append(double e);
	private:
		// storing the final iteration number
		int __iterationnumber;
		// the vector storing every error
		std::vector<double> __errorvec;
	};

	// Constructor
	VectorOptimization();

	// Deconstructor
	~VectorOptimization();

protected:

	// the target vector
	Vector __v;
	// parameter
	Parameter __p;
	// Terminal type
	Terminal __terminal;
	// Output
	OptOutput __output;

};
};

#endif