#ifndef BCDSOLVER
#define BCDSOLVER

/** 
 * BCD solver
 * 
 * 
 * general form of the problem: 
 *      
 *          f(x) + sum_{i=1}^{s}r_i(x)
 * 
 * where f(x) is differentiable while r_i(x) is not 
 * 
 */


#include "Core"

typedef COPT::KernelTrait<double> 				kernel;
typedef kernel::Vector 							Vector;
typedef kernel::Matrix 							Matrix;
typedef COPT::BasicParameter<double> 			BasicParameter;
typedef COPT::BasicOption<double> 				BasicOption;

#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <functional>
using namespace std;

typedef kernel::podscalar Float;

//////////////////////////////////////////////////////

/// r - the non-smooth part
struct FuncRInfo
{
    // type of function r
    enum RType { Norm1, Indicator };
    
    RType   type;
    int     dim;
    
    FuncRInfo(RType tp = Norm1, int dm = 1);
};

/// block control
class BlockCtrl
{
public:
    
    inline void setBlocking(const vector<int>& blocking);
    inline int  nBlock()const;
    inline int  blockDimension(int i)const;
    inline int  fullDimension()const;
    inline void blockIdexInterval(int i, int& a, int& b)const;
    inline Vector getVectorBlock(const Vector& x0, int kblock)const;
    inline void   setVectorBlock(Vector& x0, int kblock, const Vector& x)const;
    
private:
    bool __is_set = false;
    int __nblock  = 0;
    vector<int> __block;
    vector<int> __blockIdx;
};

/// solver option
struct BCDSolverOption : public BasicOption
{
    /// added member
    Float XTol  = 1e-6;   // tolerance of change in x
    Float ObjTol= 1e-6;   // tolerance of change in objective
};

/// solver parameters
class BCDSolverParam : public BasicParameter
{
public:

    /// objective
    // tau
    Float __tau;
    // f - the differentiable part
    function<Float(const Vector&)>  __fEval;
    function<Vector(const Vector&)> __fGrad;

private:    
    /// ri
    // vector of ri
    vector<FuncRInfo> __ri;
    // evaluation
    Float riEval(const FuncRInfo& r, const Vector& x)const;

public:
    /// iteration
    int __verbose;
    int __k_iter;
    
    /// intermediate variables
    Vector    __x_prev;
    Float     __obj, __obj_prev;
    BlockCtrl __block_ctrl;
    
    /// auto updated inner param
    struct ParamAuto
    {
        Float t, w, L;
    };
    vector<ParamAuto> __param_auto, __param_auto_prev;
    
public:
    
    /// defaults
    BCDSolverParam();
    BCDSolverParam(
            const function<Float(const Vector&)>& fEval,
            const function<Vector(const Vector&)>& fGrad,
            const vector<FuncRInfo>& ri,
            Float tau);
    
    void setF(
            const function<Float(const Vector&)>& fEval,
            const function<Vector(const Vector&)>& fGrad);
    void setR(const vector<FuncRInfo>& ri);
    void setTau(Float tau);
    
    inline Float    fEval(const Vector& x)const;
    inline Vector   fGrad(const Vector& x)const;
    inline Float    rEval(const Vector& x, int kblock)const;
};


/// bcd solver
typedef COPT::Solver<Float,Vector,Vector,BCDSolverOption,BCDSolverParam> tagBCDSolver;    
class BCDSolver : public tagBCDSolver
{
public:
    BCDSolver(const BCDSolverParam& bcd_par,
              const BCDSolverOption& bcd_opt);
};

#include "BCDImpl.hpp"

#endif // BCDSOLVER
