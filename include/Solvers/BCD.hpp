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


inline void get_vec_block(const Vector& x0, int id1, int id2, Vector& x);
inline void set_vec_block(Vector& x0, int id1, int id2, const Vector& x);
inline Vector get_vec_block(const Vector& x0, int id1, int id2);
inline Float shrinkage(Float z, Float tau);


//////////////////////////////////////////////////////
/// block control
class BlockCtrl
{
public:
    
    void setBlocking(const vector<int>& blocking)
    {
        m_block = blocking;
        m_nblock = m_block.size();
        m_blockIdx.resize(m_nblock+1);
        m_blockIdx[0] = 0;
        for(int i = 0 ; i < m_nblock; ++i)
        {
            m_blockIdx[i+1] = m_block[i];
            m_blockIdx[i+1] += m_blockIdx[i];
        }
        is_set = true;
    }
    
    int n_block()const {return m_nblock;}
    int block_dim(int i)const{return m_block[i];}
    int full_dim()const
    {
        int dim = 0;
        for(const auto& t:m_block) dim += t;
        return dim;
    }
    void block_id_interval(int i, int& a, int& b)const
    {
        a = m_blockIdx[i];
        b = m_blockIdx[i+1];
    }
    
    inline Vector get_vec_block(const Vector& x0, int kblock)const
    {
        assert(is_set);
        int ia, ib;
        block_id_interval(kblock, ia, ib);
        return ::get_vec_block(x0, ia, ib);
    }
    inline void set_vec_block(Vector& x0, int kblock, const Vector& x)const
    {
        assert(is_set);
        int ia, ib;
        block_id_interval(kblock, ia, ib);
        ::set_vec_block(x0, ia, ib, x);
    }
    
    bool is_set = false;
    
private:
    vector<int> m_block;
    vector<int> m_blockIdx;
    int m_nblock = 0;
};

/// solver option
struct BCDSolverOption : public BasicOption
{
///   from base class    
//    int   MaxIter;       // 
//    Float Threshold;     // 
    
    /// added member
    Float xTol;     // tolerance of change in x
    Float objTol;   // tolerance of change in objective
    
    BCDSolverOption()
        :   BasicOption()
        ,   xTol(1e-6)
        ,   objTol(1e-6)
    {}
};

// r - the non-smooth part
struct FuncRInfo
{
    // type of function r
    enum RType { Norm1, Indicator };
    
    RType type;
    int dim;
    
    FuncRInfo(RType tp = Norm1, int dm = 1)
        :   type(tp)
        ,   dim(dm)
    {}
};

/// solver parameters
class BCDSolverParam : public BasicParameter
{
public:

    /// objective
    // tau
    Float m_tau;
    // f - the differentiable part
    function<Float(const Vector&)>  m_fEval;
    function<Vector(const Vector&)> m_fGrad;

    
private:    
    /// ri
    // vector of ri
    vector<FuncRInfo> m_ri;
    // evaluation
    Float ri_eval(const FuncRInfo& r, const Vector& x)const
    {
        switch(r.type)
        {
        case FuncRInfo::RType::Norm1:
            return x.absNorm();
            break;
        case FuncRInfo::RType::Indicator:
            return 0;
            break;
        }
    }
    
public:
    /// iteration
    int verbose;
    int kIter;
    
    /// intermediate variables
    Vector xprev;
    Float obj, objprev;

    BlockCtrl m_blockCtrl;
    
    /// auto updated inner param
    struct paramAuto
    {
        Float t, w, L;
    };
    vector<paramAuto> m_paramAuto, m_paramAutoPrev;
    
public:
    
    /// defaults
    BCDSolverParam()
        :   BasicParameter()
        ,   verbose(0)
        ,   m_tau(1.0)
        ,   kIter(0)
    {}
    
    BCDSolverParam(
            const function<Float(const Vector&)>& f_eval,
            const function<Vector(const Vector&)>& f_grad,
            const vector<FuncRInfo>& ri,
            Float tau
            )
        :   BasicParameter()
        ,   verbose(0)
        ,   m_tau(1.0)
        ,   kIter(0)
    {
        set_f(f_eval, f_grad);
        set_r(ri);
        set_tau(tau);
    }
    
    void set_f(
            const function<Float(const Vector&)>& f_eval,
            const function<Vector(const Vector&)>& f_grad)
    {
        this->m_fEval = f_eval;
        this->m_fGrad = f_grad;
    }
    
    void set_r(const vector<FuncRInfo>& ri)
    {
        this->m_ri = ri;
        vector<int> blocking(ri.size());
        for(int i = 0; i < ri.size(); ++i)
        {
            blocking[i] = ri[i].dim;
        }
        m_blockCtrl.setBlocking(blocking);
    }
    
    void set_tau(Float tau)
    {
        this->m_tau = tau;
    }
    
    inline Float    f_eval(const Vector& x)const{return m_fEval(x);}
    inline Vector   f_grad(const Vector& x)const{return m_fGrad(x);}
    inline Float    r_eval(const Vector& x, int kblock)const
    {
        return ri_eval(m_ri[kblock], m_blockCtrl.get_vec_block(x,kblock));
    }
    
};

void updateL(const Vector& x, BCDSolverParam& par);
void paramInit(const Vector& x, BCDSolverParam& par);
void paramUpdate(const Vector& x, BCDSolverParam& par);
void solve_subblock(Vector& x, BCDSolverParam& par, int kblock);
Float errorEst(const Vector& x, const BCDSolverParam& par);

//////////////////////

Float bcd_objective(const Vector& x, const BCDSolverParam& par)
{
    Float obj = 0;
    // r part
    for(int i = 0; i < par.m_blockCtrl.n_block(); ++i)
        obj += par.r_eval(x, i);
    // f part
    obj += par.f_eval(x);
    return obj;
}

///
Float bcd_oneIteration(Vector& x, BCDSolverParam& par)
{
    par.xprev = x;
    par.objprev = par.obj;
    // update x block-wise
    auto nblock = par.m_blockCtrl.n_block();
    for(int kb = 0; kb < nblock; ++kb)
    {
        solve_subblock(x, par, kb);
    }
    
    //
    ++par.kIter;
    paramUpdate(x, par);
    par.obj = bcd_objective(x, par);
    
    return errorEst(x, par);
}

bool bcd_terminate(const Vector& x, const BCDSolverOption& o, const BCDSolverParam& par)
{
    // iter number
    if(par.kIter > o.MaxIter)
    {
        return true;
    }
    
    // x change
    Float e_x = (x - par.xprev).norm();
    if(e_x < o.xTol)
    {
        return true;
    }
    
    // obj change
    // ...
    
    // gradient norm, step length, relative error ...
    
    
    
    return false;
}

void bcd_solverInit(Vector& x, BCDSolverParam& par)
{
    int fd = par.m_blockCtrl.full_dim();
    if (fd > 0)
    {
        x.resize(fd);
        x.setZeros();
        par.obj = bcd_objective(x, par);
        
        //
        Float incr = 0.1*numeric_limits<Float>::max();
        par.objprev = par.obj + incr;
        par.xprev = x;
        par.xprev[0] = incr;
        paramInit(x, par);
    }
    else
    {
        cerr << "Invalid problem, please check param setting!" << endl;
    }
    
}

typedef COPT::Solver<Float,Vector,Vector,BCDSolverOption,BCDSolverParam> tagBCDSolver;    
class BCDSolver : public tagBCDSolver
{
public:
    BCDSolver(const BCDSolverParam& bcd_par,
              const BCDSolverOption& bcd_opt)
        :   tagBCDSolver(
                bcd_par,
                bcd_opt,
                bcd_objective,
                bcd_solverInit,
                bcd_oneIteration,
                bcd_terminate)
    {}
};


/////////////////////////////////////////////////////////////////

inline void get_vec_block(const Vector& x0, int id1, int id2, Vector& x)
{
    x.resize(id2-id1);
    for(int i = id1, k = 0; i < id2; ++i, ++k)
        x[k] = x0[i];
}

inline void set_vec_block(Vector& x0, int id1, int id2, const Vector& x)
{
    for(int i = id1, k = 0; i < id2; ++i, ++k)
        x0[i] = x[k];
}

inline Vector get_vec_block(const Vector& x0, int id1, int id2)
{
    Vector x(id2-id1);
    for(int i = id1, k = 0; i < id2; ++i, ++k)
        x[k] = x0[i];
    return x;
}


inline Float shrinkage(Float z, Float tau)
{
    Float t = 1.0 - tau / fabs(z);
    return t > 0.0 ? t * z : 0.0;
}

void updateL(const Vector& x, BCDSolverParam& par)
{
    auto& tpar = par.m_paramAuto;
    int nblock = par.m_blockCtrl.n_block();
    Vector g = par.f_grad(x);
    for(int i = 0; i < nblock; ++i)
    {
        Vector gi = par.m_blockCtrl.get_vec_block(g, i);
        Matrix ggt = gi.mulTrans(gi);
        Float maxe = ggt.frobeniusNorm();
        tpar[i].L = maxe;
    }
}
void paramInit(const Vector& x, BCDSolverParam& par)
{
    int nblock = par.m_blockCtrl.n_block();
    par.m_paramAuto.resize(nblock);
    // init t
    for(int i = 0;i  < nblock; ++i)
        par.m_paramAuto[i].t = 1.0;
    // init L
    updateL(x, par);
}
void paramUpdate(const Vector& x, BCDSolverParam& par)
{
    par.m_paramAutoPrev = par.m_paramAuto;
    const auto& tpar0 = par.m_paramAutoPrev;
    const auto nblock = par.m_blockCtrl.n_block();
    auto& tpar = par.m_paramAuto;
    
    // update t
    for(int i = 0; i < nblock; ++i)
    {
        auto& ti = tpar[i].t;
        ti = 0.5*sqrt( 1.0 + 4.0*ti*ti );
    }
    // update L
    updateL(x, par);
    // update w
    if(par.kIter > 0)
    {
        for(int i = 0; i < nblock; ++i)
        {
            Float wt1 = (tpar0[i].t - 1.0) / tpar0[i].t;
            Float wt2 = sqrt(tpar0[i].L / tpar[i].L);
            tpar[i].w = min(wt1, wt2);
        }
    }
    
}

void solve_subblock(Vector& x, BCDSolverParam& par, int kblock)
{
    const auto& tpar = par.m_paramAuto;
    
    Vector xi = par.m_blockCtrl.get_vec_block(x, kblock);
    if(par.kIter > 1)
    {
        Vector xi0 = par.m_blockCtrl.get_vec_block(par.xprev, kblock);
        xi = xi + tpar[kblock].w*(xi - xi0);
        par.m_blockCtrl.set_vec_block(x, kblock, xi);
    }
    
    // shrinkage for each element
    Vector g = par.f_grad(x);
    Vector gi = par.m_blockCtrl.get_vec_block(g, kblock);
    int dim = gi.size();
    Float tL = 1.0 / tpar[kblock].L;
    Vector tgi = xi - tL * gi;
    for(int i = 0; i < dim; ++i)
    {
        xi[i] = shrinkage(tgi[i], tL);
    }
    par.m_blockCtrl.set_vec_block(x, kblock, xi);
}
Float errorEst(const Vector& x, const BCDSolverParam& par)
{
    return (x - par.xprev).norm();
}


#endif // BCDSOLVER
