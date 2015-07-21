
#ifndef BCDSOLVERIMPL
#define BCDSOLVERIMPL


/**
  Implementation of BCD solver
  
*/


#include "BCD.hpp"



/// optimization module functions
Float   bcdObjective(const Vector& x, const BCDSolverParam& par);
Float   bcdOneIteration(Vector& x, BCDSolverParam& par);
bool    bcdTerminate(const Vector& x, const BCDSolverOption& o, const BCDSolverParam& par);
void    bcdSolverInit(Vector& x, BCDSolverParam& par);

/// tool functions
// vector slicing
inline void getVectorBlock(const Vector& x0, int id1, int id2, Vector& x);
inline void setVectorBlock(Vector& x0, int id1, int id2, const Vector& x);
inline Vector getVectorBlock(const Vector& x0, int id1, int id2);
// optimization related
inline Float shrinkage(Float z, Float tau);
void updateL(const Vector& x, BCDSolverParam& par);
void paramInit(const Vector& x, BCDSolverParam& par);
void paramUpdate(const Vector& x, BCDSolverParam& par);
void solveSubblock(Vector& x, BCDSolverParam& par, int kblock);
Float errorEst(const Vector& x, const BCDSolverParam& par);




//////////////////////////////////////////////////////////////////////


/// bcd solver
BCDSolver::BCDSolver(const BCDSolverParam& bcd_par,
          const BCDSolverOption& bcd_opt)
    :   tagBCDSolver(
            bcd_par,
            bcd_opt,
            bcdObjective,
            bcdSolverInit,
            bcdOneIteration,
            bcdTerminate)
{}


///
/// optimization modules
/// 
Float bcdObjective(const Vector& x, const BCDSolverParam& par)
{
    Float obj = 0;
    // r part
    for(int i = 0; i < par.__block_ctrl.nBlock(); ++i)
        obj += par.rEval(x, i);
    // f part
    obj += par.fEval(x);
    return obj;
}

Float bcdOneIteration(Vector& x, BCDSolverParam& par)
{
    par.__x_prev = x;
    par.__obj_prev = par.__obj;
    // update x block-wise
    auto nblock = par.__block_ctrl.nBlock();
    for(int kb = 0; kb < nblock; ++kb)
    {
        solveSubblock(x, par, kb);
    }
    
    //
    ++par.__k_iter;
    paramUpdate(x, par);
    par.__obj = bcdObjective(x, par);
    
    return errorEst(x, par);
}

bool bcdTerminate(const Vector& x, const BCDSolverOption& o, const BCDSolverParam& par)
{
    // iter number
    if(par.__k_iter > o.MaxIter)
    {
        return true;
    }
    
    // x change
    Float e_x = (x - par.__x_prev).norm();
    if(e_x < o.XTol)
    {
        return true;
    }
    
    // obj change
    // ...
    
    // gradient norm, step length, relative error ...
    
    
    
    return false;
}

void bcdSolverInit(Vector& x, BCDSolverParam& par)
{
    int fd = par.__block_ctrl.fullDimension();
    if (fd > 0)
    {
        x.resize(fd);
        x.setZeros();
        par.__obj = bcdObjective(x, par);
        
        //
        Float incr = 0.1*numeric_limits<Float>::max();
        par.__obj_prev = par.__obj + incr;
        par.__x_prev = x;
        par.__x_prev[0] = incr;
        paramInit(x, par);
    }
    else
    {
        cerr << "Invalid problem, please check param setting!" << endl;
    }
    
}


///
/// BlockCtrl
///
void BlockCtrl::setBlocking(const vector<int> &blocking)
{
    __block = blocking;
    __nblock = __block.size();
    __blockIdx.resize(__nblock+1);
    __blockIdx[0] = 0;
    for(int i = 0 ; i < __nblock; ++i)
    {
        __blockIdx[i+1] = __block[i];
        __blockIdx[i+1] += __blockIdx[i];
    }
    __is_set = true;
}

int BlockCtrl::nBlock()const 
{
    return __nblock;
}

int BlockCtrl::blockDimension(int i)const
{
    return __block[i];
}

int BlockCtrl::fullDimension()const
{
    int dim = 0;
    for(const auto& t:__block) dim += t;
    return dim;
}

void BlockCtrl::blockIdexInterval(int i, int& a, int& b)const
{
    a = __blockIdx[i];
    b = __blockIdx[i+1];
}

Vector BlockCtrl::getVectorBlock(const Vector& x0, int kblock)const
{
    assert(__is_set);
    int ia, ib;
    blockIdexInterval(kblock, ia, ib);
    return ::getVectorBlock(x0, ia, ib);
}

void BlockCtrl::setVectorBlock(Vector& x0, int kblock, const Vector& x)const
{
    assert(__is_set);
    int ia, ib;
    blockIdexInterval(kblock, ia, ib);
    ::setVectorBlock(x0, ia, ib, x);
}


///
/// r info
/// 
FuncRInfo::FuncRInfo(RType tp, int dm)
    :   type(tp)
    ,   dim(dm)
{}


///
/// solver parameters
/// 
BCDSolverParam::BCDSolverParam()
    :   BasicParameter()
    ,   __verbose(0)
    ,   __tau(1.0)
    ,   __k_iter(0)
{}

BCDSolverParam::BCDSolverParam(
        const function<Float(const Vector&)>& fEval,
        const function<Vector(const Vector&)>& fGrad,
        const vector<FuncRInfo>& ri,
        Float tau
        )
    :   BasicParameter()
    ,   __verbose(0)
    ,   __tau(1.0)
    ,   __k_iter(0)
{
    setF(fEval, fGrad);
    setR(ri);
    setTau(tau);
}

Float BCDSolverParam::riEval(const FuncRInfo& r, const Vector& x)const
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

void BCDSolverParam::setF(
        const function<Float(const Vector&)>& fEval,
        const function<Vector(const Vector&)>& fGrad)
{
    this->__fEval = fEval;
    this->__fGrad = fGrad;
}

void BCDSolverParam::setR(const vector<FuncRInfo>& ri)
{
    this->__ri = ri;
    vector<int> blocking(ri.size());
    for(int i = 0; i < ri.size(); ++i)
    {
        blocking[i] = ri[i].dim;
    }
    __block_ctrl.setBlocking(blocking);
}

void BCDSolverParam::setTau(Float tau)
{
    this->__tau = tau;
}

inline Float    BCDSolverParam::fEval(const Vector& x)const
{
    return __fEval(x);
}

inline Vector   BCDSolverParam::fGrad(const Vector& x)const
{
    return __fGrad(x);
}

inline Float    BCDSolverParam::rEval(const Vector& x, int kblock)const
{
    return riEval(__ri[kblock], __block_ctrl.getVectorBlock(x,kblock));
}


///
/// tool functions
/// 
inline void getVectorBlock(const Vector& x0, int id1, int id2, Vector& x)
{
    x.resize(id2-id1);
    for(int i = id1, k = 0; i < id2; ++i, ++k)
        x[k] = x0[i];
}

inline void setVectorBlock(Vector& x0, int id1, int id2, const Vector& x)
{
    for(int i = id1, k = 0; i < id2; ++i, ++k)
        x0[i] = x[k];
}

inline Vector getVectorBlock(const Vector& x0, int id1, int id2)
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
    auto& tpar = par.__param_auto;
    int nblock = par.__block_ctrl.nBlock();
    Vector g = par.fGrad(x);
    for(int i = 0; i < nblock; ++i)
    {
        Vector gi = par.__block_ctrl.getVectorBlock(g, i);
        Matrix ggt = gi.mulTrans(gi);
        Float maxe = ggt.frobeniusNorm();
        tpar[i].L = maxe;
    }
}
void paramInit(const Vector& x, BCDSolverParam& par)
{
    int nblock = par.__block_ctrl.nBlock();
    par.__param_auto.resize(nblock);
    // init t
    for(int i = 0;i  < nblock; ++i)
        par.__param_auto[i].t = 1.0;
    // init L
    updateL(x, par);
}
void paramUpdate(const Vector& x, BCDSolverParam& par)
{
    par.__param_auto_prev = par.__param_auto;
    const auto& tpar0 = par.__param_auto_prev;
    const auto nblock = par.__block_ctrl.nBlock();
    auto& tpar = par.__param_auto;
    
    // update t
    for(int i = 0; i < nblock; ++i)
    {
        auto& ti = tpar[i].t;
        ti = 0.5*sqrt( 1.0 + 4.0*ti*ti );
    }
    // update L
    updateL(x, par);
    // update w
    if(par.__k_iter > 0)
    {
        for(int i = 0; i < nblock; ++i)
        {
            Float wt1 = (tpar0[i].t - 1.0) / tpar0[i].t;
            Float wt2 = sqrt(tpar0[i].L / tpar[i].L);
            tpar[i].w = min(wt1, wt2);
        }
    }
    
}

void solveSubblock(Vector& x, BCDSolverParam& par, int kblock)
{
    const auto& tpar = par.__param_auto;
    
    Vector xi = par.__block_ctrl.getVectorBlock(x, kblock);
    if(par.__k_iter > 1)
    {
        Vector xi0 = par.__block_ctrl.getVectorBlock(par.__x_prev, kblock);
        xi = xi + tpar[kblock].w*(xi - xi0);
        par.__block_ctrl.setVectorBlock(x, kblock, xi);
    }
    
    // shrinkage for each element
    Vector g = par.fGrad(x);
    Vector gi = par.__block_ctrl.getVectorBlock(g, kblock);
    int dim = gi.size();
    Float tL = 1.0 / tpar[kblock].L;
    Vector tgi = xi - tL * gi;
    for(int i = 0; i < dim; ++i)
    {
        xi[i] = shrinkage(tgi[i], tL);
    }
    par.__block_ctrl.setVectorBlock(x, kblock, xi);
}
Float errorEst(const Vector& x, const BCDSolverParam& par)
{
    return (x - par.__x_prev).norm();
}




#endif
