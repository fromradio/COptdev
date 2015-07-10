#ifndef BCDSOLVER


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
using namespace std;

typedef kernel::podscalar Float;



/// function F: the differentiable part
/// need to be derived to overload the eval() and grad() method
class FuncF
{
public:
    
    virtual Float eval(const Vector& x) = 0;
    virtual Vector grad(const Vector& x) = 0;
    inline Float operator()(const Vector& x) {return eval(x);}
};

/// function R: the non-smooth part
/// usually in the form of 1-norm or a indicator function
class FuncR
{
public: // Function type
    
    enum RType { Norm1, Indicator };
    
public: // interface
    
    FuncR(RType type = Norm1, int dim = 1)
        :   m_type(type)
        ,   m_dim(dim)
    {
        assert(type == Norm1);
        m_funcEval = &FuncR::func_norm1;
    }

    FuncR(RType type, int dim,
          const vector<Float>& lb, const vector<Float>& ub)
        :   m_type(type)
        ,   m_dim(dim)
        ,   m_lb(lb)
        ,   m_ub(ub)
    {
        assert(type == Indicator);
        m_funcEval = &FuncR::func_indicator_rectangel;
    }
    
    Float eval(const Vector& x) {return (this->*m_funcEval)(x);}
    inline Float operator()(const Vector& x) {return eval(x); }
    
    int dim()const{return m_dim;}
    RType type()const{return m_type;}
    
    // getter & setter
    void set_type(RType t) { m_type = t; }
    void set_dim(int dim) { m_dim = dim; }
    void set_bound(const vector<Float>& lb, const vector<Float>& ub)
    { m_lb = lb; m_ub = ub; }
    
    
private: // typedef
    
    typedef Float (FuncR::*FUNC_EVAL)(const Vector& x);
    
private: // property
    
    RType m_type;
    int m_dim;
    

    FUNC_EVAL m_funcEval;
    vector<Float> m_lb, m_ub; // lower bound, upper bound
    
    
private: // inner implement
    
    Float func_norm1(const Vector& x) { return x.absNorm(); }
    
    Float func_indicator_rectangel(const Vector& x)
    {
        assert(m_lb.size() == x.size() && m_ub.size() == x.size());
        auto n = x.size();
        bool inside = true;
        for(int i = 0; i < n; ++i)
        {
            if(x[i] < m_lb[i] || x[i] > m_ub[i])
            {
                inside = false;
                break;
            }
        }
        if(inside)
            return 0;
        else
            return numeric_limits<Float>::infinity();
    }
    
};


/// solver parameters
typedef struct BCDSolverParam
{
    // objective
    Float tau;
    
    // iteration
    int maxIter;       // max iteration number
    Float errorTol;    // 
    Float xTol;
    Float objTol;
    int verbose;
    

    // defaults
    BCDSolverParam()
        :   maxIter(100)
        ,   errorTol(1e-6)
        ,   xTol(1e-6)
        ,   objTol(1e-6)
        ,   verbose(0)
        ,   tau(1.0)
    {}
    
}BCDSolverParam;


class BCDSolver : public COPT::GeneralSolver<kernel>
{
public:
    
    BCDSolver(
            BCDSolverParam param,
            FuncF& f,
            vector<FuncR>& r_vec);
    
public:
    
     void   init_x(const Vector& x);
     Float  objective(const Vector& x)const;
     
     /** get the result */
     const kernel::Vector& result() const;
     
private: // 
     
    ////
    /** vritual function that has to be implemented by derived classes */
	/** something happens in the beginning of solving */
    void solvingBegin();
   
	/** the real solving part */
//    void doSolve();
    
	/** the real computation part */
    void doCompute();
    
	/** whether the terminal satisfied */
    bool terminalSatisfied() const;
	
	/** an iteration 
	 *  the estimated error is returned for iteration
	 */
    podscalar doOneIteration();
    
    /** compute the objective function */
    podscalar objective() const;
   
    
    
private: // basic param
    
    BCDSolverParam m_param;
    FuncF& m_f;
    vector<FuncR>& m_rVec;
    
private: // opt vars
    
    Float m_obj, m_objprev;    
    Float m_step;
    Vector m_x, m_xprev, m_xprev2;
    Vector m_gradf;
    Vector m_descent, m_descentprev;
    
    vector<int> m_block;
    int m_dim;
    int m_nblock;
    
private: // iteration ctrl
    
    int m_kIter;
    bool m_isInit;
    
private: // inner tools
    
    inline Float errorEst();
    
    void solve_subblock(int kblock);
    
    // auto updated parameters
    typedef struct paramAuto
    {
        Float t;
        Float w;
        Float L;
    }paramAuto;
    vector<paramAuto> m_paramAuto, m_paramAutoPrev;
    void paramInit();
    void paramUpdate();
    void updateL();
    vector<Vector> m_xhat;
    
    // block control
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
        }
        
        int n_block()const {return m_nblock;}
        int block_dim(int i)const{return m_block[i];}
        void block_id_interval(int i, int& a, int& b)const
        {
            a = m_blockIdx[i];
            b = m_blockIdx[i+1];
        }
        
    private:
        vector<int> m_block;
        vector<int> m_blockIdx;
        int m_nblock;
    };
    BlockCtrl m_blockCtrl;
    
    static void get_block(const Vector& x0, int id1, int id2, Vector& x)
    {
        x.resize(id2-id1);
        for(int i = id1, k = 0; i < id2; ++i, ++k)
            x[k] = x0[i];
    }
    static void set_block(Vector& x0, int id1, int id2, const Vector& x)
    {
        for(int i = id1, k = 0; i < id2; ++i, ++k)
            x0[i] = x[k];
    }
    
    Float shrinkage(Float z, Float tau)
    {
        Float t = 1.0 - tau / fabs(z);
        return t > 0.0 ? t * z : 0.0;
    }
    
    void print();
    
};

///////////////////////////////////////////////////////////////////////////////////////////
/*
 * 
 *              Implementation
 * 
 * 
 * 
 *  
 */

BCDSolver::BCDSolver(BCDSolverParam param, FuncF &f, vector<FuncR> &r_vec)
    :   COPT::GeneralSolver<kernel>(param.maxIter,param.xTol)
    ,   m_param(param)
    ,   m_f(f)
    ,   m_rVec(r_vec)
    ,   m_nblock(0)
    ,   m_dim(0)
    ,   m_kIter(0)
    ,   m_isInit(false)
    
{
    m_nblock = r_vec.size();
    m_block.resize(m_nblock);
    m_paramAuto.resize(m_nblock);
    m_dim = 0;
    for(int i = 0; i < m_nblock; ++i)
    {
        m_block[i] = r_vec[i].dim();
        m_dim += m_block[i];
    }
    m_blockCtrl.setBlocking(m_block);
    m_xhat.resize(m_nblock);
    
    // resize
    m_x.resize(m_dim);
    m_xprev.resize(m_dim);
    m_xprev2.resize(m_dim);
    m_gradf.resize(m_dim);
    m_descent.resize(m_dim);
    m_descentprev.resize(m_dim);

    // init value
    m_objprev = numeric_limits<Float>::max();
    Vector tx(m_dim);
    for(auto& t:tx) t = 100.0 * m_param.xTol;
    m_xprev  = m_x      + tx;
    m_xprev2 = m_xprev  + tx;
    
    paramInit();
}

Float BCDSolver::objective(const Vector &x) const
{
    Float obj = 0;
    if (m_nblock == 1)
        obj = m_rVec[0](x);
    else
    {
        for(int i = 0; i < m_nblock; ++i)
        {
            Vector xi;
            int ia, ib;
            m_blockCtrl.block_id_interval(i, ia, ib);
            get_block(x, ia, ib, xi);
            obj += m_rVec[i](xi);
        }  
    }

    obj *= m_param.tau;
    obj += m_f(x);
    return obj;
}

void BCDSolver::init_x(const Vector& x)
{
    if (x.size() == m_dim)
    {
        m_x = x;
        m_obj = objective(m_x);
        m_isInit = true;
    }    
    else
    {
        cout << "error in init: vector dimension inconsistent" << endl;
        return;
    }
        
    if (m_param.verbose)
        print();
}



void BCDSolver::paramInit()
{
    auto& par = m_paramAuto;
    // t
    for(int k = 0; k < m_nblock; ++k)
    {
        par[k].t = 1.0;
    }
    
    // L 
    updateL();
    
}

void BCDSolver::paramUpdate()
{
    auto& par  = m_paramAuto;
    auto& par0 = m_paramAutoPrev;
    par0 = par;
    
    // t
    for(int i = 0; i < m_nblock; ++i)
    {
        par[i].t = 0.5*sqrt(1.0+4.0*par[i].t*par[i].t);
    }   
    
    // L 
    updateL();
    
    // w
    if(m_kIter > 0)
    {
        for(int i = 0; i < m_nblock; ++i)
        {
            Float wt1 = (par0[i].t - 1.0) / par0[i].t;
            Float wt2 = sqrt(par0[i].L / par[i].L);
            par[i].w = min(wt1, wt2);
        }
    }
}

void BCDSolver::updateL()
{
    Vector g = m_f.grad(m_x);
    int id1, id2;
    auto& par = m_paramAuto;
    for(int kb = 0; kb < m_nblock; ++kb)
    {
        m_blockCtrl.block_id_interval(kb, id1, id2);
        Vector gk;
        get_block(g, id1, id2, gk);
        Matrix ggt = gk.mulTrans(gk);
        Float maxe = ggt.frobeniusNorm();
        par[kb].L = maxe;
    }
    
}


void BCDSolver::solve_subblock(int kblock)
{
    const auto& ib = kblock;
    const auto& par = m_paramAuto;
    int id1, id2;
    m_blockCtrl.block_id_interval(ib,id1,id2);
    Vector xi;
    get_block(m_x, id1, id2, xi);
    if(m_kIter > 1)
    {
        Vector xi0;
        get_block(m_xprev, id1, id2, xi0);
        xi = xi + par[ib].w*(xi - xi0);
        set_block(m_x, id1, id2, xi);
        m_xhat[kblock] = xi;
    }
    
    
    // shrinkage for each element
    Vector grad = m_f.grad(m_x);
    Vector gi;
    get_block(grad, id1, id2, gi);
    
    int n = id2 - id1;
    Float tL = 1.0 / par[ib].L;
    Vector tg = xi - tL * gi;
    for(int i = 0; i < n; ++i)
    {
        xi[i] = shrinkage(tg[i], tL);
    }
    set_block(m_x, id1, id2, xi);
    
}


void BCDSolver::solvingBegin()
{
    if(!m_isInit)
    {
        m_x.setZeros();
        m_obj = objective(m_x);
    }
    
    // init xhat (x hat)
    for(int i = 0; i < m_nblock; ++i)
    {
        int ia, ib;
        m_blockCtrl.block_id_interval(i, ia, ib);
        Vector xi;
        get_block(m_x, ia, ib, xi);
        m_xhat[i] = xi;
    }
}

void BCDSolver::doCompute()
{
    
}

kernel::podscalar BCDSolver::doOneIteration()
{
    // preproc
    m_xprev  = m_x;
    m_objprev = m_obj;

    // solve 
    for(int kblock = 0; kblock < m_nblock; ++kblock)
    {
        solve_subblock(kblock);
    }
    ++m_kIter;
    paramUpdate();
    
    // postproc
    m_obj = objective(m_x);
    
    
    if (m_param.verbose)
        print();
    
    return errorEst();
}

Float BCDSolver::errorEst()
{
    return (m_x - m_xprev).norm();
}

bool BCDSolver::terminalSatisfied()const
{
    // iter number
    if (__iter_num > __max_iteration)
        return true;
    
    // x change
    Float e_x = (m_x - m_xprev).norm();
    if (e_x > m_param.xTol)
        return true;
    
    // obj change
    Float e_obj =  fabs(m_obj - m_objprev);
    if (e_obj > m_param.objTol)
        return true;
    
    // gradient norm
    // ...
    
    // step length
    // ...
    
    // relative error
    // ...
}

kernel::podscalar BCDSolver::objective()const
{
    return objective(m_x);
}

const kernel::Vector& BCDSolver::result() const
{
    return m_x;
}

void BCDSolver::print()
{
    cout << "Iter: " << m_kIter << " - x: " << m_x << " - obj: " << m_obj << endl;
}

#define BCDSOLVER
#endif // BCDSOLVER
