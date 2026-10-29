#ifndef ILPSUBSET_H
#define ILPSUBSET_H

#include "utils.h"
#include "fixed_matrix.h"
#include "sample_clonal.h"
#include <gurobi_c++.h>

typedef std::vector<GRBVar> VarArray;
typedef std::vector<VarArray> VarMatrix;
typedef std::vector<VarMatrix> Var3Matrix;
typedef std::vector<Var3Matrix> Var4Matrix;

/// A value that must be enforced by an explicit equality because it conflicts with another fixed value
struct FixedValue
{
    FIXED_MATRIX_t matrix;
    int row;
    int column;
    double value;
};

class IlpSubset
{

public:

    IlpSubset(const int n,
              const int m,
              const int k,
              const int maxc,
              const IntArray& cnstates,
              const double mu,
              const double maxHomdelSize,
              const int base,
              const bool ampdel,
              const IntSet& ampdelExempt,
              const CNMap& cn,
              const DoubleMatrix& FA,
              const DoubleMatrix& FB,
              const DoubleArray& w,
              const FixedMatrix& fixedMatrix,
              const SampleClonal& sampleClonal,
              const DoubleArray& purities,
              const VERBOSITY_t v);


    void fixC(const IntMatrix& barCA, const IntMatrix& barCB);

    void fixU(const DoubleMatrix& barU)
    {
        _barU = barU;
        _fixedMatrix.applyProportions(_barU);
        _mode = MODE_t::CARCH;
        for(int p = 0; p < _k; ++p)
        {
            double sum = 0.0;
            for(int i = 0; i < _n; ++i)
            {
                sum += _barU[i][p];
                assert( (0.0 - TOL <= _barU[i][p] && _barU[i][p] <= 0.0 + TOL) || (_barU[i][p] >= _mu - TOL) || i == 0);
            }
            assert(sum >= 1.0 - TOL && sum <= 1.0 + TOL);
        }
    }

    void init();

    /// Set a MIP start from integer copy numbers and, in FULL mode, optional clone proportions
    void hotStart(const IntMatrix& hCA, const IntMatrix& hCB, const DoubleMatrix& hU = DoubleMatrix());

    int solve(const int timeLimit, const int memoryLimit, const int nrThreads);
    int solve(const int timeLimit, const int memoryLimit, const int nrThreads, const int maxIter);

    void exportModel(const std::string& filename);

    DoubleArray getObjs()
    {
        return _objs;
    }
    Int3Matrix getACNs()
    {
        return _CAs;
    }
    Int3Matrix getBCNs()
    {
        return _CBs;
    }
    Double3Matrix getUs()
    {
        return _Us;
    }

    double runtime()
    {
        return _model.get(GRB_DoubleAttr_Runtime);
    }

    double gap()
    {
        return _model.get(GRB_DoubleAttr_MIPGap);
    }

    double LB()
    {
        return _model.get(GRB_DoubleAttr_LB);
    }

    double UB()
    {
        return _model.get(GRB_DoubleAttr_UB);
    }

    /// Model size: rows, integer and binary variables, and nonzeros
    int numConstrs()
    {
        return _model.get(GRB_IntAttr_NumConstrs);
    }

    int numDiscreteVars()
    {
        return _model.get(GRB_IntAttr_NumIntVars);
    }

    int numNonzeros()
    {
        return _model.get(GRB_IntAttr_NumNZs);
    }

    void printValues();

    static std::pair<IntMatrix, IntMatrix> firstHotStart(const DoubleMatrix &FA, const DoubleMatrix &FB, const int m, const int k, const int n,
                                                         const int cmax, const IntArray& cnstates, const int base, const bool ampdel, const IntSet& ampdelExempt,
                                                         const CNMap& cn);

    /// Maximum number of distinct tumor states of segment s, or 0 when no limit applies: clonal segments are fixed,
    /// and a limit of at least n - 1 (the number of tumor clones) cannot bind
    static int stateLimit(const IntArray& cnstates, const CNMap& cn, const int n, const int s)
    {
        const int limit = segmentCnstates(cnstates, s);
        return cn.count(s) == 0 && limit > 0 && limit < n - 1 ? limit : 0;
    }

private:

    /// Whether the amp-del constraint covers segment s: not clonal and not exempted by the user
    bool ampdelApplies(const int s) const
    {
        return _ampdel && _cn.count(s) == 0 && _ampdelExempt.count(s) == 0;
    }
    int stateLimit(const int s) const
    {
        return stateLimit(_cnstates, _cn, _n, s);
    }
    /// Whether any segment has a limit on its distinct tumor states
    bool anyStateLimit() const
    {
        for(int s = 0; s < _m; ++s)
        {
            if(stateLimit(s) > 0)
            {
                return true;
            }
        }
        return false;
    }
    /// Whether the model needs the binary expansion of the copy numbers: the vA/vB products and buildNext of FULL
    /// mode, and the pairwise state limit
    bool needsBits() const
    {
        return _mode == MODE_t::FULL ||
               (_mode == MODE_t::CARCH && g_cnstatesFormulation == CNSTATES_PAIRWISE && anyStateLimit());
    }
    /// Build main common variables
    void buildVariables();
    /// Build the variables that depend on the chosen second objective
    void buildOptionalVariables();
    /// Build main common constraints
    void buildConstraints();
    /// Build the constraints that depend on the chosen second objective
    void buildOptionalConstraints();
    /// Limit segment s to `limit` tumor states through slot states numbered by first use
    void buildSlotConstraints(const int s, const int limit);
    /// Limit segment s to `limit` tumor states through pairwise equalities of copy-number bits (benchmarks only)
    void buildPairwiseConstraints(const int s, const int limit);
    /// Start the slot variables of segment s from tumor copy numbers, numbering slots by first use
    void startSlots(const int s, const IntArray& hCA, const IntArray& hCB);
    /// Derive variable bounds from the normal clone, clonal copy numbers, fixed entries, and sample-clonal states
    void computeFixedBounds();
    /// Add explicit equalities for fixed values that conflict with another fixed value
    void addFixedValueConflicts();
    /// Enforce per-sample clonal states on present clones in FULL mode
    void fixSampleClonal();
    /// Build objective
    void buildObjective();
    /// Add the constraint to avoid current solution and to find the next one
    void buildNext();
    /// Set specific values of variables for testing
    void test();
    /// Extract from current solution the values of the variables corresponding to copy numbers of clones
    IntMatrix getCA();
    /// Extract from current solution the values of the variables corresponding to B-allele copy numbers of clones
    IntMatrix getCB();
    /// Extract from current solution the values of the variables corresponding to the proportions
    DoubleMatrix getU();
    ///
    double symmCoeff(const int index);
    /// Upper bound on the total copy number of segment s
    int upperCN(const int s) const;

    /// Number of clones
    const int _n;
    /// Number of total segments
    const int _m;
    /// Number of samples
    const int _k;
    /// Maximum copy number
    const int _cmax;
    /// Maximum number of distinct tumor copy-number states of each segment (-1 for no limit)
    const IntArray& _cnstates;
    ///
    const double _mu;
    /// Maximum cluster-weight fraction that can have a homozygous deletion
    const double _maxHomdelSize;
    ///
    const int _base;
    ///
    const bool _ampdel;
    /// Clusters exempt from the amp-del constraint
    const IntSet& _ampdelExempt;
    ///
    const CNMap& _cn;
    /// Read depth ratio r_{t, p} for each segment t in each sample p
    const DoubleMatrix& _FA;
    /// B-allele frequency d_{t, p} for each segment t in each sample p
    const DoubleMatrix& _FB;
    ///
    const DoubleArray& _w;
    /// Sparse user-specified A, B, and U entries
    const FixedMatrix& _fixedMatrix;
    /// Per-sample clonal states
    const SampleClonal& _sampleClonal;
    /// Per-sample purities derived from sample-clonal states
    const DoubleArray& _purities;
    ///
    const VERBOSITY_t _v;
    ///
    IntMatrix _barCA;
    ///
    IntMatrix _barCB;
    ///
    DoubleMatrix _barU;
    ///
    MODE_t _mode;
    /// Maximum B-allele copy number
    const int _M;

    /// Bounds on cA[s][i] and cB[s][i] implied by fixed values
    IntMatrix _lbA;
    IntMatrix _ubA;
    IntMatrix _lbB;
    IntMatrix _ubB;
    /// Bounds on u[i][p] implied by fixed values
    DoubleMatrix _lbU;
    DoubleMatrix _ubU;
    /// Fixed values that could not be expressed as bounds
    std::vector<FixedValue> _conflicts;

    /// Gurobi environment
    GRBEnv _env;
    /// Gurobi model
    GRBModel _model;
    /// Variable yA[s][p] is the error between observed and estimated fractional A-allele copy number of segment s in sample p
    VarMatrix _yA;
    /// Variable yB[s][p] is the error between observed and estimated fractional B-allele copy number of segment s in sample p
    VarMatrix _yB;
    /// Variable fA[t][p] is the fractional A-allele copy number of chosen segment t in sample p
    VarMatrix _fA;
    /// Variable fB[t][p] is the fractional B-allele copy number of chosen segment t in sample p
    VarMatrix _fB;
    /// Variable cA[t][i] is the integer A-allele copy number of the chosen segment t in clone i
    VarMatrix _cA;
    /// Variable cB[t][i] is the integer B-allele copy number of the chosen segment t in clone i
    VarMatrix _cB;
    /// Variable cA[b][t][i] is the b-th bit in the binary representation of cA[t][i]
    Var3Matrix _bitcA;
    /// Variable cB[b][t][i] is the b-th bit in the binary representation of cB[t][i]
    Var3Matrix _bitcB;
    /// Variable u[i][p] is the proportion of clone i in sample p
    VarMatrix _u;
    ///
    VarMatrix _x;
    ///
    VarArray _adA;
    ///
    VarArray _adB;
    /// Variable z[b][t][i][p] is the product between c[b][t][i] and u[i][p]
    Var4Matrix _vA;
    /// Variable zB[b][t][i][p] is the product between cB[b][t][i] and u[i][p]
    Var4Matrix _vB;
    /// Variable z[t][i][l] = 1 iff tumor clone i of segment t takes the state of slot l
    Var3Matrix _z;
    /// Variable sA[t][l] is the A-allele copy number of slot l of segment t (slot formulation)
    VarMatrix _sA;
    /// Variable sB[t][l] is the B-allele copy number of slot l of segment t (slot formulation)
    VarMatrix _sB;
    /// Variable barcA[b][t][i][j] = 1 iff _cA[b][t][i] and _cA[b][t][j] are different
    Var4Matrix _barcIJA;
    /// Variable barcB[b][t][i][j] = 1 iff _cB[b][t][i] and _cB[b][t][j] are different
    Var4Matrix _barcIJB;

    /// Objective value for each found solution
    DoubleArray _objs;
    /// Integer A-allele copy numbers of the clones for each solution
    Int3Matrix _CAs;
    /// Integer B-allele copy numbers of the clones for each solution
    Int3Matrix _CBs;
    /// Proportions of the clones in the samples for each solution
    Double3Matrix _Us;

};

#endif // ILPSUBSET_H
