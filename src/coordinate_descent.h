#ifndef _COORDINATE_H_
#define _COORDINATE_H_

#include "utils.h"
#include "worker.h"


class CoordinateDescent
{
public:
    CoordinateDescent(const DoubleMatrix& FA,
                      const DoubleMatrix& FB,
                      const DoubleArray& bins,
                      const int m,
                      const int k,
                      const int n,
                      const int cmax,
                      const int d,
                      const double mu,
                      const double maxHomdelSize,
                      const int base,
                      const bool ampdel,
                      const IntSet& ampdelExempt,
                      const CNMap& cn,
                      const FixedMatrix& fixedMatrix,
                      const SampleClonal& sampleClonal,
                      const DoubleArray& purities,
                      const unsigned int iterConvergence,
                      const unsigned int maxIter,
                      const int nrSeeds,
                      const int nrWorkers,
                      const int nrILPthreads,
                      const int timeLimit,
                      const int memoryLimit,
                      const VERBOSITY_t v);

public:
    /// Run the method
    void run();
    /// Get best objective value
    double getObjValue() const
    {
        return _objs[_best];
    }

    IntMatrix getACN() const
    {
        return _CAs[_best];
    }

    IntMatrix getBCN() const
    {
        return _CBs[_best];
    }

    DoubleMatrix getU() const
    {
        return _Us[_best];
    }

private:
    ///
    const DoubleMatrix& _FA;
    ///
    const DoubleMatrix& _FB;
    ///
    const DoubleArray& _bins;
    ///
    const int _m;
    ///
    const int _k;
    /// Number of clones
    const int _n;
    ///
    const int _cmax;
    ///
    const int _d;
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
    /// Sparse user-specified A, B, and U entries
    const FixedMatrix& _fixedMatrix;
    /// Per-sample clonal states
    const SampleClonal& _sampleClonal;
    /// Per-sample purities
    const DoubleArray& _purities;
    /// Number of iterations for checking convergence
    const unsigned int _iterConvergence;
    /// Maximum number of iterations for each seed
    const unsigned int _maxIter;
    /// Number of seeds
    const int _nrSeeds;
    /// Number of workers
    const int _nrWorkers;
    /// Number of LPthreads
    const int _nrILPthreads;
    /// Time limit (seconds)
    const int _timeLimit;
    /// Memory limit (MB)
    const int _memoryLimit;

    /// Dirichlet parameter for generating initial M
    int _sizeBubbles;
    /// All initial seeds
    Double3Matrix _seeds;
    ///
    int _seedidx;
    /// Thread group
    std::vector<std::thread> _threads;
    /// Mutex
    std::mutex _mutexTasks;
    /// Mutex
    std::mutex _mutexResults;
    ///
    DoubleArray _objs;
    ///
    Int3Matrix _CAs;
    ///
    Int3Matrix _CBs;
    ///
    Double3Matrix _Us;
    ///
    int _best;
    ///
    const VERBOSITY_t _v;

    /// Construct random M
    DoubleMatrix buildRandomU();
    /// Construct a sample-clonal-compatible random M
    DoubleMatrix buildSampleClonalU();
    /// Randomly split a sample's residual purity into shares of at least mu each
    DoubleArray splitSamplePurity(const double residual, const int numParts);
    /// Construct random vector summing up to 1
    DoubleArray buildPartitionVector(const int num_parts);
    ///
    void runWorker();
    ///
    void runInstance(const int seedIdx);
};

#endif // _COORDINATE_H_
