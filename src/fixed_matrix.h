#ifndef FIXED_MATRIX_H
#define FIXED_MATRIX_H

#include "utils.h"

class SampleClonal;

typedef enum {
    MATRIX_A = 0,
    MATRIX_B = 1,
    MATRIX_U = 2
} FIXED_MATRIX_t;

struct FixedMatrixEntry
{
    FIXED_MATRIX_t matrix;
    int row;
    int column;
    double value;
    int line;
    std::string rowLabel;
    std::string columnLabel;
};

class FixedMatrix
{
public:
    FixedMatrix();

    static FixedMatrix read(const std::string& filename,
                            const int clones,
                            const std::map<std::string, int>& clusterToIdx,
                            const std::map<std::string, int>& sampleToIdx);

    bool empty() const
    {
        return _entries.empty();
    }

    int clones() const
    {
        return _clones;
    }

    const std::string& source() const
    {
        return _source;
    }

    const std::vector<FixedMatrixEntry>& entries() const
    {
        return _entries;
    }

    void validate(const int cmax,
                  const double mu,
                  const int base,
                  const bool ampdel,
                  const int cnstates,
                  const CNMap& copyNumbers,
                  const DoubleArray& weights,
                  const double maxHomdelSize) const;

    void validate(const int cmax,
                  const double mu,
                  const int base,
                  const bool ampdel,
                  const int cnstates,
                  const CNMap& copyNumbers,
                  const DoubleArray& weights,
                  const SampleClonal& sampleClonal,
                  const DoubleArray& purities,
                  const double maxHomdelSize) const;

    void applyCopyNumbers(IntMatrix& A, IntMatrix& B) const;
    void applyProportions(DoubleMatrix& U) const;
    void conditionProportions(DoubleMatrix& U, const double mu) const;
    void verify(const IntMatrix& A, const IntMatrix& B, const DoubleMatrix& U) const;

private:
    int _clones;
    std::string _source;
    std::vector<FixedMatrixEntry> _entries;
};

#endif // FIXED_MATRIX_H
