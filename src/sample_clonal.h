#ifndef _SAMPLE_CLONAL_H_
#define _SAMPLE_CLONAL_H_

#include "utils.h"


struct SampleClonalEntry
{
    int sample;
    int cluster;
    int cnA;
    int cnB;
    int line;
    std::string sampleLabel;
    std::string clusterLabel;

    int total() const
    {
        return cnA + cnB;
    }
};


struct SampleClonalCalibration
{
    DoubleArray gamma;
    DoubleArray purity;
};


class SampleClonal
{
public:
    SampleClonal();

    static SampleClonal read(const std::string& filename,
                             const std::map<std::string, int>& clusterToIdx,
                             const std::map<std::string, int>& sampleToIdx);

    bool empty() const
    {
        return _entries.empty();
    }

    const std::vector<SampleClonalEntry>& entries() const
    {
        return _entries;
    }

    const std::string& source() const
    {
        return _source;
    }

    std::vector<SampleClonalEntry> entriesForSample(const int sample) const;
    SampleClonalCalibration calibrate(const DoubleMatrix& R, const int samples) const;
    int maxTotalCopyNumber() const;
    bool samplesCompatible(const int left, const int right) const;
    void validate(const int cmax,
                  const int base,
                  const bool ampdel,
                  const IntSet& ampdelExempt,
                  const IntArray& cnstates,
                  const DoubleArray& weights,
                  const DoubleArray& purities,
                  const double maxHomdelSize) const;
    void verify(const IntMatrix& A,
                const IntMatrix& B,
                const DoubleMatrix& U,
                const DoubleArray& purities) const;

private:
    std::string _source;
    std::vector<SampleClonalEntry> _entries;
};


#endif
