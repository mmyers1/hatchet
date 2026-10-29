#include "sample_clonal.h"

#include <cerrno>

namespace
{
const std::string FORMAT_MARKER = "# HATCHET_SAMPLE_CLONAL 1";
const std::string HEADER = "sample\tcluster\tcn_a\tcn_b";
const double PURITY_TOLERANCE = 0.05;

std::string formatError(const std::string& source, const int line, const std::string& message)
{
    return source + ":" + std::to_string(line) + ": " + message;
}

std::vector<std::string> splitTabs(const std::string& line)
{
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while(std::getline(stream, field, '\t'))
    {
        fields.push_back(field);
    }
    if(!line.empty() && line.back() == '\t')
    {
        fields.push_back("");
    }
    return fields;
}

int parseCopyNumber(const std::string& rawValue, const std::string& source, const int line)
{
    const std::string value = trim(rawValue);
    errno = 0;
    char* end = NULL;
    const long parsed = std::strtol(value.c_str(), &end, 10);
    if(value.empty() || errno != 0 || end == value.c_str() || *end != '\0' || parsed < 0 ||
       parsed > std::numeric_limits<int>::max())
    {
        throw formatError(source, line, "copy-number value '" + value + "' must be a nonnegative integer");
    }
    return static_cast<int>(parsed);
}
}


SampleClonal::SampleClonal()
    : _source()
    , _entries()
{
}


SampleClonal SampleClonal::read(const std::string& filename,
                                const std::map<std::string, int>& clusterToIdx,
                                const std::map<std::string, int>& sampleToIdx)
{
    std::ifstream input(filename.c_str());
    if(!input.good())
    {
        throw "The specified sample-clonal file \"" + filename + "\" does not exist or cannot be opened!";
    }

    SampleClonal result;
    result._source = filename;
    std::string line;
    if(!std::getline(input, line))
    {
        throw formatError(filename, 1, "expected format marker '" + FORMAT_MARKER + "'");
    }
    if(line.size() >= 3 && static_cast<unsigned char>(line[0]) == 0xEF &&
       static_cast<unsigned char>(line[1]) == 0xBB && static_cast<unsigned char>(line[2]) == 0xBF)
    {
        line.erase(0, 3);
    }
    if(!line.empty() && line.back() == '\r')
    {
        line.pop_back();
    }
    if(line != FORMAT_MARKER)
    {
        throw formatError(filename, 1, "expected format marker '" + FORMAT_MARKER + "'");
    }

    bool headerSeen = false;
    std::set<std::pair<int, int> > coordinates;
    std::map<int, std::vector<SampleClonalEntry> > bySample;
    int lineNumber = 1;
    while(std::getline(input, line))
    {
        ++lineNumber;
        if(!line.empty() && line.back() == '\r')
        {
            line.pop_back();
        }
        if(trim(line).empty())
        {
            continue;
        }
        if(line[0] == '#')
        {
            if(line.find("# HATCHET_SAMPLE_CLONAL") == 0)
            {
                throw formatError(filename, lineNumber, "the format marker may only appear on line 1");
            }
            continue;
        }
        if(!headerSeen)
        {
            if(line != HEADER)
            {
                throw formatError(filename, lineNumber, "expected tab-separated header '" + HEADER + "'");
            }
            headerSeen = true;
            continue;
        }
        if(line == HEADER)
        {
            throw formatError(filename, lineNumber, "duplicate header");
        }

        const std::vector<std::string> fields = splitTabs(line);
        if(fields.size() != 4)
        {
            throw formatError(filename, lineNumber, "expected 4 tab-separated fields");
        }
        if(fields[0].empty())
        {
            throw formatError(filename, lineNumber, "sample identifier must not be empty");
        }
        if(fields[1].empty())
        {
            throw formatError(filename, lineNumber, "cluster identifier must not be empty");
        }

        const std::map<std::string, int>::const_iterator sample = sampleToIdx.find(fields[0]);
        if(sample == sampleToIdx.end())
        {
            throw formatError(filename, lineNumber, "unknown SEG sample identifier '" + fields[0] + "'");
        }
        const std::map<std::string, int>::const_iterator cluster = clusterToIdx.find(fields[1]);
        if(cluster == clusterToIdx.end())
        {
            throw formatError(filename, lineNumber, "unknown SEG cluster identifier '" + fields[1] + "'");
        }
        const std::pair<int, int> coordinate(sample->second, cluster->second);
        if(coordinates.count(coordinate) > 0)
        {
            throw formatError(filename, lineNumber, "duplicate sample-clonal entry");
        }
        coordinates.insert(coordinate);

        SampleClonalEntry entry;
        entry.sample = sample->second;
        entry.cluster = cluster->second;
        entry.cnA = parseCopyNumber(fields[2], filename, lineNumber);
        entry.cnB = parseCopyNumber(fields[3], filename, lineNumber);
        entry.line = lineNumber;
        entry.sampleLabel = fields[0];
        entry.clusterLabel = fields[1];
        result._entries.push_back(entry);
        bySample[entry.sample].push_back(entry);
    }

    if(!headerSeen)
    {
        throw formatError(filename, lineNumber, "missing sample-clonal header");
    }
    if(result._entries.empty())
    {
        throw formatError(filename, lineNumber, "the file must contain at least one sample-clonal entry");
    }
    for(std::map<std::string, int>::const_iterator sample = sampleToIdx.begin(); sample != sampleToIdx.end(); ++sample)
    {
        const std::vector<SampleClonalEntry>& entries = bySample[sample->second];
        if(entries.size() < 2)
        {
            throw formatError(filename, lineNumber,
                              "sample '" + sample->first + "' must have at least two sample-clonal entries");
        }
        // only the first two entries calibrate the sample; any further entries are only enforced
        if(entries[0].total() == entries[1].total())
        {
            throw formatError(filename, entries[1].line,
                              "sample '" + sample->first + "' must have two distinct total copy-number states");
        }
    }
    return result;
}


std::vector<SampleClonalEntry> SampleClonal::entriesForSample(const int sample) const
{
    std::vector<SampleClonalEntry> result;
    for(std::vector<SampleClonalEntry>::const_iterator entry = _entries.begin(); entry != _entries.end(); ++entry)
    {
        if(entry->sample == sample)
        {
            result.push_back(*entry);
        }
    }
    return result;
}


SampleClonalCalibration SampleClonal::calibrate(const DoubleMatrix& R, const int samples) const
{
    SampleClonalCalibration result;
    result.gamma = DoubleArray(samples, 0.0);
    result.purity = DoubleArray(samples, 0.0);
    for(int sample = 0; sample < samples; ++sample)
    {
        // entries after the first two are enforced during solving but do not calibrate the sample
        const std::vector<SampleClonalEntry> entries = entriesForSample(sample);
        if(entries.size() < 2)
        {
            throw _source + ": missing sample-clonal entries for sample index " + std::to_string(sample);
        }
        const double firstRdr = R[entries[0].cluster][sample];
        const double secondRdr = R[entries[1].cluster][sample];
        const int firstTotal = entries[0].total();
        const int secondTotal = entries[1].total();
        const double denominator = firstRdr * (secondTotal - 2) + secondRdr * (2 - firstTotal);
        if(!std::isfinite(firstRdr) || !std::isfinite(secondRdr) || std::fabs(denominator) <= 1e-12)
        {
            throw formatError(_source, entries[1].line,
                              "sample '" + entries[0].sampleLabel + "' has a singular calibration system");
        }

        const double gamma = 2.0 * (secondTotal - firstTotal) / denominator;
        const SampleClonalEntry& purityEntry = firstTotal != 2 ? entries[0] : entries[1];
        const double purityRdr = firstTotal != 2 ? firstRdr : secondRdr;
        double purity = (gamma * purityRdr - 2.0) / (purityEntry.total() - 2);
        if(1.0 <= purity && purity <= 1.0 + PURITY_TOLERANCE)
        {
            purity = 1.0;
        }
        else if(-PURITY_TOLERANCE <= purity && purity <= 0.0)
        {
            purity = 0.0;
        }
        if(!std::isfinite(gamma) || gamma <= 0.0)
        {
            throw formatError(_source, entries[1].line,
                              "sample '" + entries[0].sampleLabel + "' produced an invalid scaling factor");
        }
        if(!std::isfinite(purity) || purity < 0.0 || purity > 1.0)
        {
            throw formatError(_source, entries[1].line,
                              "sample '" + entries[0].sampleLabel + "' produced an invalid purity");
        }
        result.gamma[sample] = gamma;
        result.purity[sample] = purity;
    }
    return result;
}


int SampleClonal::maxTotalCopyNumber() const
{
    int result = 0;
    for(std::vector<SampleClonalEntry>::const_iterator entry = _entries.begin(); entry != _entries.end(); ++entry)
    {
        result = std::max(result, entry->total());
    }
    return result;
}


bool SampleClonal::samplesCompatible(const int left, const int right) const
{
    const std::vector<SampleClonalEntry> leftEntries = entriesForSample(left);
    const std::vector<SampleClonalEntry> rightEntries = entriesForSample(right);
    for(std::vector<SampleClonalEntry>::const_iterator leftEntry = leftEntries.begin();
        leftEntry != leftEntries.end(); ++leftEntry)
    {
        for(std::vector<SampleClonalEntry>::const_iterator rightEntry = rightEntries.begin();
            rightEntry != rightEntries.end(); ++rightEntry)
        {
            if(leftEntry->cluster == rightEntry->cluster &&
               (leftEntry->cnA != rightEntry->cnA || leftEntry->cnB != rightEntry->cnB))
            {
                return false;
            }
        }
    }
    return true;
}


void SampleClonal::validate(const int cmax,
                            const int base,
                            const bool ampdel,
                            const IntSet& ampdelExempt,
                            const IntArray& cnstates,
                            const DoubleArray& weights,
                            const DoubleArray& purities,
                            const double maxHomdelSize) const
{
    const double totalWeight = std::accumulate(weights.begin(), weights.end(), 0.0);
    std::map<int, std::set<CNState> > statesByCluster;
    std::map<int, std::vector<int> > alleleAByCluster;
    std::map<int, std::vector<int> > alleleBByCluster;
    for(std::vector<SampleClonalEntry>::const_iterator entry = _entries.begin(); entry != _entries.end(); ++entry)
    {
        if(entry->sample < 0 || entry->sample >= static_cast<int>(purities.size()) || purities[entry->sample] <= TOL)
        {
            continue;
        }
        if(entry->total() > cmax)
        {
            throw formatError(_source, entry->line,
                              "sample-clonal total copy number exceeds copy-number bound " +
                              std::to_string(cmax));
        }
        if(entry->total() == 0 && totalWeight > 0.0 && weights[entry->cluster] > 0.0 &&
           weights[entry->cluster] / totalWeight >= maxHomdelSize)
        {
            throw formatError(_source, entry->line,
                              "large cluster '" + entry->clusterLabel +
                              "' cannot have required tumor copy number (0, 0)");
        }
        statesByCluster[entry->cluster].insert(CNState(entry->cnA, entry->cnB));
        alleleAByCluster[entry->cluster].push_back(entry->cnA);
        alleleBByCluster[entry->cluster].push_back(entry->cnB);
    }

    for(std::map<int, std::set<CNState> >::const_iterator states = statesByCluster.begin();
        states != statesByCluster.end(); ++states)
    {
        const int limit = segmentCnstates(cnstates, states->first);
        if(limit > 0 && static_cast<int>(states->second.size()) > limit)
        {
            throw _source + ": sample-clonal states require more copy-number states than allowed";
        }
    }
    if(ampdel)
    {
        const std::map<int, std::vector<int> >* alleles[2] = {&alleleAByCluster, &alleleBByCluster};
        for(int allele = 0; allele < 2; ++allele)
        {
            for(std::map<int, std::vector<int> >::const_iterator values = alleles[allele]->begin();
                values != alleles[allele]->end(); ++values)
            {
                if(ampdelExempt.count(values->first) > 0)
                {
                    continue;
                }
                bool below = false;
                bool above = false;
                for(std::vector<int>::const_iterator value = values->second.begin(); value != values->second.end(); ++value)
                {
                    below = below || *value < base;
                    above = above || *value > base;
                }
                if(below && above)
                {
                    throw _source + ": sample-clonal states force both amplification and deletion for one allele";
                }
            }
        }
    }
}


void SampleClonal::verify(const IntMatrix& A,
                          const IntMatrix& B,
                          const DoubleMatrix& U,
                          const DoubleArray& purities) const
{
    for(int sample = 0; sample < static_cast<int>(purities.size()); ++sample)
    {
        const double expectedNormal = 1.0 - purities[sample];
        if(std::fabs(U[0][sample] - expectedNormal) > TOL)
        {
            throw _source + ": inferred normal proportion violates sample-clonal purity for sample index " +
                  std::to_string(sample);
        }
        const std::vector<SampleClonalEntry> sampleEntries = entriesForSample(sample);
        for(int clone = 1; clone < static_cast<int>(U.size()); ++clone)
        {
            if(U[clone][sample] <= TOL)
            {
                continue;
            }
            for(std::vector<SampleClonalEntry>::const_iterator entry = sampleEntries.begin();
                entry != sampleEntries.end(); ++entry)
            {
                if(A[entry->cluster][clone] != entry->cnA || B[entry->cluster][clone] != entry->cnB)
                {
                    throw _source + ": inferred present clone violates sample-clonal state for sample '" +
                          entry->sampleLabel + "', cluster '" + entry->clusterLabel + "'";
                }
            }
        }
    }
}
