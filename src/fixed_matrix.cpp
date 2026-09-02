#include "fixed_matrix.h"

#include <cerrno>

namespace
{
const std::string FORMAT_MARKER = "# HATCHET_FIXED_ENTRIES 1";
const std::string HEADER = "matrix\trow\tcolumn\tvalue";

std::string formatError(const std::string& source, const int line, const std::string& message)
{
    return source + ":" + std::to_string(line) + ": " + message;
}

std::string strip(const std::string& value)
{
    return rtrim(ltrim(value));
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

int parsePositiveInteger(const std::string& value,
                         const std::string& source,
                         const int line,
                         const std::string& description)
{
    if(value.empty())
    {
        throw formatError(source, line, description + " must be an integer");
    }

    errno = 0;
    char* end = NULL;
    const long parsed = std::strtol(value.c_str(), &end, 10);
    if(errno != 0 || end == value.c_str() || *end != '\0' || parsed < 0 || parsed > std::numeric_limits<int>::max())
    {
        throw formatError(source, line, description + " must be a nonnegative integer");
    }
    return static_cast<int>(parsed);
}

int parseClone(const std::string& label, const int clones, const std::string& source, const int line)
{
    if(label == "normal")
    {
        return 0;
    }
    if(label.size() <= 5 || label.substr(0, 5) != "clone" || label[5] < '1' || label[5] > '9')
    {
        throw formatError(source, line, "invalid clone label '" + label + "'; expected normal or cloneN");
    }
    for(std::string::size_type index = 6; index < label.size(); ++index)
    {
        if(label[index] < '0' || label[index] > '9')
        {
            throw formatError(source, line, "invalid clone label '" + label + "'; expected normal or cloneN");
        }
    }

    const int clone = parsePositiveInteger(label.substr(5), source, line, "clone index");
    if(clone < 1 || clone >= clones)
    {
        throw formatError(source, line, "clone label '" + label + "' is outside the declared clone count");
    }
    return clone;
}

double parseProportion(const std::string& value, const std::string& source, const int line)
{
    if(value.empty())
    {
        throw formatError(source, line, "proportion value must be numeric");
    }

    errno = 0;
    char* end = NULL;
    const double parsed = std::strtod(value.c_str(), &end);
    if(errno != 0 || end == value.c_str() || *end != '\0')
    {
        throw formatError(source, line, "proportion value '" + value + "' must be numeric");
    }
    if(!std::isfinite(parsed))
    {
        throw formatError(source, line, "proportion values must be finite");
    }
    if(parsed < 0.0 || parsed > 1.0)
    {
        throw formatError(source, line, "proportion values must be in [0, 1]");
    }
    return parsed;
}
}

FixedMatrix::FixedMatrix()
    : _clones(0)
    , _source()
    , _entries()
{
}

FixedMatrix FixedMatrix::read(const std::string& filename,
                              const int clones,
                              const std::map<std::string, int>& clusterToIdx,
                              const std::map<std::string, int>& sampleToIdx)
{
    std::ifstream input(filename.c_str());
    if(!input.good())
    {
        throw "The specified fixed-matrix file \"" + filename + "\" does not exist or cannot be opened!";
    }

    FixedMatrix result;
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
    bool clonesSeen = false;
    std::set<std::tuple<int, int, int> > coordinates;
    int lineNumber = 1;
    while(std::getline(input, line))
    {
        ++lineNumber;
        if(!line.empty() && line.back() == '\r')
        {
            line.pop_back();
        }
        if(strip(line).empty())
        {
            continue;
        }

        const std::string stripped = strip(line);
        if(stripped[0] == '#')
        {
            const std::string clonesPrefix = "# clones";
            if(stripped.substr(0, clonesPrefix.size()) == clonesPrefix)
            {
                if(headerSeen)
                {
                    throw formatError(filename, lineNumber, "clone-count metadata must appear before the header");
                }
                if(clonesSeen)
                {
                    throw formatError(filename, lineNumber, "duplicate clone-count metadata");
                }
                const std::string cloneCount = strip(stripped.substr(clonesPrefix.size()));
                if(cloneCount.empty() ||
                   std::find_if(cloneCount.begin(), cloneCount.end(),
                                [](char value) { return value < '0' || value > '9'; }) != cloneCount.end())
                {
                    throw formatError(filename, lineNumber, "clone count must be a nonnegative integer");
                }
                result._clones = parsePositiveInteger(cloneCount, filename, lineNumber, "clone count");
                if(result._clones < 2)
                {
                    throw formatError(filename, lineNumber, "the clone count must be at least 2, including normal");
                }
                if(result._clones != clones)
                {
                    throw formatError(filename, lineNumber,
                                      "declared clone count does not match the solver clone count " +
                                      std::to_string(clones));
                }
                clonesSeen = true;
            }
            continue;
        }

        if(!headerSeen)
        {
            if(!clonesSeen)
            {
                throw formatError(filename, lineNumber, "clone-count metadata must appear before the header");
            }
            if(line != HEADER)
            {
                throw formatError(filename, lineNumber, "expected tab-separated header '" + HEADER + "'");
            }
            headerSeen = true;
            continue;
        }

        const std::vector<std::string> fields = splitTabs(line);
        if(fields.size() != 4)
        {
            throw formatError(filename, lineNumber, "expected 4 tab-separated fields");
        }

        FixedMatrixEntry entry;
        entry.line = lineNumber;
        entry.rowLabel = fields[1];
        entry.columnLabel = fields[2];
        const std::string value = strip(fields[3]);
        if(fields[0] == "A" || fields[0] == "B")
        {
            entry.matrix = fields[0] == "A" ? MATRIX_A : MATRIX_B;
            std::map<std::string, int>::const_iterator cluster = clusterToIdx.find(entry.rowLabel);
            if(cluster == clusterToIdx.end())
            {
                throw formatError(filename, lineNumber, "unknown SEG cluster identifier '" + entry.rowLabel + "'");
            }
            entry.row = cluster->second;
            entry.column = parseClone(entry.columnLabel, clones, filename, lineNumber);
            entry.value = parsePositiveInteger(value, filename, lineNumber, "copy-number value");
            if(entry.column == 0 && entry.value != 1.0)
            {
                throw formatError(filename, lineNumber,
                                  fields[0] + " entries for the normal clone must equal 1");
            }
        }
        else if(fields[0] == "U")
        {
            entry.matrix = MATRIX_U;
            entry.row = parseClone(entry.rowLabel, clones, filename, lineNumber);
            std::map<std::string, int>::const_iterator sample = sampleToIdx.find(entry.columnLabel);
            if(sample == sampleToIdx.end())
            {
                throw formatError(filename, lineNumber, "unknown SEG sample identifier '" + entry.columnLabel + "'");
            }
            entry.column = sample->second;
            entry.value = parseProportion(value, filename, lineNumber);
        }
        else
        {
            throw formatError(filename, lineNumber, "unknown matrix '" + fields[0] + "'; expected A, B, or U");
        }

        const std::tuple<int, int, int> coordinate(entry.matrix, entry.row, entry.column);
        if(coordinates.count(coordinate) > 0)
        {
            throw formatError(filename, lineNumber, "duplicate fixed matrix entry");
        }
        coordinates.insert(coordinate);
        result._entries.push_back(entry);
    }

    if(!clonesSeen)
    {
        throw formatError(filename, lineNumber, "missing clone-count metadata");
    }
    if(!headerSeen)
    {
        throw formatError(filename, lineNumber, "missing matrix entry header");
    }
    if(result._entries.empty())
    {
        throw formatError(filename, lineNumber, "the file must contain at least one fixed entry");
    }
    return result;
}

void FixedMatrix::validate(const int cmax,
                           const double mu,
                           const int base,
                           const bool ampdel,
                           const int cnstates,
                           const CNMap& copyNumbers,
                           const DoubleArray& weights) const
{
    std::map<std::tuple<int, int, int>, const FixedMatrixEntry*> values;
    for(std::vector<FixedMatrixEntry>::const_iterator it = _entries.begin(); it != _entries.end(); ++it)
    {
        values[std::make_tuple(it->matrix, it->row, it->column)] = &(*it);
        if(it->matrix == MATRIX_U)
        {
            if(it->row > 0 && TOL < it->value && it->value < mu - TOL)
            {
                throw formatError(_source, it->line,
                                  "fixed tumor proportion is between 0 and the minimum clone proportion");
            }
            continue;
        }

        int upper = cmax;
        CNMap::const_iterator clonal = copyNumbers.find(it->row);
        if(clonal != copyNumbers.end())
        {
            upper = std::max(upper, clonal->second.first + clonal->second.second);
            if(it->column > 0)
            {
                const int expected = it->matrix == MATRIX_A ? clonal->second.first : clonal->second.second;
                if(static_cast<int>(it->value) != expected)
                {
                    throw formatError(_source, it->line, "fixed copy number conflicts with a clonal copy number");
                }
            }
        }
        if(it->value > upper)
        {
            throw formatError(_source, it->line, "fixed copy number exceeds the copy-number bound");
        }
    }

    for(int segment = 0; segment < static_cast<int>(weights.size()); ++segment)
    {
        int upper = cmax;
        CNMap::const_iterator clonal = copyNumbers.find(segment);
        if(clonal != copyNumbers.end())
        {
            upper = std::max(upper, clonal->second.first + clonal->second.second);
        }
        std::set<CNState> states;
        for(int clone = 0; clone < _clones; ++clone)
        {
            std::map<std::tuple<int, int, int>, const FixedMatrixEntry*>::const_iterator a =
                values.find(std::make_tuple(MATRIX_A, segment, clone));
            std::map<std::tuple<int, int, int>, const FixedMatrixEntry*>::const_iterator b =
                values.find(std::make_tuple(MATRIX_B, segment, clone));
            if(a == values.end() || b == values.end())
            {
                continue;
            }

            const int total = static_cast<int>(a->second->value + b->second->value);
            if(total > upper)
            {
                throw formatError(_source, b->second->line, "fixed A+B total exceeds the copy-number bound");
            }
            if(clone > 0)
            {
                states.insert(CNState(static_cast<int>(a->second->value), static_cast<int>(b->second->value)));
                const double totalWeight = std::accumulate(weights.begin(), weights.end(), 0.0);
                if(totalWeight > 0.0 && weights[segment] / totalWeight >= 0.005 && total == 0)
                {
                    throw formatError(_source, b->second->line,
                                      "a large cluster cannot have fixed tumor copy number (0, 0)");
                }
            }
        }
        if(cnstates > 0 && static_cast<int>(states.size()) > cnstates)
        {
            throw formatError(_source, 1, "fixed entries exceed the maximum number of copy-number states");
        }

        if(ampdel && clonal == copyNumbers.end())
        {
            for(int matrix = MATRIX_A; matrix <= MATRIX_B; ++matrix)
            {
                bool below = false;
                bool above = false;
                int errorLine = 1;
                for(int clone = 1; clone < _clones; ++clone)
                {
                    std::map<std::tuple<int, int, int>, const FixedMatrixEntry*>::const_iterator entry =
                        values.find(std::make_tuple(matrix, segment, clone));
                    if(entry != values.end())
                    {
                        below = below || entry->second->value < base;
                        above = above || entry->second->value > base;
                        errorLine = std::max(errorLine, entry->second->line);
                    }
                }
                if(below && above)
                {
                    throw formatError(_source, errorLine,
                                      "fixed entries force both amplification and deletion for one allele");
                }
            }
        }
    }

    int sampleCount = 0;
    for(std::vector<FixedMatrixEntry>::const_iterator it = _entries.begin(); it != _entries.end(); ++it)
    {
        if(it->matrix == MATRIX_U)
        {
            sampleCount = std::max(sampleCount, it->column + 1);
        }
    }
    for(int sample = 0; sample < sampleCount; ++sample)
    {
        double fixedSum = 0.0;
        int fixedCount = 0;
        int errorLine = 1;
        bool normalFree = true;
        for(int clone = 0; clone < _clones; ++clone)
        {
            std::map<std::tuple<int, int, int>, const FixedMatrixEntry*>::const_iterator entry =
                values.find(std::make_tuple(MATRIX_U, clone, sample));
            if(entry != values.end())
            {
                fixedSum += entry->second->value;
                ++fixedCount;
                normalFree = normalFree && clone != 0;
                errorLine = std::max(errorLine, entry->second->line);
            }
        }
        if(fixedSum > 1.0 + TOL)
        {
            throw formatError(_source, errorLine, "fixed U entries for one sample sum to more than 1");
        }
        const double residual = 1.0 - fixedSum;
        if(fixedCount == _clones && std::fabs(residual) > TOL)
        {
            throw formatError(_source, errorLine, "fully fixed U entries for one sample must sum to 1");
        }
        if(fixedCount < _clones && !normalFree && TOL < residual && residual < mu - TOL)
        {
            throw formatError(_source, errorLine,
                              "remaining U proportion is below the minimum clone proportion");
        }
    }

    for(int clone = 1; clone < _clones - 1; ++clone)
    {
        bool complete = true;
        double left = 0.0;
        double right = 0.0;
        int errorLine = 1;
        for(int segment = 0; segment < static_cast<int>(weights.size()); ++segment)
        {
            for(int matrix = MATRIX_A; matrix <= MATRIX_B; ++matrix)
            {
                std::map<std::tuple<int, int, int>, const FixedMatrixEntry*>::const_iterator leftEntry =
                    values.find(std::make_tuple(matrix, segment, clone));
                std::map<std::tuple<int, int, int>, const FixedMatrixEntry*>::const_iterator rightEntry =
                    values.find(std::make_tuple(matrix, segment, clone + 1));
                if(leftEntry == values.end() || rightEntry == values.end())
                {
                    complete = false;
                    break;
                }
                left += leftEntry->second->value * (segment + 1);
                right += rightEntry->second->value * (segment + 1);
                errorLine = std::max(errorLine, rightEntry->second->line);
            }
            if(!complete)
            {
                break;
            }
        }
        if(complete && left > right)
        {
            throw formatError(_source, errorLine, "fully fixed clone burdens violate canonical clone order");
        }
    }
}

void FixedMatrix::applyCopyNumbers(IntMatrix& A, IntMatrix& B) const
{
    for(std::vector<FixedMatrixEntry>::const_iterator it = _entries.begin(); it != _entries.end(); ++it)
    {
        if(it->matrix == MATRIX_A)
        {
            A[it->row][it->column] = static_cast<int>(it->value);
        }
        else if(it->matrix == MATRIX_B)
        {
            B[it->row][it->column] = static_cast<int>(it->value);
        }
    }
}

void FixedMatrix::applyProportions(DoubleMatrix& U) const
{
    for(std::vector<FixedMatrixEntry>::const_iterator it = _entries.begin(); it != _entries.end(); ++it)
    {
        if(it->matrix == MATRIX_U)
        {
            U[it->row][it->column] = it->value;
        }
    }
}

void FixedMatrix::conditionProportions(DoubleMatrix& U, const double mu) const
{
    if(_entries.empty())
    {
        return;
    }

    const int samples = U.empty() ? 0 : static_cast<int>(U[0].size());
    for(int sample = 0; sample < samples; ++sample)
    {
        std::set<int> fixedClones;
        double fixedSum = 0.0;
        for(std::vector<FixedMatrixEntry>::const_iterator it = _entries.begin(); it != _entries.end(); ++it)
        {
            if(it->matrix == MATRIX_U && it->column == sample)
            {
                fixedClones.insert(it->row);
                fixedSum += it->value;
            }
        }
        if(fixedClones.empty())
        {
            continue;
        }

        std::vector<int> eligible;
        const double residual = 1.0 - fixedSum;
        for(int clone = 0; clone < _clones; ++clone)
        {
            if(fixedClones.count(clone) == 0)
            {
                U[clone][sample] = 0.0;
                if(residual >= mu - TOL || clone == 0)
                {
                    eligible.push_back(clone);
                }
            }
        }
        applyProportions(U);
        if(residual > TOL)
        {
            if(eligible.empty())
            {
                throw formatError(_source, 1, "cannot allocate the remaining U proportion");
            }
            std::uniform_int_distribution<int> select(0, static_cast<int>(eligible.size()) - 1);
            U[eligible[select(g_rng)]][sample] = residual;
        }
    }
}

void FixedMatrix::verify(const IntMatrix& A, const IntMatrix& B, const DoubleMatrix& U) const
{
    for(std::vector<FixedMatrixEntry>::const_iterator it = _entries.begin(); it != _entries.end(); ++it)
    {
        double actual = 0.0;
        if(it->matrix == MATRIX_A)
        {
            actual = A[it->row][it->column];
        }
        else if(it->matrix == MATRIX_B)
        {
            actual = B[it->row][it->column];
        }
        else
        {
            actual = U[it->row][it->column];
        }
        if(std::fabs(actual - it->value) > TOL)
        {
            throw formatError(_source, it->line,
                              "solver returned " + std::to_string(actual) + " for a fixed entry; expected " +
                              std::to_string(it->value));
        }
    }
}
