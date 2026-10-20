#include "coordinate_descent.h"

CoordinateDescent::CoordinateDescent(const DoubleMatrix& FA,
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
                                     const VERBOSITY_t v)
    : _FA(FA)
    , _FB(FB)
    , _bins(bins)
    , _m(m)
    , _k(k)
    , _n(n)
    , _cmax(cmax)
    , _d(d)
    , _mu(mu)
    , _maxHomdelSize(maxHomdelSize)
    , _base(base)
    , _ampdel(ampdel)
    , _cn(cn)
    , _fixedMatrix(fixedMatrix)
    , _sampleClonal(sampleClonal)
    , _purities(purities)
    , _iterConvergence(iterConvergence)
    , _maxIter(maxIter)
    , _nrSeeds(nrSeeds)
    , _nrWorkers(nrWorkers)
    , _nrILPthreads(nrILPthreads)
    , _timeLimit(timeLimit)
    , _memoryLimit(memoryLimit)
    , _sizeBubbles(-1)
    , _seeds()
    , _seedidx(0)
    , _threads()
    , _mutexTasks()
    , _mutexResults()
    , _objs()
    , _CAs()
    , _CBs()
    , _Us()
    , _best(-1)
    , _v(v)
{
    if(mu <= 0.1) {
        _sizeBubbles = 10;
    } else if(mu <= 0.15) {
        _sizeBubbles = 6;
    } else if(mu <= 0.2) {
        _sizeBubbles = 5;
    } else {
        _sizeBubbles = 3;
    }
}


void CoordinateDescent::run()
{
    for (int i = 0; i < _nrSeeds; ++i)
    {
        _seeds.push_back(buildRandomU());
    }

    log("Detailed data about coordinate descence is going to be written in: coordinate_descence.tsv", VERBOSITY_t::DEBUGGING, _v);
    if(_v >= VERBOSITY_t::VERBOSE)
    {
        std::ofstream coordinate_out("coordinate_descence.tsv", std::ofstream::out);
        coordinate_out << "";
    }
    log("Coordinate Descence {", VERBOSITY_t::VERBOSE, _v);
    for(int i = 0; i < _nrWorkers; ++i)
    {
        _threads.emplace_back(&CoordinateDescent::runWorker, this);
    }

    for(auto& thread : _threads)
    {
        thread.join();
    }

    log("}", VERBOSITY_t::VERBOSE, _v);
    if(_objs.empty())
    {
        throw "Not a single feasible sample-clonal coordinate-descent seed was found!";
    }
    _best = std::distance(_objs.begin(), std::min_element(_objs.begin(), _objs.end()));
}


DoubleMatrix CoordinateDescent::buildRandomU()
{
    if(!_sampleClonal.empty())
    {
        return buildSampleClonalU();
    }

    DoubleMatrix U;
    for(unsigned int p = 0; p < _k; ++p)
    {
        std::uniform_int_distribution<> uni(1, _n);
        int non_zeros = std::max(uni(g_rng), uni(g_rng));
        non_zeros = std::min(non_zeros, _sizeBubbles);
        U.push_back(buildPartitionVector(non_zeros));
    }
    DoubleMatrix T(transpose(U));
    _fixedMatrix.conditionProportions(T, _mu);
    log("SEED " + std::to_string(_seeds.size()) + ":\n" + toString(T), VERBOSITY_t::DEBUGGING, _v);
    return T;
}


DoubleMatrix CoordinateDescent::buildSampleClonalU()
{
    DoubleMatrix U(_n, DoubleArray(_k, 0.0));
    std::vector<std::set<int> > fixedClones(_k);
    std::vector<std::set<int> > assignedSamples(_n);
    for(int sample = 0; sample < _k; ++sample)
    {
        U[0][sample] = 1.0 - _purities[sample];
    }

    const std::vector<FixedMatrixEntry>& fixedEntries = _fixedMatrix.entries();
    for(std::vector<FixedMatrixEntry>::const_iterator entry = fixedEntries.begin(); entry != fixedEntries.end(); ++entry)
    {
        if(entry->matrix != MATRIX_U)
        {
            continue;
        }
        if(entry->row == 0)
        {
            if(std::fabs(entry->value - U[0][entry->column]) > TOL)
            {
                throw "Fixed normal proportion conflicts with sample-clonal purity";
            }
            continue;
        }
        fixedClones[entry->column].insert(entry->row);
        U[entry->row][entry->column] = entry->value;
        if(entry->value > TOL)
        {
            assignedSamples[entry->row].insert(entry->column);
        }
    }

    const std::function<bool(int, int)> cloneCompatible =
        [this, &fixedEntries](const int clone, const int sample)
        {
            const std::vector<SampleClonalEntry> entries = _sampleClonal.entriesForSample(sample);
            for(std::vector<SampleClonalEntry>::const_iterator required = entries.begin();
                required != entries.end(); ++required)
            {
                for(std::vector<FixedMatrixEntry>::const_iterator fixed = fixedEntries.begin();
                    fixed != fixedEntries.end(); ++fixed)
                {
                    if(fixed->row != required->cluster || fixed->column != clone)
                    {
                        continue;
                    }
                    if((fixed->matrix == MATRIX_A && static_cast<int>(fixed->value) != required->cnA) ||
                       (fixed->matrix == MATRIX_B && static_cast<int>(fixed->value) != required->cnB))
                    {
                        return false;
                    }
                }
            }
            return true;
        };

    DoubleArray residual(_k, 0.0);
    std::vector<int> sampleOrder;
    for(int sample = 0; sample < _k; ++sample)
    {
        double fixedTumor = 0.0;
        for(int clone = 1; clone < _n; ++clone)
        {
            fixedTumor += U[clone][sample];
        }
        residual[sample] = _purities[sample] - fixedTumor;
        if(residual[sample] < -TOL)
        {
            throw "Fixed tumor proportions exceed sample-clonal purity";
        }
        if(residual[sample] > TOL)
        {
            if(_mu > 0.0 && residual[sample] < _mu - TOL)
            {
                throw "Remaining sample-clonal purity is below the minimum clone proportion";
            }
            sampleOrder.push_back(sample);
        }
    }
    std::shuffle(sampleOrder.begin(), sampleOrder.end(), g_rng);
    std::sort(sampleOrder.begin(), sampleOrder.end(),
              [this](const int left, const int right)
              {
                  int leftConflicts = 0;
                  int rightConflicts = 0;
                  for(int sample = 0; sample < _k; ++sample)
                  {
                      leftConflicts += _sampleClonal.samplesCompatible(left, sample) ? 0 : 1;
                      rightConflicts += _sampleClonal.samplesCompatible(right, sample) ? 0 : 1;
                  }
                  return leftConflicts > rightConflicts;
              });

    const std::function<bool(int, int)> cloneAvailable =
        [this, &fixedClones, &assignedSamples, &cloneCompatible](const int clone, const int sample)
        {
            bool compatible = fixedClones[sample].count(clone) == 0 && cloneCompatible(clone, sample);
            for(std::set<int>::const_iterator assigned = assignedSamples[clone].begin();
                compatible && assigned != assignedSamples[clone].end(); ++assigned)
            {
                compatible = _sampleClonal.samplesCompatible(sample, *assigned);
            }
            return compatible;
        };

    std::map<int, int> assignment;
    std::function<bool(std::size_t)> assignSample;
    assignSample =
        [this, &sampleOrder, &assignedSamples, &assignment, &cloneAvailable, &assignSample]
        (const std::size_t position) -> bool
        {
            if(position == sampleOrder.size())
            {
                return true;
            }
            const int sample = sampleOrder[position];
            IntArray candidates;
            for(int clone = 1; clone < _n; ++clone)
            {
                if(cloneAvailable(clone, sample))
                {
                    candidates.push_back(clone);
                }
            }
            std::shuffle(candidates.begin(), candidates.end(), g_rng);
            for(IntArray::const_iterator clone = candidates.begin(); clone != candidates.end(); ++clone)
            {
                assignedSamples[*clone].insert(sample);
                assignment[sample] = *clone;
                if(assignSample(position + 1))
                {
                    return true;
                }
                assignedSamples[*clone].erase(sample);
                assignment.erase(sample);
            }
            return false;
        };

    if(!assignSample(0))
    {
        throw "No feasible tumor-clone assignment exists for the sample-clonal states";
    }

    // Spread each sample over a random number of compatible clones, as buildRandomU does: with one clone per
    // sample, clones beyond the number of samples would be absent from every seed and never optimized.
    std::map<int, IntArray> clonesBySample;
    IntArray extensionOrder;
    for(std::map<int, int>::const_iterator assigned = assignment.begin(); assigned != assignment.end(); ++assigned)
    {
        clonesBySample[assigned->first].push_back(assigned->second);
        extensionOrder.push_back(assigned->first);
    }
    std::shuffle(extensionOrder.begin(), extensionOrder.end(), g_rng);
    for(IntArray::const_iterator sample = extensionOrder.begin(); sample != extensionOrder.end(); ++sample)
    {
        std::uniform_int_distribution<> uni(1, _n - 1);
        int numParts = std::max(uni(g_rng), uni(g_rng));
        if(_mu > 0.0)
        {
            numParts = std::min(numParts, std::max(1, static_cast<int>(std::floor((residual[*sample] + TOL) / _mu))));
        }
        IntArray& clones = clonesBySample[*sample];
        IntArray candidates;
        for(int clone = 1; clone < _n; ++clone)
        {
            if(std::find(clones.begin(), clones.end(), clone) == clones.end() && cloneAvailable(clone, *sample))
            {
                candidates.push_back(clone);
            }
        }
        std::shuffle(candidates.begin(), candidates.end(), g_rng);
        for(IntArray::const_iterator clone = candidates.begin();
            clone != candidates.end() && static_cast<int>(clones.size()) < numParts; ++clone)
        {
            assignedSamples[*clone].insert(*sample);
            clones.push_back(*clone);
        }
    }

    for(std::map<int, IntArray>::const_iterator sample = clonesBySample.begin();
        sample != clonesBySample.end(); ++sample)
    {
        const DoubleArray shares = splitSamplePurity(residual[sample->first], static_cast<int>(sample->second.size()));
        for(std::size_t part = 0; part < shares.size(); ++part)
        {
            U[sample->second[part]][sample->first] = shares[part];
        }
    }
    return U;
}


DoubleArray CoordinateDescent::splitSamplePurity(const double residual, const int numParts)
{
    const double minShare = std::max(_mu, 0.0);
    const double free = std::max(residual - numParts * minShare, 0.0);
    std::uniform_real_distribution<double> uniform(0.0, free);
    DoubleArray cuts(numParts - 1);
    for(double& cut : cuts)
    {
        cut = uniform(g_rng);
    }
    std::sort(cuts.begin(), cuts.end());

    DoubleArray shares(numParts, minShare);
    double previous = 0.0;
    double assigned = 0.0;
    for(int part = 0; part < numParts - 1; ++part)
    {
        shares[part] += cuts[part] - previous;
        previous = cuts[part];
        assigned += shares[part];
    }
    shares[numParts - 1] = residual - assigned;
    return shares;
}


DoubleArray CoordinateDescent::buildPartitionVector(const int num_parts)
{
    assert(num_parts <= _sizeBubbles);
    IntArray positions(_n);
    std::iota(positions.begin(), positions.end(), 0);
    std::shuffle(positions.begin(), positions.end(), g_rng);

    DoubleArray result(_n, 0.0);
    if(num_parts > 1)
    {
        IntArray samples(_sizeBubbles - 1);
        std::iota(samples.begin(), samples.end(), 1);
        assert(samples.front() == 1 && samples.back() == _sizeBubbles - 1);
        std::shuffle(samples.begin(), samples.end(), g_rng);

        DoubleArray bubbles(num_parts - 1, 0.0);
        for(unsigned int i = 0; i < num_parts - 1; ++i)
        {
            bubbles[i] = (double)samples.back() / (double)_sizeBubbles;
            samples.pop_back();
        }
        std::sort(bubbles.begin(), bubbles.end());

        for(unsigned int i = 0; i < num_parts; ++i)
        {
            if(i == 0)
            {
                result[positions[i]] = bubbles[i];
            }
            else if(i == (num_parts-1))
            {
                result[positions[i]] = 1.0 - bubbles.back();
            }
            else
            {
                result[positions[i]] = bubbles[i] - bubbles[i-1];
            }
            if(_mu - TOL <= result[positions[i]] && result[positions[i]] < _mu)
            {
                result[positions[i]] = _mu;
            }
            assert(result[positions[i]] >= _mu);
        }
    }
    else {
        result[positions[0]] = 1.0;
    }
    return result;
}


void CoordinateDescent::runWorker()
{
    while(true)
    {
        int seed = -1;
        {
            std::lock_guard<std::mutex> lock(_mutexTasks);
            seed = _seedidx < _seeds.size() ? _seedidx++ : -1;
        }
        std::this_thread::sleep_for(std::chrono::seconds(1));
        if(seed > -1) {
            Worker worker(_FA, _FB, _bins, _m, _k, _n, _cmax, _d, _mu, _maxHomdelSize, _base, _ampdel, _cn, _fixedMatrix,
                          _sampleClonal, _purities,
                          _iterConvergence, _maxIter, _timeLimit, _memoryLimit, _nrILPthreads, _seeds[seed], seed, _v);
            double obj = worker.solve();

            if(!std::isfinite(obj))
            {
                continue;
            }

            {
                std::lock_guard<std::mutex> lock(_mutexResults);
                _objs.push_back(obj);
                _CAs.push_back(worker.getCA());
                _CBs.push_back(worker.getCB());
                _Us.push_back(worker.getU());
            }

            {
                std::lock_guard<std::mutex> lock(g_output_mutex);
                std::locale::global(std::locale("C"));
                log(std::to_string(obj) + "; ", VERBOSITY_t::VERBOSE, _v);
            }
        } else {
            break;
        }
    }
}
