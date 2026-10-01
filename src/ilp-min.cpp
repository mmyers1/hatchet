#include "ilp-min.h"

namespace
{
    /// Narrow [lower, upper] to a single value; returns false when the value lies outside the current bounds
    template<typename T>
    bool fixBounds(T& lower, T& upper, const T value)
    {
        if(value < lower || value > upper)
        {
            return false;
        }
        lower = value;
        upper = value;
        return true;
    }
}

IlpSubset::IlpSubset(const int n,
                     const int m,
                     const int k,
                     const int cmax,
                     const int d,
                     const double mu,
                     const double maxHomdelSize,
                     const int base,
                     const bool ampdel,
                     const CNMap& cn,
                     const DoubleMatrix& FA,
                     const DoubleMatrix& FB,
                     const DoubleArray& w,
                     const FixedMatrix& fixedMatrix,
                     const SampleClonal& sampleClonal,
                     const DoubleArray& purities,
                     const VERBOSITY_t v)
    : _n(n)
    , _m(m)
    , _k(k)
    , _cmax(cmax)
    , _d(d)
    , _mu(mu)
    , _maxHomdelSize(maxHomdelSize)
    , _base(base)
    , _ampdel(ampdel)
    , _cn(cn)
    , _FA(FA)
    , _FB(FB)
    , _w(w)
    , _fixedMatrix(fixedMatrix)
    , _sampleClonal(sampleClonal)
    , _purities(purities)
    , _v(v)
    , _barCA()
    , _barCB()
    , _barU()
    , _mode(MODE_t::FULL)
    , _M(floor(log2(cmax)) + 1)
    , _env()
    , _model(_env)
    , _yA()
    , _yB()
    , _fA()
    , _fB()
    , _cA()
    , _cB()
    , _bitcA()
    , _bitcB()
    , _u()
    , _x()
    , _adA()
    , _adB()
    , _vA()
    , _vB()
    , _barcIJA()
    , _barcIJB()
    , _objs()
    , _CAs()
    , _CBs()
    , _Us()
{
}


void IlpSubset::fixC(const IntMatrix& barCA, const IntMatrix& barCB)
{
    _barCA = barCA;
    _barCB = barCB;
    _fixedMatrix.applyCopyNumbers(_barCA, _barCB);
    _mode = MODE_t::UARCH;
    for(int s = 0; s < _m; ++s)
    {
        int upper = 0;
        if(_cn.count(s) > 0)
        {
            upper = _cn.at(s).first + _cn.at(s).second;
        } else {
            upper = _cmax;
        }

        int adA = 0;
        int adB = 0;
        for(int i = 0; i < _n; ++i)
        {
            if(i > 0)
            {
                assert(barCA[s][i] + barCB[s][i] <= upper);
                assert(barCA[s][i] >= 0 && barCB[s][i] >= 0);
                if(_ampdel)
                {
                    adA = adA == 0 && barCA[s][i] > _base ? 1 : adA;
                    adA = adA == 0 && barCA[s][i] < _base ? -1 : adA;
                    adB = adB == 0 && barCB[s][i] > _base ? 1 : adB;
                    adB = adB == 0 && barCB[s][i] < _base ? -1 : adB;

                    assert(adA <= 0 || barCA[s][i] >= _base);
                    assert(adA >= 0 || barCA[s][i] <= _base);
                    assert(adB <= 0 || barCB[s][i] >= _base);
                    assert(adB >= 0 || barCB[s][i] <= _base);
                }
            } else {
                assert(barCA[s][i] == 1);
                assert(barCB[s][i] == 1);
            }
        }
    }
}


void IlpSubset::init()
{
    computeFixedBounds();
    buildVariables();
    if(_d > 0 && (_mode == MODE_t::FULL || _mode == MODE_t::CARCH))
    {
        buildOptionalVariables();
    }

    buildConstraints();
    if(_d > 0 && (_mode == MODE_t::FULL || _mode == MODE_t::CARCH))
    {
        buildOptionalConstraints();
    }
    fixSampleClonal();
    addFixedValueConflicts();

    buildObjective();
    //test();
    _model.update();
}


void IlpSubset::hotStart(const IntMatrix& hCA, const IntMatrix& hCB, const DoubleMatrix& hU)
{
    if(_mode == MODE_t::UARCH)
    {
        throw "A Hot Start cannot be provided in UARCH mode when mixing proportions U are fixed!";
    }
    if(hCA.size() != _m || hCB.size() != _m)
    {
        throw "The rows of both hot starts must be equal to the current number m of segments!";
    }
    for(int s = 0; s < _m; ++s)
    {
        if(hCA[s].size() != _n || hCB[s].size() != _n)
        {
            throw "The columns of both hot starts must be equal to the current number n of clones!";
        }
        if(hCA[s][0] != 1 || hCB[s][0] != 1)
        {
            std::cout << toString(hCB) << std::endl;
            throw "The first clone must have only allele-specific copy numbers equal to 1";
        }
    }

    // Clone columns are used as given: the model has no clone-order constraint, and in CARCH mode
    // (and in FULL mode when hU is supplied) the column identities are tied to the proportions.
    for(int s = 0; s < _m; ++s)
    {
        for(int i = 0; i < _n; ++i)
        {
            _cA[s][i].set(GRB_DoubleAttr_Start, hCA[s][i]);
            _cB[s][i].set(GRB_DoubleAttr_Start, hCB[s][i]);
            if(!_bitcA.empty())
            {
                for(int b = 0; b < _M; ++b)
                {
                    _bitcA[b][s][i].set(GRB_DoubleAttr_Start, (hCA[s][i] >> b) & 1);
                    _bitcB[b][s][i].set(GRB_DoubleAttr_Start, (hCB[s][i] >> b) & 1);
                }
            }
        }

        if(_ampdel && _cn.count(s) == 0)
        {
            bool ampA = false;
            bool ampB = false;
            for(int i = 1; i < _n; ++i)
            {
                ampA = ampA || hCA[s][i] > _base;
                ampB = ampB || hCB[s][i] > _base;
            }
            _adA[s].set(GRB_DoubleAttr_Start, ampA ? 1.0 : 0.0);
            _adB[s].set(GRB_DoubleAttr_Start, ampB ? 1.0 : 0.0);
        }
    }

    if(_mode == MODE_t::FULL && !hU.empty())
    {
        if(hU.size() != _n)
        {
            throw "The rows of the proportion hot start must be equal to the current number n of clones!";
        }
        for(int i = 0; i < _n; ++i)
        {
            if(hU[i].size() != _k)
            {
                throw "The columns of the proportion hot start must be equal to the current number k of samples!";
            }
            for(int p = 0; p < _k; ++p)
            {
                _u[i][p].set(GRB_DoubleAttr_Start, hU[i][p]);
                if(i > 0 && !_x.empty())
                {
                    _x[i][p].set(GRB_DoubleAttr_Start, hU[i][p] > TOL ? 1.0 : 0.0);
                }
            }
        }
    }
}


int IlpSubset::solve(const int timeLimit, const int memoryLimit, const int nrThreads)
{
    return solve(timeLimit, memoryLimit, nrThreads, 1);
}


int IlpSubset::solve(const int timeLimit, const int memoryLimit, const int nrThreads, const int maxIter)
{
    int count = 0;
    try
    {
        if (timeLimit > 0)
        {
            _model.getEnv().set(GRB_DoubleParam_TimeLimit, timeLimit);
        }
        if (memoryLimit > 0)
        {
            _model.getEnv().set(GRB_DoubleParam_NodeLimit, memoryLimit);
        }
        if (nrThreads > 0)
        {
            _model.getEnv().set(GRB_IntParam_Threads, nrThreads);
        }
        if(_mode == MODE_t::CARCH || _mode == MODE_t::UARCH || _v <= VERBOSITY_t::ESSENTIAL)
        {
            _model.getEnv().set(GRB_IntParam_LogToConsole, 0);
        }

        bool proceed = true;
        do {
            _model.optimize();
            //printValues();
            int status = _model.get(GRB_IntAttr_Status);

            if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT)
            {
                _objs.push_back(_model.get(GRB_DoubleAttr_ObjVal));
                _CAs.push_back(getCA());
                _CBs.push_back(getCB());
                _Us.push_back(getU());
                if(_mode == MODE_t::FULL)
                {
                    buildNext();
                }
                proceed = true;
                ++count;
            }
            else if (status == GRB_INF_OR_UNBD)
            {
                std::cerr << "Model is infeasible or unbounded" << std::endl;
                _model.computeIIS();
                _model.write("iis.ilp");
                proceed = false;
            }
            else if (status == GRB_INFEASIBLE)
            {
                std::cerr << "Model is infeasible" << std::endl;
                _model.computeIIS();
                _model.write("iis.ilp");
                proceed = false;
            }
            else if (status == GRB_UNBOUNDED)
            {
                std::cerr << "Model is unbounded" << std::endl;
                _model.computeIIS();
                _model.write("iis.ilp");
                proceed = false;
            }
        } while(proceed && (maxIter == 0 || count < maxIter));

        return count;
    }
    catch (const GRBException& e)
    {
        std::cerr << "Error code = " << e.getErrorCode() << std::endl;
        std::cerr << e.getMessage() << std::endl;

        return -1;
    }
    catch (...)
    {
        return -1;
    }

    return 0;
}


void IlpSubset::exportModel(const std::string& filename)
{
    _model.write(filename);
}


void IlpSubset::printValues()
{
    for(int s = 0; s < _m; ++s)
    {
        for(int p = 0; p < _k; ++p)
        {
            std::cerr << _yA[s][p].get(GRB_StringAttr_VarName) << " = " << _yA[s][p].get(GRB_DoubleAttr_X) << std::endl;
            std::cerr << _yB[s][p].get(GRB_StringAttr_VarName) << " = " << _yB[s][p].get(GRB_DoubleAttr_X) << std::endl;
        }
    }

    for(int s = 0; s < _m; ++s)
    {
        for(int p = 0; p < _k; ++p)
        {
            std::cerr << _fA[s][p].get(GRB_StringAttr_VarName) << " = " << _fA[s][p].get(GRB_DoubleAttr_X) << std::endl;
            std::cerr << _fB[s][p].get(GRB_StringAttr_VarName) << " = " << _fB[s][p].get(GRB_DoubleAttr_X) << std::endl;
        }
    }

    if(_mode == MODE_t::FULL || _mode == MODE_t::CARCH)
    {
        for(int s = 0; s < _m; ++s)
        {
            for(int i = 0; i < _n; ++i)
            {
                std::cerr << _cA[s][i].get(GRB_StringAttr_VarName) << " = " << _cA[s][i].get(GRB_DoubleAttr_X) << std::endl;
                std::cerr << _cB[s][i].get(GRB_StringAttr_VarName) << " = " << _cB[s][i].get(GRB_DoubleAttr_X) << std::endl;
            }
        }

        for(int b = 0; b < _M; ++b)
        {
            for(int s = 0; s < _m; ++s)
            {
                for(int i = 0; i < _n; ++i)
                {
                    std::cerr << _bitcA[b][s][i].get(GRB_StringAttr_VarName) << " = " << _bitcA[b][s][i].get(GRB_DoubleAttr_X) << std::endl;
                    std::cerr << _bitcB[b][s][i].get(GRB_StringAttr_VarName) << " = " << _bitcB[b][s][i].get(GRB_DoubleAttr_X) << std::endl;
                }
            }
        }
    }

    if(_mode == MODE_t::FULL || _mode == MODE_t::UARCH)
    {
        for(int i =0; i < _n; ++i)
        {
            for(int p = 0; p < _k; ++p)
            {
                std::cerr << _u[i][p].get(GRB_StringAttr_VarName) << " = " << _u[i][p].get(GRB_DoubleAttr_X) << std::endl;
            }
        }
    }

    if(_mode == MODE_t::FULL)
    {
        for(int b = 0; b < _M; ++b)
        {
            for(int s = 0; s < _m; ++s)
            {
                for(int i = 0; i < _n; ++i)
                {
                    for(int p = 0; p < _k; ++p)
                    {
                        std::cerr << _vA[b][s][i][p].get(GRB_StringAttr_VarName) << " = " << _vA[b][s][i][p].get(GRB_DoubleAttr_X) << std::endl;
                        std::cerr << _vB[b][s][i][p].get(GRB_StringAttr_VarName) << " = " << _vB[b][s][i][p].get(GRB_DoubleAttr_X) << std::endl;
                    }
                }
            }
        }
    }
}


std::pair<IntMatrix, IntMatrix> IlpSubset::firstHotStart(const DoubleMatrix &FA, const DoubleMatrix &FB, const int m, const int k, const int n,
                                                         const int cmax, const int d, const int base, const bool ampdel, const CNMap& cn)
{
    IntMatrix targetA;
    IntMatrix targetB;
    if(d > 0)
    {
        IntMatrix optionsA(m);
        IntMatrix optionsB(m);
        for(int s = 0; s < m; ++s)
        {
            optionsA[s] = IntArray();
            optionsB[s] = IntArray();
            for(int p = 0; p < k; ++p)
            {
                optionsA[s].push_back((int)std::round(FA[s][p]));
                optionsB[s].push_back((int)std::round(FB[s][p]));
            }
            std::shuffle(optionsA[s].begin(), optionsA[s].end(), g_rng);
            while(optionsA[s].size() > d)
            {
                optionsA[s].pop_back();
                std::shuffle(optionsA[s].begin(), optionsA[s].end(), g_rng);
            }
            std::shuffle(optionsB[s].begin(), optionsB[s].end(), g_rng);
            while(optionsB[s].size() > d)
            {
                optionsB[s].pop_back();
                std::shuffle(optionsB[s].begin(), optionsB[s].end(), g_rng);
            }
        }
        targetA = IntMatrix(m);
        targetB = IntMatrix(m);
        for(int s = 0; s < m; ++s)
        {
            targetA[s] = IntArray(d);
            targetB[s] = IntArray(d);
            for(int l = 0; l < d; ++l)
            {
                targetA[s][l] = optionsA[s][l % optionsA[s].size()];
                targetB[s][l] = optionsB[s][l % optionsB[s].size()];
            }
        }
    } else {
        targetA = IntMatrix(m);
        targetB = IntMatrix(m);
        for(int s = 0; s < m; ++s)
        {
            targetA[s] = IntArray(k);
            targetB[s] = IntArray(k);
            for(int p = 0; p < k; ++p)
            {
                targetA[s][p] = std::round(FA[s][p]);
                targetB[s][p] = std::round(FB[s][p]);
            }
        }
    }

    IntMatrix hCA(m);
    IntMatrix hCB(m);
    for(int s = 0; s < m; ++s)
    {
        hCA[s] = IntArray(n);
        hCB[s] = IntArray(n);
        int adA = 0;
        int adB = 0;
        for(int i = 0; i < n; ++i)
        {
            if(i > 0)
            {
                if(cn.count(s) == 0)
                {
                    int mod = std::min(n, k);
                    int a = std::min(targetA[s][i % mod], cmax);
                    int b = std::min(targetB[s][i % mod], cmax);

                    if(ampdel)
                    {
                        adA = adA == 0 && a > base ? 1 : adA;
                        adA = adA == 0 && a < base ? -1 : adA;
                        adB = adB == 0 && b > base ? 1 : adB;
                        adB = adB == 0 && b < base ? -1 : adB;

                        a = adA >= 0 ? std::max(a, base) : std::min(a, base);
                        b = adB >= 0 ? std::max(b, base) : std::min(b, base);

                        if(a + b > cmax)
                        {
                            a = cmax - base;
                            b = cmax - a;
                        }
                    }

                    if(a + b <= cmax)
                    {
                        hCA[s][i] = a;
                        hCB[s][i] = b;
                    } else {
                        hCA[s][i] = a;
                        hCB[s][i] = std::max(cmax - a, 0);
                    }
                    assert(hCA[s][i] + hCB[s][i] <= cmax);
                    assert(hCA[s][i] >= 0);
                    assert(hCB[s][i] >= 0);
                } else {
                    hCA[s][i] = cn.at(s).first;
                    hCB[s][i] = cn.at(s).second;
                }
            } else {
                hCA[s][0] = 1;
                hCB[s][0] = 1;
            }
        }
    }

    return std::make_pair(hCA, hCB);
}


void IlpSubset::buildVariables()
{
    char buf[1024];

    _yA = VarMatrix(_m);
    _yB = VarMatrix(_m);
    for(int s = 0; s < _m; ++s)
    {
        _yA[s] = VarArray(_k);
        _yB[s] = VarArray(_k);
        for(int p = 0; p < _k; ++p)
        {
            snprintf(buf, 1024, "yA_%d_%d", s+1, p+1);
            _yA[s][p] = _model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, buf);

            snprintf(buf, 1024, "yB_%d_%d", s+1, p+1);
            _yB[s][p] = _model.addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, buf);
        }
    }

    _fA = VarMatrix(_m);
    _fB = VarMatrix(_m);
    for(int s = 0; s < _m; ++s)
    {
        int upper = 0;
        if(_cn.count(s) > 0)
        {
            upper = std::max(_cn.at(s).first + _cn.at(s).second, _cmax);
        } else {
            upper = _cmax;
        }

        _fA[s] = VarArray(_k);
        _fB[s] = VarArray(_k);
        for(int p = 0; p < _k; ++p)
        {
            snprintf(buf, 1024, "fA_%d_%d", s+1, p+1);
            _fA[s][p] = _model.addVar(0.0, upper, 0, GRB_CONTINUOUS, buf);

            snprintf(buf, 1024, "fB_%d_%d", s+1, p+1);
            _fB[s][p] = _model.addVar(0.0, upper, 0, GRB_CONTINUOUS, buf);
        }
    }

    if(_mode == MODE_t::FULL || _mode == MODE_t::CARCH)
    {
        _cA = VarMatrix(_m);
        _cB = VarMatrix(_m);
        for(int s = 0; s < _m; ++s)
        {
            _cA[s] = VarArray(_n);
            _cB[s] = VarArray(_n);
            for(int i = 0; i < _n; ++i)
            {
                snprintf(buf, 1024, "cA_%d_%d", s+1, i+1);
                _cA[s][i] = _model.addVar(_lbA[s][i], _ubA[s][i], 0, GRB_INTEGER, buf);

                snprintf(buf, 1024, "cB_%d_%d", s+1, i+1);
                _cB[s][i] = _model.addVar(_lbB[s][i], _ubB[s][i], 0, GRB_INTEGER, buf);
            }
        }

        if(_ampdel)
        {
            _adA = VarArray(_m);
            _adB = VarArray(_m);
            for(int s = 0; s < _m; ++s)
            {
                if(_cn.count(s) == 0)
                {
                    snprintf(buf, 1024, "adA_%d", s+1);
                    _adA[s] = _model.addVar(0, 1, 0, GRB_BINARY, buf);

                    snprintf(buf, 1024, "adB_%d", s+1);
                    _adB[s] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
                }
            }
        }
    }

    if(_mode == MODE_t::FULL || (_d > 0 && _mode == MODE_t::CARCH))
    {
        _bitcA = Var3Matrix(_M);
        _bitcB = Var3Matrix(_M);
        for(int b = 0; b < _M; ++b)
        {
            _bitcA[b] = VarMatrix(_m);
            _bitcB[b] = VarMatrix(_m);
            for(int s = 0; s < _m; ++s)
            {
                _bitcA[b][s] = VarArray(_n);
                _bitcB[b][s] = VarArray(_n);
                for(int i = 0; i < _n; ++i)
                {
                    snprintf(buf, 1024, "bitcA_%d_%d_%d", b+1, s+1, i+1);
                    _bitcA[b][s][i] = _model.addVar(0, 1, 0, GRB_BINARY, buf);

                    snprintf(buf, 1024, "bitcB_%d_%d_%d", b+1, s+1, i+1);
                    _bitcB[b][s][i] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
                }
            }
        }
    }

    if(_mode == MODE_t::FULL || _mode == MODE_t::UARCH)
    {
        _u = VarMatrix(_n);
        for(int i = 0; i < _n; ++i)
        {
            _u[i] = VarArray(_k);
            for(int p = 0; p < _k; ++p)
            {
                snprintf(buf, 1024, "u_%d_%d", i+1, p+1);
                _u[i][p] = _model.addVar(_lbU[i][p], _ubU[i][p], 0, GRB_CONTINUOUS, buf);
            }
        }
    }

    if(_mode == MODE_t::FULL)
    {
        _vA = Var4Matrix(_M);
        _vB = Var4Matrix(_M);
        for(int b = 0; b < _M; ++b)
        {
            _vA[b] = Var3Matrix(_m);
            _vB[b] = Var3Matrix(_m);
            for(int s = 0; s < _m; ++s)
            {
                _vA[b][s] = VarMatrix(_n);
                _vB[b][s] = VarMatrix(_n);
                for(int i = 0; i < _n; ++i)
                {
                    _vA[b][s][i] = VarArray(_k);
                    _vB[b][s][i] = VarArray(_k);
                    for(int p = 0; p < _k; ++p)
                    {
                        snprintf(buf, 1024, "vA_%d_%d_%d_%d", b+1, s+1, i+1, p+1);
                        _vA[b][s][i][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);

                        snprintf(buf, 1024, "vB_%d_%d_%d_%d", b+1, s+1, i+1, p+1);
                        _vB[b][s][i][p] = _model.addVar(0, 1, 0, GRB_CONTINUOUS, buf);
                    }
                }
            }
        }
    }

    if((_mode == MODE_t::FULL || _mode == MODE_t::UARCH) && (_mu > 0.0 || !_sampleClonal.empty()))
    {
        _x = VarMatrix(_n);
        for(int i = 0; i < _n; ++i)
        {
            _x[i] = VarArray(_k);
            for(int p = 0; p < _k; ++p)
            {
                snprintf(buf, 1024, "x_%d_%d", i+1, p+1);
                _x[i][p] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
            }
        }
    }

    _model.update();
}


void IlpSubset::buildOptionalVariables()
{
    char buf[1024];

    _z = Var3Matrix(_m);
    for(int s = 0; s < _m; ++s)
    {
        _z[s] = VarMatrix(_n);
        for(int i = 1; i < _n; ++i)
        {
            _z[s][i] = VarArray(_d);
            for(int l = 0; l < _d; ++l)
            {
                snprintf(buf, 1024, "z_%d_%d_%d", s+1, i+1, l+1);
                _z[s][i][l] = _model.addVar(0, 1, 0, GRB_BINARY, buf);
            }
        }
    }

    _model.update();
}


void IlpSubset::buildConstraints()
{
    for(int s = 0; s < _m; ++s)
    {
        for(int p = 0; p < _k; ++p)
        {
            _model.addConstr(_FA[s][p] - _fA[s][p] <= _yA[s][p]);
            _model.addConstr(_fA[s][p] - _FA[s][p] <= _yA[s][p]);

            _model.addConstr(_FB[s][p] - _fB[s][p] <= _yB[s][p]);
            _model.addConstr(_fB[s][p] - _FB[s][p] <= _yB[s][p]);
        }
    }

    if(_mode == MODE_t::FULL)
    {
        for(int s = 0; s < _m; ++s)
        {
            for(int p = 0; p < _k; ++p)
            {
                GRBLinExpr sum;
                sum = 0;
                for(int i = 0; i < _n; ++i)
                {
                    for(int b = 0; b < _M; ++b)
                    {
                        sum += _vA[b][s][i][p] * pow(2, b);
                        _model.addConstr(_vA[b][s][i][p] <= _bitcA[b][s][i]);
                        _model.addConstr(_vA[b][s][i][p] <= _u[i][p]);
                        _model.addConstr(_vA[b][s][i][p] >= _bitcA[b][s][i] + _u[i][p] - 1);
                    }
                }
                _model.addConstr(_fA[s][p] == sum);
                sum.clear();
            }
        }

        for(int s = 0; s < _m; ++s)
        {
            for(int p = 0; p < _k; ++p)
            {
                GRBLinExpr sumB;
                sumB = 0;
                for(int i = 0; i < _n; ++i)
                {
                    for(int b = 0; b < _M; ++b)
                    {
                        sumB += _vB[b][s][i][p] * pow(2, b);
                        _model.addConstr(_vB[b][s][i][p] <= _bitcB[b][s][i]);
                        _model.addConstr(_vB[b][s][i][p] <= _u[i][p]);
                        _model.addConstr(_vB[b][s][i][p] >= _bitcB[b][s][i] + _u[i][p] - 1);
                    }
                }
                _model.addConstr(_fB[s][p] == sumB);
                sumB.clear();
            }
        }

        for(int i = 0; i < _n; ++i)
        {
            for(int p = 0; p < _k; ++p)
            {
                GRBLinExpr sum;
                sum = 0;
                for(int s = 0; s < _m; ++s)
                {
                    for(int b = 0; b < _M; ++b)
                    {
                        sum += _bitcA[b][s][i] + _bitcB[b][s][i];
                    }
                }
                _model.addConstr(sum >= _u[i][p]);
                sum.clear();
            }
        }

    }

    if(_mode == MODE_t::FULL || (_d > 0 && _mode == MODE_t::CARCH))
    {
        for(int s = 0; s < _m; ++s)
        {
            for(int i = 0; i < _n; ++i)
            {
                GRBLinExpr sumA;
                sumA = 0;
                for(int b = 0; b < _M; ++b)
                {
                    sumA += _bitcA[b][s][i] * pow(2, b);
                }
                GRBLinExpr sumB;
                sumB = 0;
                for (int b = 0; b < _M; ++b) {
                    sumB += _bitcB[b][s][i] * pow(2, b);
                }
                _model.addConstr(_cA[s][i] == sumA);
                _model.addConstr(_cB[s][i] == sumB);
                sumA.clear();
                sumB.clear();
            }
        }

    }

    if (_mode == MODE_t::CARCH) {
        for(int s = 0; s < _m; ++s)
        {
            for(int p = 0; p < _k; ++p)
            {
                GRBLinExpr sum;
                sum = 0;
                for(int i = 0; i < _n; ++i)
                {
                    if(i == 0 || _barU[i][p] > TOL)
                    {
                        sum += _cA[s][i] * _barU[i][p];
                    }
                }
                _model.addConstr(_fA[s][p] == sum);
                sum.clear();
            }
        }

        for(int s = 0; s < _m; ++s)
        {
            for(int p = 0; p < _k; ++p)
            {
                GRBLinExpr sumB;
                sumB = 0;
                for(int i = 0; i < _n; ++i)
                {
                    if(i == 0 || _barU[i][p] > TOL)
                    {
                        sumB += _cB[s][i] * _barU[i][p];
                    }
                }
                _model.addConstr(_fB[s][p] == sumB);
                sumB.clear();
            }
        }
    }

    if(_mode == MODE_t::FULL || _mode == MODE_t::CARCH)
    {
        // The normal clone and other fixed cells are enforced through variable bounds, so rows that
        // those bounds already imply are omitted.
        const double totalWeight = std::accumulate(_w.begin(), _w.end(), 0.0);
        for(int s = 0; s < _m; ++s)
        {
            const int upper = upperCN(s);
            const bool noHomdel = totalWeight > 0.0 && _w[s] > 0.0 && _w[s] / totalWeight >= _maxHomdelSize;
            for(int i = 0; i < _n; ++i)
            {
                if(_ubA[s][i] + _ubB[s][i] > upper)
                {
                    _model.addConstr(_cA[s][i] + _cB[s][i] <= upper);
                }
                if(i > 0 && noHomdel && _lbA[s][i] + _lbB[s][i] < 1)
                {
                    _model.addConstr(_cA[s][i] + _cB[s][i] >= 1);
                }
            }
        }

        if(_ampdel)
        {
            for(int s = 0; s < _m; ++s)
            {
                if(_cn.count(s) == 0)
                {
                    for(int i = 1; i < _n; ++i)
                    {
                        _model.addConstr(_cA[s][i] <= _cmax * _adA[s] + _base - _base * _adA[s]);
                        _model.addConstr(_cA[s][i] >= _base * _adA[s]);

                        _model.addConstr(_cB[s][i] <= _cmax * _adB[s] + _base - _base * _adB[s]);
                        _model.addConstr(_cB[s][i] >= _base * _adB[s]);
                    }
                }
            }
        }
    }

    if (_mode == MODE_t::UARCH) {
        for(int s = 0; s < _m; ++s)
        {
            for(int p = 0; p < _k; ++p)
            {
                GRBLinExpr sum;
                sum = 0;
                for(int i = 0; i < _n; ++i)
                {
                    sum += _barCA[s][i] * _u[i][p];
                }
                _model.addConstr(_fA[s][p] == sum);
                sum.clear();
            }
        }

        for(int s = 0; s < _m; ++s)
        {
            for(int p = 0; p < _k; ++p)
            {
                GRBLinExpr sumB;
                sumB = 0;
                for(int i = 0; i < _n; ++i)
                {
                    sumB += _barCB[s][i] * _u[i][p];
                }
                _model.addConstr(_fB[s][p] == sumB);
                sumB.clear();
            }
        }
    }

    if(_mode == MODE_t::FULL || _mode == MODE_t::UARCH)
    {
        for(int p = 0; p < _k; ++p)
        {
            GRBLinExpr sum;
            sum = 0;
            for(int i = 0; i < _n; ++i)
            {
                sum += _u[i][p];
            }
            _model.addConstr(sum == 1);
            sum.clear();
        }
    }

    if((_mode == MODE_t::FULL || _mode == MODE_t::UARCH) && (_mu > 0.0 || !_sampleClonal.empty()))
    {
        for(int p = 0; p < _k; ++p)
        {
            for(int i = 1; i < _n; ++i)
            {
                _model.addConstr(_x[i][p] >= _u[i][p]);
                if(_mu > 0.0)
                {
                    _model.addConstr(_u[i][p] >= _mu * _x[i][p]);
                }
            }
        }
    }

    _model.update();
}


void IlpSubset::buildOptionalConstraints()
{
    for(int s = 0; s < _m; ++s)
    {
        for(int i = 1; i < _n; ++i)
        {
            GRBLinExpr sum;
            sum = 0;
            for(int l = 0; l < _d; ++l)
            {
                sum += _z[s][i][l];
            }
            _model.addConstr(sum == 1);
        }
    }

    for(int s = 0; s < _m; ++s)
    {
        for(int b = 0; b < _M; ++b)
        {
            for(int l = 0; l < _d; ++l)
            {
                for(int i = 1; i < _n - 1; ++i)
                {
                    for(int j = i; j < _n; ++j)
                    {
                        _model.addConstr(_bitcA[b][s][i] - _bitcA[b][s][j] <= 2 - _z[s][i][l] - _z[s][j][l]);
                        _model.addConstr(_bitcA[b][s][j] - _bitcA[b][s][i] <= 2 - _z[s][i][l] - _z[s][j][l]);

                        _model.addConstr(_bitcB[b][s][i] - _bitcB[b][s][j] <= 2 - _z[s][i][l] - _z[s][j][l]);
                        _model.addConstr(_bitcB[b][s][j] - _bitcB[b][s][i] <= 2 - _z[s][i][l] - _z[s][j][l]);
                    }
                }
            }
        }
    }

    for(int s = 0; s < _m; ++s)
    {
        for(int l = 0; l < _d - 1; ++l)
        {
            GRBLinExpr sumL;
            GRBLinExpr sumL1;
            sumL = 0;
            sumL1 = 0;
            for(int i = 1; i < _n; ++i)
            {
                sumL += _z[s][i][l] * symmCoeff(i);
                sumL1 += _z[s][i][l+1] * symmCoeff(i);
            }
            _model.addConstr(sumL <= sumL1);
            sumL.clear();
            sumL1.clear();
        }
    }

    _model.update();
}


void IlpSubset::computeFixedBounds()
{
    _conflicts.clear();
    const bool copyNumbers = _mode == MODE_t::FULL || _mode == MODE_t::CARCH;
    const bool proportions = _mode == MODE_t::FULL || _mode == MODE_t::UARCH;

    const auto fixCN = [this](const FIXED_MATRIX_t matrix, const int s, const int i, const int value)
    {
        IntMatrix& lower = matrix == MATRIX_A ? _lbA : _lbB;
        IntMatrix& upper = matrix == MATRIX_A ? _ubA : _ubB;
        if(!fixBounds(lower[s][i], upper[s][i], value))
        {
            _conflicts.push_back(FixedValue{matrix, s, i, static_cast<double>(value)});
        }
    };
    const auto fixProportion = [this](const int i, const int p, const double value)
    {
        if(!fixBounds(_lbU[i][p], _ubU[i][p], value))
        {
            _conflicts.push_back(FixedValue{MATRIX_U, i, p, value});
        }
    };

    if(copyNumbers)
    {
        _lbA = IntMatrix(_m, IntArray(_n, 0));
        _lbB = _lbA;
        _ubA = IntMatrix(_m);
        for(int s = 0; s < _m; ++s)
        {
            _ubA[s] = IntArray(_n, upperCN(s));
        }
        _ubB = _ubA;

        for(int s = 0; s < _m; ++s)
        {
            fixCN(MATRIX_A, s, 0, 1);
            fixCN(MATRIX_B, s, 0, 1);
        }
        for(auto const& cluster : _cn)
        {
            for(int i = 1; i < _n; ++i)
            {
                fixCN(MATRIX_A, cluster.first, i, cluster.second.first);
                fixCN(MATRIX_B, cluster.first, i, cluster.second.second);
            }
        }
    }

    if(proportions)
    {
        _lbU = DoubleMatrix(_n, DoubleArray(_k, 0.0));
        _ubU = DoubleMatrix(_n, DoubleArray(_k, 1.0));
    }

    const std::vector<FixedMatrixEntry>& fixedEntries = _fixedMatrix.entries();
    for(std::vector<FixedMatrixEntry>::const_iterator it = fixedEntries.begin(); it != fixedEntries.end(); ++it)
    {
        if((it->matrix == MATRIX_A || it->matrix == MATRIX_B) && copyNumbers)
        {
            fixCN(it->matrix, it->row, it->column, static_cast<int>(std::lround(it->value)));
        }
        else if(it->matrix == MATRIX_U && proportions)
        {
            fixProportion(it->row, it->column, it->value);
        }
    }

    if(_sampleClonal.empty())
    {
        return;
    }

    if(proportions)
    {
        for(int sample = 0; sample < _k; ++sample)
        {
            fixProportion(0, sample, 1.0 - _purities[sample]);
        }
    }

    const std::vector<SampleClonalEntry>& entries = _sampleClonal.entries();
    if(_mode == MODE_t::CARCH)
    {
        for(std::vector<SampleClonalEntry>::const_iterator entry = entries.begin(); entry != entries.end(); ++entry)
        {
            for(int clone = 1; clone < _n; ++clone)
            {
                if(_barU[clone][entry->sample] > TOL)
                {
                    fixCN(MATRIX_A, entry->cluster, clone, entry->cnA);
                    fixCN(MATRIX_B, entry->cluster, clone, entry->cnB);
                }
            }
        }
    }
    else if(_mode == MODE_t::UARCH)
    {
        for(int sample = 0; sample < _k; ++sample)
        {
            const std::vector<SampleClonalEntry> sampleEntries = _sampleClonal.entriesForSample(sample);
            for(int clone = 1; clone < _n; ++clone)
            {
                bool compatible = true;
                for(std::vector<SampleClonalEntry>::const_iterator entry = sampleEntries.begin();
                    entry != sampleEntries.end(); ++entry)
                {
                    compatible = compatible &&
                                 _barCA[entry->cluster][clone] == entry->cnA &&
                                 _barCB[entry->cluster][clone] == entry->cnB;
                }
                if(!compatible)
                {
                    fixProportion(clone, sample, 0.0);
                }
            }
        }
    }
}


void IlpSubset::addFixedValueConflicts()
{
    // Equalities for values that disagree with an earlier bound make the model infeasible, as before
    for(std::vector<FixedValue>::const_iterator it = _conflicts.begin(); it != _conflicts.end(); ++it)
    {
        if(it->matrix == MATRIX_A)
        {
            _model.addConstr(_cA[it->row][it->column] == it->value);
        }
        else if(it->matrix == MATRIX_B)
        {
            _model.addConstr(_cB[it->row][it->column] == it->value);
        }
        else
        {
            _model.addConstr(_u[it->row][it->column] == it->value);
        }
    }
}


void IlpSubset::fixSampleClonal()
{
    // CARCH states, UARCH exclusions, and derived purities are variable bounds (see computeFixedBounds)
    if(_sampleClonal.empty() || _mode != MODE_t::FULL)
    {
        return;
    }

    const std::vector<SampleClonalEntry>& entries = _sampleClonal.entries();
    for(std::vector<SampleClonalEntry>::const_iterator entry = entries.begin(); entry != entries.end(); ++entry)
    {
        for(int clone = 1; clone < _n; ++clone)
        {
            _model.addConstr(_cA[entry->cluster][clone] - entry->cnA <=
                             _cmax * (1.0 - _x[clone][entry->sample]));
            _model.addConstr(entry->cnA - _cA[entry->cluster][clone] <=
                             _cmax * (1.0 - _x[clone][entry->sample]));
            _model.addConstr(_cB[entry->cluster][clone] - entry->cnB <=
                             _cmax * (1.0 - _x[clone][entry->sample]));
            _model.addConstr(entry->cnB - _cB[entry->cluster][clone] <=
                             _cmax * (1.0 - _x[clone][entry->sample]));
        }
    }
}

void IlpSubset::buildObjective()
{
    GRBLinExpr obj;
    obj = 0;

    double norm = std::accumulate(_w.begin(), _w.end(), 0.0);
    DoubleArray ws;
    if(_mode == MODE_t::FULL) {log("Cluster objective weights:  ", VERBOSITY_t::DEBUGGING, _v);}
    for(auto& w : _w)
    {
        ws.push_back((double)w * 100.0 / (double)norm);
        if(_mode == MODE_t::FULL) {log(std::to_string(w) + "=  " + std::to_string(ws.back()) + ", ", VERBOSITY_t::DEBUGGING, _v);}
    }
    if(_mode == MODE_t::FULL) {log("\n", VERBOSITY_t::DEBUGGING, _v);}

    for(int s = 0; s < _m; ++s)
    {
        for(int p = 0; p < _k; ++p)
        {
            obj += _yA[s][p] * ws[s] + _yB[s][p] * ws[s];
        }
    }

    _model.setObjective(obj, GRB_MINIMIZE);
    _model.update();
}

/// Based on (Balas and Jeroslow, J. of App. Math., 1972)
void IlpSubset::buildNext()
{
    GRBLinExpr sum1;
    GRBLinExpr sum0;
    sum1 = 0;
    sum0 = 0;
    int count1 = 0;

    for(int s = 0; s < _m; ++s)
    {
        for(int i = 0; i < _n; ++i)
        {
            for(int b = 0; b < _M; ++b)
            {
                if(round(_bitcA[b][s][i].get(GRB_DoubleAttr_X)) == 1)
                {
                    sum1 += _bitcA[b][s][i];
                    ++count1;
                } else {
                    sum0 += _bitcA[b][s][i];
                }
            }

            for(int b = 0; b < _M; ++b)
            {
                if(round(_bitcB[b][s][i].get(GRB_DoubleAttr_X)) == 1)
                {
                    sum1 += _bitcB[b][s][i];
                    ++count1;
                } else {
                    sum0 += _bitcB[b][s][i];
                }
            }
        }
    }

    _model.addConstr(sum1 - sum0 <= count1 - 1);
}


void IlpSubset::test()
{
//    _model.addConstr(_cA[0][0] == 1);
//    _model.addConstr(_cA[1][0] == 1);
//    _model.addConstr(_cA[2][0] == 1);
//    _model.addConstr(_cA[3][0] == 1);
//    _model.addConstr(_cA[4][0] == 1);
//
//    _model.addConstr(_cB[0][0] == 1);
//    _model.addConstr(_cB[1][0] == 1);
//    _model.addConstr(_cB[2][0] == 1);
//    _model.addConstr(_cB[3][0] == 1);
//    _model.addConstr(_cB[4][0] == 1);

    _model.computeIIS();
    _model.write("iis.ilp");
}


IntMatrix IlpSubset::getCA()
{
    if(_mode == MODE_t::FULL || _mode == MODE_t::CARCH)
    {
        IntMatrix CA(_m);
        for(int s = 0; s < _m; ++s)
        {
            CA[s] = IntArray(_n);
            for(int i = 0; i < _n; ++i)
            {
                CA[s][i] = round(_cA[s][i].get(GRB_DoubleAttr_X));
            }
        }
        return CA;
    } else {
        return _barCA;
    }
}


IntMatrix IlpSubset::getCB()
{
    if(_mode == MODE_t::FULL || _mode == MODE_t::CARCH)
    {
        IntMatrix CB(_m);
        for(int s = 0; s < _m; ++s)
        {
            CB[s] = IntArray(_n);
            for(int i = 0; i < _n; ++i)
            {
                CB[s][i] = round(_cB[s][i].get(GRB_DoubleAttr_X));
            }
        }
        return CB;
    } else {
        return _barCB;
    }
}


DoubleMatrix IlpSubset::getU()
{
    if(_mode == MODE_t::FULL || _mode == MODE_t::UARCH)
    {
        DoubleMatrix U(_n);
        for(int i = 0; i < _n; ++i)
        {
            U[i] = DoubleArray(_k);
            for(int p = 0; p < _k; ++p)
            {
                U[i][p] = _u[i][p].get(GRB_DoubleAttr_X);
            }
        }
        return U;
    } else {
        return _barU;
    }
}


int IlpSubset::upperCN(const int s) const
{
    if(_cn.count(s) > 0)
    {
        return std::max(_cn.at(s).first + _cn.at(s).second, _cmax);
    }
    return _cmax;
}


double IlpSubset::symmCoeff(const int index)
{
    return std::pow(index+1, 1);
}
