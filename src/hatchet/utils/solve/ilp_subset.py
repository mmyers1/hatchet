import textwrap
import math
import numpy as np
from pyomo import environ as pe
from pyomo.opt import SolverStatus, TerminationCondition

from hatchet.utils.solve.utils import Random


class ILPSubset:
    def __init__(
        self,
        n,
        cn_max,
        d,
        mu,
        ampdel,
        copy_numbers,
        f_a,
        f_b,
        w,
        purities,
        fixed_matrix=None,
        base=None,
        sample_clonal=None,
    ):

        # Each ILPSubset maintains its own data, so make a deep-copy of passed-in DataFrames
        f_a, f_b = f_a.copy(deep=True), f_b.copy(deep=True)

        assert f_a.shape == f_b.shape
        assert np.all(f_a.index == f_b.index)
        assert np.all(f_a.columns == f_b.columns)

        self.m, self.k = f_a.shape
        self.f_a = f_a
        self.f_b = f_b
        self.cluster_ids = f_a.index
        self.sample_ids = f_a.columns

        self.n = n
        self.cn_max = cn_max
        self.d = d
        self.mu = mu
        self.ampdel = ampdel
        self.copy_numbers = copy_numbers
        self.w = w
        self.purities = purities
        self.fixed_matrix = fixed_matrix
        self._base = base
        self.sample_clonal = sample_clonal

        self.tol = 0.001

        self.mode = 'FULL'

        # Values we want to optimize for, as dataframes
        self.cA = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self.cB = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self.u = [[np.nan for _ in range(self.k)] for _ in range(self.n)]

        # Fixed values of cA/cB/u
        self._fixed_cA = [[np.nan for _ in range(n)] for _ in range(self.m)]
        self._fixed_cB = [[np.nan for _ in range(n)] for _ in range(self.m)]
        self._fixed_u = [[np.nan for _ in range(self.k)] for _ in range(self.n)]

        self.warmstart = False  # set on hot_start()
        self.model = None  # initialized on create_model()

    def __copy__(self):
        return ILPSubset(
            n=self.n,
            cn_max=self.cn_max,
            d=self.d,
            mu=self.mu,
            ampdel=self.ampdel,
            copy_numbers=self.copy_numbers,
            f_a=self.f_a,
            f_b=self.f_b,
            w=self.w,
            purities=self.purities,
            fixed_matrix=self.fixed_matrix,
            base=self._base,
            sample_clonal=self.sample_clonal,
        )

    def __str__(self):
        # Pyomo pprint gives us too much information - too unwieldy for large models
        # This method is implemented to supply the bare-minimum but useful model information.
        if self.model is None:
            return ''
        else:
            return textwrap.dedent(
                f"""
                # ------------------------------------------
                #   Problem Information
                # ------------------------------------------
                #     Number of constraints: {self.model.nconstraints()}
                #     Number of variables: {self.model.nvariables()}
                # ------------------------------------------
            """
            )

    @property
    def M(self):
        return math.floor(math.log2(self.cn_max)) + 1

    @property
    def base(self):
        return self._base if self._base is not None else min(2, len(self.copy_numbers))

    @staticmethod
    def symmCoeff(i):
        return math.pow(i + 1, 1)

    @property
    def optimized_cA(self):
        if self.mode == 'UARCH':
            return self._fixed_cA
        else:
            return self.cA

    @property
    def optimized_cB(self):
        if self.mode == 'UARCH':
            return self._fixed_cB
        else:
            return self.cB

    @property
    def optimized_u(self):
        if self.mode == 'CARCH':
            return self._fixed_u
        else:
            return self.u

    def create_model(self, pprint=False):

        m, n, k = self.m, self.n, self.k
        f_a, f_b = self.f_a, self.f_b
        cn_max = self.cn_max
        ampdel = self.ampdel
        copy_numbers = self.copy_numbers
        mode_t = self.mode
        d = self.d
        _M = self.M
        _base = self.base
        purities = self.purities

        model = pe.ConcreteModel()

        fA = {}
        fB = {}
        yA = {}
        yB = {}
        adA = {}
        adB = {}

        for _m in range(m):

            cluster_id = f_a.index[_m]

            # upper bound for solver
            ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

            for _k in range(k):
                yA[(_m, _k)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                model.add_component(f'yA_{_m + 1}_{_k + 1}', yA[(_m, _k)])
                yB[(_m, _k)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                model.add_component(f'yB_{_m + 1}_{_k + 1}', yB[(_m, _k)])
                fA[(_m, _k)] = pe.Var(bounds=(0, ub), domain=pe.Reals)
                model.add_component(f'fA_{_m + 1}_{_k + 1}', fA[(_m, _k)])
                fB[(_m, _k)] = pe.Var(bounds=(0, ub), domain=pe.Reals)
                model.add_component(f'fB_{_m + 1}_{_k + 1}', fB[(_m, _k)])

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):

                cluster_id = f_a.index[_m]

                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

                for _n in range(n):
                    self.cA[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model.add_component(f'cA_{_m + 1}_{_n + 1}', self.cA[_m][_n])
                    self.cB[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model.add_component(f'cB_{_m + 1}_{_n + 1}', self.cB[_m][_n])

            if ampdel:
                for _m in range(m):
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        adA[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'adA_{_m + 1}', adA[_m])
                        adB[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'adB_{_m + 1}', adB[_m])

        bitcA = {}
        bitcB = {}
        if (mode_t == 'FULL') or (d > 0 and mode_t == 'CARCH'):
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        bitcA[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(
                            f'bitcA_{_b + 1}_{_m + 1}_{_n + 1}',
                            bitcA[(_b, _m, _n)],
                        )
                        bitcB[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(
                            f'bitcB_{_b + 1}_{_m + 1}_{_n + 1}',
                            bitcB[(_b, _m, _n)],
                        )

        if mode_t in ('FULL', 'UARCH'):
            for _n in range(n):
                for _k in range(k):
                    self.u[_n][_k] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                    model.add_component(f'u_{_n + 1}_{_k + 1}', self.u[_n][_k])

        vA = {}
        vB = {}
        if mode_t == 'FULL':
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        for _k in range(k):
                            vA[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model.add_component(
                                f'vA_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}',
                                vA[(_b, _m, _n, _k)],
                            )
                            vB[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model.add_component(
                                f'vB_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}',
                                vB[(_b, _m, _n, _k)],
                            )

        x = {}
        if (mode_t in ('FULL', 'UARCH')) and (self.mu > 0 or self.sample_clonal is not None):
            for _n in range(n):
                for _k in range(k):
                    x[(_n, _k)] = pe.Var(domain=pe.Binary)
                    model.add_component(f'x_{_n + 1}_{_k + 1}', x[(_n, _k)])

        # buildOptionalVariables
        z = {}
        if (mode_t in ('FULL', 'CARCH')) and d > 0:
            for _m in range(self.m):
                for _n in range(1, self.n):
                    for _d in range(d):
                        z[(_m, _n, _d)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'z_{_m + 1}_{_n + 1}_{_d + 1}', z[(_m, _n, _d)])

        # CONSTRAINTS
        model.constraints = pe.ConstraintList()

        for _m in range(m):

            cluster_id = f_a.index[_m]
            f_a_values = f_a.loc[cluster_id].values
            f_b_values = f_b.loc[cluster_id].values

            for _k in range(k):
                _yA, _yB = yA[(_m, _k)], yB[(_m, _k)]
                _fA, _fB = fA[(_m, _k)], fB[(_m, _k)]

                model.constraints.add(float(f_a_values[_k]) - _fA <= _yA)
                model.constraints.add(_fA - float(f_a_values[_k]) <= _yA)
                model.constraints.add(float(f_b_values[_k]) - _fB <= _yB)
                model.constraints.add(_fB - float(f_b_values[_k]) <= _yB)

        if mode_t == 'FULL':
            for _m in range(m):
                for _k in range(k):

                    sum_a = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum_a += vA[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model.constraints.add(vA[(_b, _m, _n, _k)] <= bitcA[(_b, _m, _n)])
                            model.constraints.add(vA[(_b, _m, _n, _k)] <= self.u[_n][_k])
                            model.constraints.add(vA[(_b, _m, _n, _k)] >= bitcA[(_b, _m, _n)] + self.u[_n][_k] - 1)

                    model.constraints.add(fA[(_m, _k)] == sum_a)

            for _m in range(m):
                for _k in range(k):

                    sum_b = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum_b += vB[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model.constraints.add(vB[(_b, _m, _n, _k)] <= bitcB[(_b, _m, _n)])
                            model.constraints.add(vB[(_b, _m, _n, _k)] <= self.u[_n][_k])
                            model.constraints.add(vB[(_b, _m, _n, _k)] >= bitcB[(_b, _m, _n)] + self.u[_n][_k] - 1)

                    model.constraints.add(fB[(_m, _k)] == sum_b)

            for _n in range(n):
                for _k in range(k):
                    _sum = 0
                    for _m in range(m):
                        for _b in range(_M):
                            _sum += bitcA[(_b, _m, _n)] + bitcB[(_b, _m, _n)]
                    model.constraints.add(_sum >= self.u[_n][_k])

        if (mode_t == 'FULL') or (d > 0 and mode_t == 'CARCH'):
            for _m in range(m):
                cluster_id = f_a.index[_m]
                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

                for _n in range(n):
                    sum_a = 0
                    sum_b = 0
                    for _b in range(_M):
                        sum_a += bitcA[(_b, _m, _n)] * math.pow(2, _b)
                        sum_b += bitcB[(_b, _m, _n)] * math.pow(2, _b)

                    model.constraints.add(self.cA[_m][_n] == sum_a)
                    model.constraints.add(self.cB[_m][_n] == sum_b)
                    model.constraints.add(self.cA[_m][_n] + self.cB[_m][_n] <= ub)

        #  this is where we fix the purities if user provides purity values for (tumor) samples
        # in each sample s, clone 0 (healty & normal) should have the purity 1 - purity_s
        if mode_t in ('FULL', 'UARCH') and purities:
            for i, p in enumerate(purities):
                model.constraints.add(self.u[0][i] == 1 - p)

        if mode_t == 'CARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sumA = 0
                    _sumB = 0
                    for _n in range(n):
                        if _n == 0 or self._fixed_u[_n][_k] > self.tol:
                            _sumA += self.cA[_m][_n] * self._fixed_u[_n][_k]
                            _sumB += self.cB[_m][_n] * self._fixed_u[_n][_k]
                    model.constraints.add(fA[(_m, _k)] == _sumA)
                    model.constraints.add(fB[(_m, _k)] == _sumB)

                cluster_id = f_a.index[_m]
                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)
                for _n in range(n):
                    model.constraints.add(self.cA[_m][_n] + self.cB[_m][_n] <= ub)

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):
                model.constraints.add(self.cA[_m][0] == 1)
                model.constraints.add(self.cB[_m][0] == 1)

                cluster_id = self.cluster_ids[_m]

                # if the fraction of the cluster size is larger than or equal to 0.005 in the whole dataset
                # then its integer copy number cannot have (0,0) copy number
                if self.w[cluster_id] / sum(self.w) >= 0.005:
                    for _n in range(1, n):
                        model.constraints.add(self.cA[_m][_n] + self.cB[_m][_n] >= 1)

            if ampdel:
                for _m in range(m):
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        for _n in range(1, n):
                            model.constraints.add(self.cA[_m][_n] <= cn_max * adA[_m] + _base - _base * adA[_m])
                            model.constraints.add(self.cA[_m][_n] >= _base * adA[_m])
                            model.constraints.add(self.cB[_m][_n] <= cn_max * adB[_m] + _base - _base * adB[_m])
                            model.constraints.add(self.cB[_m][_n] >= _base * adB[_m])

        if mode_t == 'UARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sumA = 0
                    _sumB = 0
                    for _n in range(n):
                        _sumA += int(self._fixed_cA[_m][_n]) * self.u[_n][_k]
                        _sumB += int(self._fixed_cB[_m][_n]) * self.u[_n][_k]
                    model.constraints.add(fA[(_m, _k)] == _sumA)
                    model.constraints.add(fB[(_m, _k)] == _sumB)

        if mode_t in ('FULL', 'UARCH'):
            for _k in range(k):
                _sum = 0
                for _n in range(n):
                    _sum += self.u[_n][_k]
                model.constraints.add(_sum == 1)

        if (mode_t in ('FULL', 'UARCH')) and (self.mu > 0 or self.sample_clonal is not None):
            for _k in range(k):
                for _n in range(1, n):
                    model.constraints.add(x[(_n, _k)] >= self.u[_n][_k])
                    if self.mu > 0:
                        model.constraints.add(self.u[_n][_k] >= self.mu * x[(_n, _k)])

        self.apply_sample_clonal_constraints(model, x)

        # buildOptionalConstraints
        if (mode_t in ('FULL', 'CARCH')) and d > 0:
            for _m in range(self.m):
                for _n in range(1, self.n):
                    _sum = 0
                    for _d in range(self.d):
                        _sum += z[(_m, _n, _d)]
                    model.constraints.add(_sum == 1)

            for _m in range(self.m):
                for _M in range(self.M):
                    for _d in range(d):
                        for _i in range(1, self.n - 1):
                            for _j in range(1, self.n):
                                model.constraints.add(
                                    bitcA[(_M, _m, _i)] - bitcA[(_M, _m, _j)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcA[(_M, _m, _j)] - bitcA[(_M, _m, _i)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcB[(_M, _m, _i)] - bitcB[(_M, _m, _j)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcB[(_M, _m, _j)] - bitcB[(_M, _m, _i)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )

            for _m in range(self.m):
                for _d in range(d - 1):
                    _sum_l = _sum_l1 = 0
                    for _n in range(1, self.n):
                        _sum_l += z[(_m, _n, _d)] * self.symmCoeff(_n)
                        _sum_l1 += z[(_m, _n, _d + 1)] * self.symmCoeff(_n)
                        model.constraints.add(_sum_l <= _sum_l1)

        if mode_t in ('FULL', 'CARCH'):
            self.build_symmetry_breaking(model)
            self.fix_given_cn(model)

        self.apply_fixed_matrix_constraints(model)

        if mode_t == 'FULL':
            self.hot_start()

        obj = 0
        for _m in range(m):
            for _k in range(k):
                cluster_id = self.cluster_ids[_m]
                obj += (yA[(_m, _k)] + yB[(_m, _k)]) * self.w[cluster_id]

        model.obj = pe.Objective(expr=obj, sense=pe.minimize)
        self.model = model

        if pprint:
            print(str(self))

    def build_symmetry_breaking(self, model):
        for i in range(1, self.n - 1):
            _sum1 = 0
            _sum2 = 0
            for _m in range(self.m):
                _symCoeff = self.symmCoeff(_m)
                _sum1 += self.cA[_m][i] * _symCoeff + self.cB[_m][i] * _symCoeff
                _sum2 += self.cA[_m][i + 1] * _symCoeff + self.cB[_m][i + 1] * _symCoeff
            model.constraints.add(_sum1 <= _sum2)

    def apply_fixed_matrix_constraints(self, model):
        if self.fixed_matrix is None:
            return

        model.fixed_matrix_constraints = pe.ConstraintList()
        for entry in self.fixed_matrix.entries:
            if entry.matrix == 'A' and self.mode in ('FULL', 'CARCH'):
                model.fixed_matrix_constraints.add(self.cA[entry.row_index][entry.column_index] == entry.value)
            elif entry.matrix == 'B' and self.mode in ('FULL', 'CARCH'):
                model.fixed_matrix_constraints.add(self.cB[entry.row_index][entry.column_index] == entry.value)
            elif entry.matrix == 'U' and self.mode in ('FULL', 'UARCH'):
                model.fixed_matrix_constraints.add(self.u[entry.row_index][entry.column_index] == entry.value)

    def apply_sample_clonal_constraints(self, model, presence):
        if self.sample_clonal is None:
            return

        model.sample_clonal_constraints = pe.ConstraintList()
        if self.mode == 'FULL':
            for entry in self.sample_clonal.entries:
                for clone_index in range(1, self.n):
                    clone_present = presence[(clone_index, entry.sample_index)]
                    model.sample_clonal_constraints.add(
                        self.cA[entry.cluster_index][clone_index] - entry.cn_a <= self.cn_max * (1 - clone_present)
                    )
                    model.sample_clonal_constraints.add(
                        entry.cn_a - self.cA[entry.cluster_index][clone_index] <= self.cn_max * (1 - clone_present)
                    )
                    model.sample_clonal_constraints.add(
                        self.cB[entry.cluster_index][clone_index] - entry.cn_b <= self.cn_max * (1 - clone_present)
                    )
                    model.sample_clonal_constraints.add(
                        entry.cn_b - self.cB[entry.cluster_index][clone_index] <= self.cn_max * (1 - clone_present)
                    )
        elif self.mode == 'CARCH':
            for entry in self.sample_clonal.entries:
                for clone_index in range(1, self.n):
                    if self._fixed_u[clone_index][entry.sample_index] > self.tol:
                        model.sample_clonal_constraints.add(self.cA[entry.cluster_index][clone_index] == entry.cn_a)
                        model.sample_clonal_constraints.add(self.cB[entry.cluster_index][clone_index] == entry.cn_b)
        elif self.mode == 'UARCH':
            for sample_index in range(self.k):
                entries = self.sample_clonal.entries_for_sample(sample_index)
                for clone_index in range(1, self.n):
                    if any(
                        int(self._fixed_cA[entry.cluster_index][clone_index]) != entry.cn_a
                        or int(self._fixed_cB[entry.cluster_index][clone_index]) != entry.cn_b
                        for entry in entries
                    ):
                        model.sample_clonal_constraints.add(self.u[clone_index][sample_index] == 0)

    def overlay_fixed_matrix(self, matrix, values):
        if self.fixed_matrix is None:
            return
        for entry in self.fixed_matrix.entries_for(matrix):
            values[entry.row_index][entry.column_index] = entry.value

    def fix_given_cn(self, model):
        for _m in range(self.m):
            cluster_id = self.f_a.index[_m]
            if cluster_id in self.copy_numbers:
                _cnA, _cnB = self.copy_numbers[cluster_id]
                for _n in range(1, self.n):
                    model.constraints.add(self.cA[_m][_n] == _cnA)
                    model.constraints.add(self.cB[_m][_n] == _cnB)

    def first_hot_start(self):

        if self.d > 0:
            targetA = np.empty((self.m, self.d))
            targetB = np.empty_like(targetA)
            for _m in range(self.m):
                targetA[_m, :] = np.random.choice(self.f_a.iloc[_m].values, self.d, replace=False)
                targetB[_m, :] = np.random.choice(self.f_b.iloc[_m].values, self.d, replace=False)
            targetA = targetA.round()
            targetB = targetB.round()
        else:
            targetA = np.round(self.f_a.values)
            targetB = np.round(self.f_b.values)

        hcA = np.zeros((self.m, self.n))
        hcB = np.zeros((self.m, self.n))

        for _m in range(self.m):
            cluster_id = self.f_a.index[_m]

            adA = 0
            adB = 0
            for _n in range(self.n):
                if _n == 0:
                    hcA[_m][0] = 1
                    hcB[_m][0] = 1
                else:
                    if cluster_id in self.copy_numbers:
                        hcA[_m][_n] = self.copy_numbers[cluster_id][0]
                        hcB[_m][_n] = self.copy_numbers[cluster_id][1]
                    else:
                        mod = min(self.n, self.k)
                        a = min(targetA[_m][_n % mod], self.cn_max)
                        b = min(targetB[_m][_n % mod], self.cn_max)

                        if self.ampdel:
                            base = self.base
                            if adA == 0 and a > base:
                                adA = 1
                            if adA == 0 and a < base:
                                adA = -1
                            if adB == 0 and b > base:
                                adB = 1
                            if adB == 0 and b < base:
                                adB = -1

                            a = max(a, base) if adA >= 0 else min(a, base)
                            b = max(b, base) if adB >= 0 else min(b, base)

                            if a + b > self.cn_max:
                                a = self.cn_max - base
                                b = self.cn_max - a

                        if a + b <= self.cn_max:
                            hcA[_m][_n] = a
                            hcB[_m][_n] = b
                        else:
                            hcA[_m][_n] = a
                            hcB[_m][_n] = max(self.cn_max - a, 0)

                        assert hcA[_m][_n] + hcB[_m][_n] <= self.cn_max
                        assert hcA[_m][_n] >= 0
                        assert hcB[_m][_n] >= 0

                self.overlay_fixed_matrix('A', hcA)
                self.overlay_fixed_matrix('B', hcB)
        return hcA, hcB

    def hot_start(self, f_a=None, f_b=None):
        def get_rank(hcA, hcB):
            # m, n = hcA.shape
            m, n = len(hcA), len(hcA[0])
            rank = np.zeros(n)
            rank[0] = -1
            for _m in range(m):
                for _n in range(1, n):
                    rank[_n] += hcA[_m][_n] * self.symmCoeff(_m) + hcB[_m][_n] * self.symmCoeff(_m)

            return rank

        if f_a is None and f_b is None:
            f_a, f_b = self.first_hot_start()
        if self.fixed_matrix is None:
            rank = get_rank(f_a, f_b)
            rank_indices = np.argsort(rank)
        else:
            rank_indices = np.arange(self.n)

        for _m in range(self.m):
            for _n in range(0, self.n):
                self.cA[_m][rank_indices[_n]].value = f_a[_m][_n]
                self.cB[_m][rank_indices[_n]].value = f_b[_m][_n]
        self.warmstart = True

    def fix_u(self, u):
        if len(u) != self.n or any(len(row) != self.k for row in u):
            raise ValueError(f'Fixed U must have shape ({self.n}, {self.k})')
        self._fixed_u = [list(row) for row in u]
        self.overlay_fixed_matrix('U', self._fixed_u)
        for sample_index in range(self.k):
            sample_sum = sum(self._fixed_u[clone_index][sample_index] for clone_index in range(self.n))
            if abs(sample_sum - 1) > self.tol:
                raise ValueError(f'Fixed U column {sample_index} sums to {sample_sum}, not 1')
            for clone_index in range(1, self.n):
                value = self._fixed_u[clone_index][sample_index]
                if self.tol < value < self.mu - self.tol:
                    raise ValueError(
                        f'Fixed U[{clone_index}, {sample_index}]={value} is between 0 and minprop {self.mu}'
                    )
        self.mode = 'CARCH'
        # TODO: sanity checks in fixU

    def fix_c(self, cA, cB):
        if len(cA) != self.m or any(len(row) != self.n for row in cA):
            raise ValueError(f'Fixed A must have shape ({self.m}, {self.n})')
        if len(cB) != self.m or any(len(row) != self.n for row in cB):
            raise ValueError(f'Fixed B must have shape ({self.m}, {self.n})')
        self._fixed_cA = [list(row) for row in cA]
        self._fixed_cB = [list(row) for row in cB]
        self.overlay_fixed_matrix('A', self._fixed_cA)
        self.overlay_fixed_matrix('B', self._fixed_cB)
        self.mode = 'UARCH'
        # TODO: sanity checks in fixC

    def build_random_u(self, random_seed=None):
        if self.sample_clonal is not None:
            with Random(random_seed):
                return self.build_sample_clonal_u()

        def _calculate_size_bubbles(mu):
            if mu <= 0.1:
                return 10
            elif mu <= 0.15:
                return 6
            elif mu <= 0.2:
                return 5
            else:
                return 3

        def _build_partition_vector(n, n_parts, size_bubbles, mu=0.03):
            """
            Return fractional components of n components as an ndarray of size n
            """
            # initialize all fractions to 0
            result = np.zeros(n)
            # generate n_parts random positions from all available positions
            positions = np.random.choice(np.arange(n), n_parts, replace=False)

            bubbles = np.sort(np.random.choice(np.arange(1, size_bubbles), n_parts - 1, replace=False)) / size_bubbles
            _result = np.diff(bubbles, prepend=0, append=1)

            # set fractions at selected positions to successive difference between bubbles
            result[positions] = _result
            # for fractions that are within tol of mu, clamp them up to mu
            result[(mu - self.tol <= result) & (result < mu)] = mu

            return result

        size_bubbles = _calculate_size_bubbles(self.mu)
        U = np.empty((self.n, self.k))

        with Random(random_seed):
            for _k in range(self.k):
                # Generate 2 random ints between [1, n] and take the max - biasing towards larger numbers
                _n0 = np.random.randint(1, self.n + 1)
                _n1 = np.random.randint(1, self.n + 1)
                n_parts = min(max(_n0, _n1), size_bubbles)  # no. of clones with non-zero proportion in the mix
                v = _build_partition_vector(self.n, n_parts, size_bubbles, mu=self.mu)
                U[:, _k] = v
        self.condition_fixed_u(U)
        return U

    def build_sample_clonal_u(self):
        if self.purities is None or len(self.purities) != self.k:
            raise ValueError('Sample-clonal seeds require one derived purity per sample')

        U = np.zeros((self.n, self.k))
        fixed_tumor_clones = {sample_index: set() for sample_index in range(self.k)}
        assigned_samples = {clone_index: set() for clone_index in range(1, self.n)}

        for sample_index, purity in enumerate(self.purities):
            U[0][sample_index] = 1 - purity

        if self.fixed_matrix is not None:
            for entry in self.fixed_matrix.entries_for('U'):
                clone_index = entry.row_index
                sample_index = entry.column_index
                if clone_index == 0:
                    if abs(entry.value - U[0][sample_index]) > self.tol:
                        raise ValueError(
                            f'Fixed normal proportion for sample {sample_index} conflicts with sample-clonal purity'
                        )
                    continue
                fixed_tumor_clones[sample_index].add(clone_index)
                U[clone_index][sample_index] = entry.value
                if entry.value > self.tol:
                    assigned_samples[clone_index].add(sample_index)

        states_by_sample = {}
        for sample_index in range(self.k):
            states_by_sample[sample_index] = {
                entry.cluster_index: (entry.cn_a, entry.cn_b)
                for entry in self.sample_clonal.entries_for_sample(sample_index)
            }

        def compatible(left, right):
            shared_clusters = set(states_by_sample[left]).intersection(states_by_sample[right])
            return all(
                states_by_sample[left][cluster] == states_by_sample[right][cluster] for cluster in shared_clusters
            )

        for clone_index, sample_indices in assigned_samples.items():
            sample_indices = tuple(sample_indices)
            if any(
                not compatible(sample_indices[left], sample_indices[right])
                for left in range(len(sample_indices))
                for right in range(left + 1, len(sample_indices))
            ):
                raise ValueError(f'Fixed U entries assign clone {clone_index} to incompatible samples')

        residual_by_sample = {}
        for sample_index, purity in enumerate(self.purities):
            fixed_tumor_sum = sum(U[clone_index][sample_index] for clone_index in range(1, self.n))
            residual = purity - fixed_tumor_sum
            if residual < -self.tol:
                raise ValueError(f'Fixed tumor proportions exceed sample-clonal purity for sample {sample_index}')
            if residual > self.tol:
                if self.mu > 0 and residual < self.mu - self.tol:
                    raise ValueError(
                        f'Remaining tumor purity {residual} for sample {sample_index} is below minprop {self.mu}'
                    )
                residual_by_sample[sample_index] = residual

        sample_order = list(residual_by_sample)
        np.random.shuffle(sample_order)
        sample_order.sort(
            key=lambda sample_index: sum(
                not compatible(sample_index, other_index) for other_index in residual_by_sample
            ),
            reverse=True,
        )

        assignment = {}

        def assign_sample(position):
            if position == len(sample_order):
                return True
            sample_index = sample_order[position]
            candidates = [
                clone_index
                for clone_index in range(1, self.n)
                if clone_index not in fixed_tumor_clones[sample_index]
                and all(compatible(sample_index, assigned) for assigned in assigned_samples[clone_index])
            ]
            np.random.shuffle(candidates)
            for clone_index in candidates:
                assigned_samples[clone_index].add(sample_index)
                assignment[sample_index] = clone_index
                if assign_sample(position + 1):
                    return True
                assigned_samples[clone_index].remove(sample_index)
                del assignment[sample_index]
            return False

        if not assign_sample(0):
            raise ValueError('No feasible tumor-clone assignment exists for the sample-clonal states')

        for sample_index, clone_index in assignment.items():
            U[clone_index][sample_index] = residual_by_sample[sample_index]
        return U

    def condition_fixed_u(self, u):
        if self.fixed_matrix is None:
            return

        for sample_index in range(self.k):
            entries = [entry for entry in self.fixed_matrix.entries_for('U') if entry.column_index == sample_index]
            if not entries:
                continue

            fixed_clones = {entry.row_index for entry in entries}
            fixed_sum = sum(entry.value for entry in entries)
            free_clones = [clone_index for clone_index in range(self.n) if clone_index not in fixed_clones]
            for clone_index in free_clones:
                u[clone_index][sample_index] = 0
            for entry in entries:
                u[entry.row_index][sample_index] = entry.value

            residual = 1 - fixed_sum
            if residual > self.tol:
                eligible_clones = free_clones
                if residual < self.mu - self.tol:
                    eligible_clones = [clone_index for clone_index in free_clones if clone_index == 0]
                if not eligible_clones:
                    raise ValueError(f'Cannot allocate remaining U proportion {residual} in sample {sample_index}')
                selected_clone = np.random.choice(eligible_clones)
                u[selected_clone][sample_index] = residual

    def run(self, solver_type='gurobi', timelimit=None, write_path=None):
        if solver_type == 'gurobipy':
            solver = pe.SolverFactory('gurobi', solver_io='python')
        else:
            solver = pe.SolverFactory(solver_type)

        kwargs = {'report_timing': False}
        if timelimit is not None:
            kwargs['timelimit'] = int(timelimit)
        if solver.warm_start_capable():
            kwargs['warmstart'] = self.warmstart

        results = solver.solve(self.model, **kwargs)
        if (
            (results.solver.status == SolverStatus.ok)
            and (
                results.solver.termination_condition
                in (
                    TerminationCondition.optimal,
                    TerminationCondition.feasible,
                )
            )
            or (
                results.solver.status == SolverStatus.aborted
                and results.solver.termination_condition == TerminationCondition.maxTimeLimit
            )
        ):
            pass
        else:
            return None

        if write_path is not None:
            self.model.write(write_path)

        return (
            self.model.obj(),
            [[int(getattr(x, 'value', x)) for x in row] for row in self.optimized_cA],
            [[int(getattr(x, 'value', x)) for x in row] for row in self.optimized_cB],
            [[getattr(x, 'value', x) for x in row] for row in self.optimized_u],
            self.cluster_ids,
            self.sample_ids,
        )


class ILPSubsetSplit(ILPSubset):
    def __init__(
        self,
        n,
        cn_max,
        d,
        mu,
        ampdel,
        copy_numbers,
        f_a,
        f_b,
        binsA,
        binsB,
        lengths,
        purities=None,
        fixed_matrix=None,
        base=None,
        sample_clonal=None,
    ):

        # Each ILPSubset maintains its own data, so make a deep-copy of passed-in DataFrames
        f_a, f_b = f_a.copy(deep=True), f_b.copy(deep=True)
        self.binsA = {k: v.copy(deep=True) for k, v in binsA.items()}
        self.binsB = {k: v.copy(deep=True) for k, v in binsB.items()}
        self.lengths = {k: np.copy(v) for k, v in lengths.items()}

        assert f_a.shape == f_b.shape
        assert np.all(f_a.index == f_b.index)
        assert np.all(f_a.columns == f_b.columns)

        self.m, self.k = f_a.shape
        self.f_a = f_a
        self.f_b = f_b
        self.cluster_ids = f_a.index
        self.sample_ids = f_a.columns

        self.n = n
        self.cn_max = cn_max
        self.d = d
        self.mu = mu
        self.ampdel = ampdel
        self.copy_numbers = copy_numbers
        self.purities = purities
        self.fixed_matrix = fixed_matrix
        self._base = base
        self.sample_clonal = sample_clonal

        self.tol = 0.001

        self.mode = 'FULL'

        # Values we want to optimize for, as dataframes
        self.cA = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self.cB = [[np.nan for _ in range(self.n)] for _ in range(self.m)]
        self.u = [[np.nan for _ in range(self.k)] for _ in range(self.n)]

        # Fixed values of cA/cB/u
        self._fixed_cA = [[np.nan for _ in range(n)] for _ in range(self.m)]
        self._fixed_cB = [[np.nan for _ in range(n)] for _ in range(self.m)]
        self._fixed_u = [[np.nan for _ in range(self.k)] for _ in range(self.n)]

        self.warmstart = False  # set on hot_start()
        self.model = None  # initialized on create_model()

    def __copy__(self):
        return ILPSubsetSplit(
            n=self.n,
            cn_max=self.cn_max,
            d=self.d,
            mu=self.mu,
            ampdel=self.ampdel,
            copy_numbers=self.copy_numbers,
            f_a=self.f_a,
            f_b=self.f_b,
            binsA=self.binsA,
            binsB=self.binsB,
            lengths=self.lengths,
            purities=self.purities,
            fixed_matrix=self.fixed_matrix,
            base=self._base,
            sample_clonal=self.sample_clonal,
        )

    def create_model(self, pprint=False):

        m, n, k = self.m, self.n, self.k
        f_a = self.f_a
        cn_max = self.cn_max
        ampdel = self.ampdel
        copy_numbers = self.copy_numbers
        mode_t = self.mode
        d = self.d
        _M = self.M
        _base = self.base
        purities = self.purities

        model = pe.ConcreteModel()

        fA = {}
        fB = {}
        yA = {}
        yB = {}
        adA = {}
        adB = {}

        for _m in range(m):

            cluster_id = f_a.index[_m]

            # upper bound for solver
            ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

            for _k in range(k):

                fA[(_m, _k)] = pe.Var(bounds=(0, ub), domain=pe.Reals)
                model.add_component(f'fA_{_m + 1}_{_k + 1}', fA[(_m, _k)])
                fB[(_m, _k)] = pe.Var(bounds=(0, ub), domain=pe.Reals)
                model.add_component(f'fB_{_m + 1}_{_k + 1}', fB[(_m, _k)])

                # Need an objective function term for each bin for absolute vlaue
                for _i in range(len(self.binsA[cluster_id])):
                    yA[(_m, _k, _i)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                    model.add_component(f'yA_{_m + 1}_{_k + 1}_{_i + 1}', yA[(_m, _k, _i)])
                    yB[(_m, _k, _i)] = pe.Var(bounds=(0, np.inf), domain=pe.Reals)
                    model.add_component(f'yB_{_m + 1}_{_k + 1}_{_i + 1}', yB[(_m, _k, _i)])

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):

                cluster_id = f_a.index[_m]

                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

                for _n in range(n):
                    self.cA[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model.add_component(f'cA_{_m + 1}_{_n + 1}', self.cA[_m][_n])
                    self.cB[_m][_n] = pe.Var(bounds=(0, ub), domain=pe.Integers)
                    model.add_component(f'cB_{_m + 1}_{_n + 1}', self.cB[_m][_n])

            if ampdel:
                for _m in range(m):
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        adA[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'adA_{_m + 1}', adA[_m])
                        adB[_m] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'adB_{_m + 1}', adB[_m])

        bitcA = {}
        bitcB = {}
        if (mode_t == 'FULL') or (d > 0 and mode_t == 'CARCH'):
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        bitcA[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(
                            f'bitcA_{_b + 1}_{_m + 1}_{_n + 1}',
                            bitcA[(_b, _m, _n)],
                        )
                        bitcB[(_b, _m, _n)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(
                            f'bitcB_{_b + 1}_{_m + 1}_{_n + 1}',
                            bitcB[(_b, _m, _n)],
                        )

        if mode_t in ('FULL', 'UARCH'):
            for _n in range(n):
                for _k in range(k):
                    self.u[_n][_k] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                    model.add_component(f'u_{_n + 1}_{_k + 1}', self.u[_n][_k])

        vA = {}
        vB = {}
        if mode_t == 'FULL':
            for _b in range(_M):
                for _m in range(m):
                    for _n in range(n):
                        for _k in range(k):
                            vA[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model.add_component(
                                f'vA_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}',
                                vA[(_b, _m, _n, _k)],
                            )
                            vB[(_b, _m, _n, _k)] = pe.Var(bounds=(0, 1), domain=pe.Reals)
                            model.add_component(
                                f'vB_{_b + 1}_{_m + 1}_{_n + 1}_{_k + 1}',
                                vB[(_b, _m, _n, _k)],
                            )

        x = {}
        if (mode_t in ('FULL', 'UARCH')) and (self.mu > 0 or self.sample_clonal is not None):
            for _n in range(n):
                for _k in range(k):
                    x[(_n, _k)] = pe.Var(domain=pe.Binary)
                    model.add_component(f'x_{_n + 1}_{_k + 1}', x[(_n, _k)])

        # buildOptionalVariables
        z = {}
        if (mode_t in ('FULL', 'CARCH')) and d > 0:
            for _m in range(self.m):
                for _n in range(1, self.n):
                    for _d in range(d):
                        z[(_m, _n, _d)] = pe.Var(bounds=(0, 1), domain=pe.Binary)
                        model.add_component(f'z_{_m + 1}_{_n + 1}_{_d + 1}', z[(_m, _n, _d)])

        # CONSTRAINTS
        model.constraints = pe.ConstraintList()

        for _m in range(m):

            cluster_id = f_a.index[_m]
            f_a_values = self.binsA[cluster_id].values
            f_b_values = self.binsB[cluster_id].values

            for _k in range(k):
                _fA, _fB = fA[(_m, _k)], fB[(_m, _k)]

                for _i in range(f_a_values.shape[0]):
                    _yA, _yB = yA[(_m, _k, _i)], yB[(_m, _k, _i)]

                    # Instead of 1 term per inferred CN, use all values in bin
                    model.constraints.add(f_a_values[_i, _k] - _fA <= _yA)
                    model.constraints.add(_fA - f_a_values[_i, _k] <= _yA)

                    # model.constraints.add(myA >= _yA)
                    model.constraints.add(f_b_values[_i, _k] - _fB <= _yB)
                    model.constraints.add(_fB - f_b_values[_i, _k] <= _yB)
                    # model.constraints.add(myB >= _yB)

        if mode_t == 'FULL':
            for _m in range(m):
                for _k in range(k):

                    sum_a = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum_a += vA[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model.constraints.add(vA[(_b, _m, _n, _k)] <= bitcA[(_b, _m, _n)])
                            model.constraints.add(vA[(_b, _m, _n, _k)] <= self.u[_n][_k])
                            model.constraints.add(vA[(_b, _m, _n, _k)] >= bitcA[(_b, _m, _n)] + self.u[_n][_k] - 1)

                    model.constraints.add(fA[(_m, _k)] == sum_a)

            for _m in range(m):
                for _k in range(k):

                    sum_b = 0
                    for _n in range(n):
                        for _b in range(_M):
                            sum_b += vB[(_b, _m, _n, _k)] * math.pow(2, _b)
                            model.constraints.add(vB[(_b, _m, _n, _k)] <= bitcB[(_b, _m, _n)])
                            model.constraints.add(vB[(_b, _m, _n, _k)] <= self.u[_n][_k])
                            model.constraints.add(vB[(_b, _m, _n, _k)] >= bitcB[(_b, _m, _n)] + self.u[_n][_k] - 1)

                    model.constraints.add(fB[(_m, _k)] == sum_b)

            for _n in range(n):
                for _k in range(k):
                    _sum = 0
                    for _m in range(m):
                        for _b in range(_M):
                            _sum += bitcA[(_b, _m, _n)] + bitcB[(_b, _m, _n)]
                    model.constraints.add(_sum >= self.u[_n][_k])

        if (mode_t == 'FULL') or (d > 0 and mode_t == 'CARCH'):
            for _m in range(m):
                cluster_id = f_a.index[_m]
                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)

                for _n in range(n):
                    sum_a = 0
                    sum_b = 0
                    for _b in range(_M):
                        sum_a += bitcA[(_b, _m, _n)] * math.pow(2, _b)
                        sum_b += bitcB[(_b, _m, _n)] * math.pow(2, _b)

                    model.constraints.add(self.cA[_m][_n] == sum_a)
                    model.constraints.add(self.cB[_m][_n] == sum_b)
                    model.constraints.add(self.cA[_m][_n] + self.cB[_m][_n] <= ub)

        if mode_t in ('FULL', 'UARCH') and purities:
            for sample_index, purity in enumerate(purities):
                model.constraints.add(self.u[0][sample_index] == 1 - purity)

        if mode_t == 'CARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sumA = 0
                    _sumB = 0
                    for _n in range(n):
                        if _n == 0 or self._fixed_u[_n][_k] > self.tol:
                            _sumA += self.cA[_m][_n] * self._fixed_u[_n][_k]
                            _sumB += self.cB[_m][_n] * self._fixed_u[_n][_k]
                    model.constraints.add(fA[(_m, _k)] == _sumA)
                    model.constraints.add(fB[(_m, _k)] == _sumB)

                cluster_id = f_a.index[_m]
                # upper bound for solver
                ub = max(sum(copy_numbers.get(cluster_id, (0, 0))), cn_max)
                for _n in range(n):
                    model.constraints.add(self.cA[_m][_n] + self.cB[_m][_n] <= ub)

        if mode_t in ('FULL', 'CARCH'):
            for _m in range(m):
                model.constraints.add(self.cA[_m][0] == 1)
                model.constraints.add(self.cB[_m][0] == 1)

            if ampdel:
                for _m in range(m):
                    cluster_id = f_a.index[_m]
                    if cluster_id not in copy_numbers:
                        for _n in range(1, n):
                            model.constraints.add(self.cA[_m][_n] <= cn_max * adA[_m] + _base - _base * adA[_m])
                            model.constraints.add(self.cA[_m][_n] >= _base * adA[_m])
                            model.constraints.add(self.cB[_m][_n] <= cn_max * adB[_m] + _base - _base * adB[_m])
                            model.constraints.add(self.cB[_m][_n] >= _base * adB[_m])

        if mode_t == 'UARCH':
            # TODO: These loops can be collapsed once validation against C++ is complete
            for _m in range(m):
                for _k in range(k):
                    _sumA = 0
                    _sumB = 0
                    for _n in range(n):
                        _sumA += int(self._fixed_cA[_m][_n]) * self.u[_n][_k]
                        _sumB += int(self._fixed_cB[_m][_n]) * self.u[_n][_k]
                    model.constraints.add(fA[(_m, _k)] == _sumA)
                    model.constraints.add(fB[(_m, _k)] == _sumB)

        if mode_t in ('FULL', 'UARCH'):
            for _k in range(k):
                _sum = 0
                for _n in range(n):
                    _sum += self.u[_n][_k]
                model.constraints.add(_sum == 1)

        if (mode_t in ('FULL', 'UARCH')) and (self.mu > 0 or self.sample_clonal is not None):
            for _k in range(k):
                for _n in range(1, n):
                    model.constraints.add(x[(_n, _k)] >= self.u[_n][_k])
                    if self.mu > 0:
                        model.constraints.add(self.u[_n][_k] >= self.mu * x[(_n, _k)])

        self.apply_sample_clonal_constraints(model, x)

        # buildOptionalConstraints
        if (mode_t in ('FULL', 'CARCH')) and d > 0:
            for _m in range(self.m):
                for _n in range(1, self.n):
                    _sum = 0
                    for _d in range(self.d):
                        _sum += z[(_m, _n, _d)]
                    model.constraints.add(_sum == 1)

            for _m in range(self.m):
                for _M in range(self.M):
                    for _d in range(d):
                        for _i in range(1, self.n - 1):
                            for _j in range(1, self.n):
                                model.constraints.add(
                                    bitcA[(_M, _m, _i)] - bitcA[(_M, _m, _j)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcA[(_M, _m, _j)] - bitcA[(_M, _m, _i)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcB[(_M, _m, _i)] - bitcB[(_M, _m, _j)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )
                                model.constraints.add(
                                    bitcB[(_M, _m, _j)] - bitcB[(_M, _m, _i)] <= 2 - z[(_m, _i, _d)] - z[(_m, _j, _d)]
                                )

            for _m in range(self.m):
                for _d in range(d - 1):
                    _sum_l = _sum_l1 = 0
                    for _n in range(1, self.n):
                        _sum_l += z[(_m, _n, _d)] * self.symmCoeff(_n)
                        _sum_l1 += z[(_m, _n, _d + 1)] * self.symmCoeff(_n)
                        model.constraints.add(_sum_l <= _sum_l1)

        if mode_t in ('FULL', 'CARCH'):
            self.build_symmetry_breaking(model)
            self.fix_given_cn(model)

        self.apply_fixed_matrix_constraints(model)

        if mode_t == 'FULL':
            self.hot_start()

        obj = 0
        for _m in range(m):
            cluster_id = self.cluster_ids[_m]
            n_bins = len(self.binsA[cluster_id])

            for _k in range(k):
                for _i in range(n_bins):
                    obj += (yA[(_m, _k, _i)] + yB[(_m, _k, _i)]) * self.lengths[cluster_id][_i]

        model.obj = pe.Objective(expr=obj, sense=pe.minimize)
        self.model = model

        if pprint:
            print(str(self))
