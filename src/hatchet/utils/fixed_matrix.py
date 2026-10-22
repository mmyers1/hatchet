import math
import os
import re
from dataclasses import dataclass
from typing import Tuple, Union


FORMAT_MARKER = '# HATCHET_FIXED_ENTRIES 1'
HEADER = 'matrix\trow\tcolumn\tvalue'
_CLONES_PATTERN = re.compile(r'^#\s+clones\s+(\d+)\s*$')
_CLONE_PATTERN = re.compile(r'^clone([1-9]\d*)$')


class FixedMatrixFormatError(ValueError):
    pass


@dataclass(frozen=True)
class FixedMatrixEntry:
    matrix: str
    row: str
    column: str
    value: Union[int, float]
    line_number: int

    @property
    def coordinate(self):
        return self.matrix, self.row, self.column


@dataclass(frozen=True)
class FixedMatrix:
    clones: int
    entries: Tuple[FixedMatrixEntry, ...]
    source: str

    def entries_for(self, matrix):
        return tuple(entry for entry in self.entries if entry.matrix == matrix)


@dataclass(frozen=True)
class ResolvedFixedMatrixEntry:
    matrix: str
    row_index: int
    column_index: int
    value: Union[int, float]
    source_entry: FixedMatrixEntry


@dataclass(frozen=True)
class ResolvedFixedMatrix:
    clones: int
    entries: Tuple[ResolvedFixedMatrixEntry, ...]
    cluster_ids: tuple
    sample_ids: tuple
    source: str

    def entries_for(self, matrix):
        return tuple(entry for entry in self.entries if entry.matrix == matrix)

    @property
    def addresses_tumor_clones(self):
        """True when any entry refers to a tumor clone, which pins that clone's column identity."""
        return any((entry.row_index if entry.matrix == 'U' else entry.column_index) > 0 for entry in self.entries)


def _format_error(source, line_number, message):
    return FixedMatrixFormatError(f'{source}:{line_number}: {message}')


def _clone_index(label, clones, source, line_number):
    if label == 'normal':
        return 0

    match = _CLONE_PATTERN.fullmatch(label)
    if match is None:
        raise _format_error(
            source,
            line_number,
            f'invalid clone label {label!r}; expected "normal" or "cloneN" with N >= 1',
        )

    clone_index = int(match.group(1))
    if clone_index >= clones:
        raise _format_error(
            source,
            line_number,
            f'clone label {label!r} is outside the declared clone count {clones}',
        )
    return clone_index


def _parse_copy_number(value, source, line_number):
    try:
        parsed = int(value)
    except ValueError:
        raise _format_error(source, line_number, f'copy-number value {value!r} must be an integer')

    if str(parsed) != value and f'+{parsed}' != value:
        raise _format_error(source, line_number, f'copy-number value {value!r} must be an integer')
    if parsed < 0:
        raise _format_error(source, line_number, 'copy-number values must be nonnegative')
    return parsed


def _parse_proportion(value, source, line_number):
    try:
        parsed = float(value)
    except ValueError:
        raise _format_error(source, line_number, f'proportion value {value!r} must be numeric')

    if not math.isfinite(parsed):
        raise _format_error(source, line_number, 'proportion values must be finite')
    if not 0 <= parsed <= 1:
        raise _format_error(source, line_number, 'proportion values must be in [0, 1]')
    return parsed


def parse_fixed_matrix(path):
    source = os.fspath(path)
    try:
        with open(source, encoding='utf-8-sig') as fixed_file:
            return _parse_fixed_matrix(fixed_file, source)
    except UnicodeDecodeError as exc:
        raise FixedMatrixFormatError(f'{source}: fixed-matrix files must use UTF-8 encoding') from exc


def _parse_fixed_matrix(fixed_file, source):
    first_line = fixed_file.readline()
    if first_line.rstrip('\r\n') != FORMAT_MARKER:
        raise _format_error(source, 1, f'expected format marker {FORMAT_MARKER!r}')

    clones = None
    header_seen = False
    entries = []
    coordinates = set()

    for line_number, raw_line in enumerate(fixed_file, start=2):
        line = raw_line.rstrip('\r\n')
        if not line.strip():
            continue

        stripped_line = line.strip()
        if stripped_line.startswith('#'):
            if stripped_line.startswith('# HATCHET_FIXED_ENTRIES'):
                raise _format_error(source, line_number, 'the format marker may only appear on line 1')

            clones_match = _CLONES_PATTERN.fullmatch(stripped_line)
            if clones_match is not None:
                if header_seen:
                    raise _format_error(source, line_number, 'clone-count metadata must appear before the header')
                if clones is not None:
                    raise _format_error(source, line_number, 'duplicate clone-count metadata')
                clones = int(clones_match.group(1))
                if clones < 2:
                    raise _format_error(source, line_number, 'the clone count must be at least 2, including normal')
            continue

        if not header_seen:
            if clones is None:
                raise _format_error(source, line_number, 'clone-count metadata must appear before the header')
            if line != HEADER:
                raise _format_error(source, line_number, f'expected tab-separated header {HEADER!r}')
            header_seen = True
            continue

        if line == HEADER:
            raise _format_error(source, line_number, 'duplicate header')

        fields = line.split('\t')
        if len(fields) != 4:
            raise _format_error(source, line_number, f'expected 4 tab-separated fields, found {len(fields)}')

        matrix, row, column, raw_value = fields
        if matrix not in ('A', 'B', 'U'):
            raise _format_error(source, line_number, f'unknown matrix {matrix!r}; expected A, B, or U')
        if not row:
            raise _format_error(source, line_number, 'row identifier must not be empty')
        if not column:
            raise _format_error(source, line_number, 'column identifier must not be empty')

        clone_label = column if matrix in ('A', 'B') else row
        clone_index = _clone_index(clone_label, clones, source, line_number)
        value_text = raw_value.strip()
        if matrix in ('A', 'B'):
            value = _parse_copy_number(value_text, source, line_number)
            if clone_index == 0 and value != 1:
                raise _format_error(source, line_number, f'{matrix} entries for the normal clone must equal 1')
        else:
            value = _parse_proportion(value_text, source, line_number)

        coordinate = matrix, row, column
        if coordinate in coordinates:
            raise _format_error(source, line_number, f'duplicate fixed entry for {matrix}[{row}, {column}]')
        coordinates.add(coordinate)
        entries.append(FixedMatrixEntry(matrix, row, column, value, line_number))

    final_line = max(1, line_number if 'line_number' in locals() else 1)
    if clones is None:
        raise _format_error(source, final_line, 'missing clone-count metadata')
    if not header_seen:
        raise _format_error(source, final_line, 'missing matrix entry header')
    if not entries:
        raise _format_error(source, final_line, 'the file must contain at least one fixed entry')

    return FixedMatrix(clones=clones, entries=tuple(entries), source=source)


def _label_lookup(labels, label_type, source):
    lookup = {}
    for index, label in enumerate(labels):
        text = str(label)
        if text in lookup:
            raise FixedMatrixFormatError(
                f'{source}: input {label_type} identifiers are ambiguous when represented as text: {text!r}'
            )
        lookup[text] = index
    return lookup


def resolve_fixed_matrix(fixed_matrix, cluster_ids, sample_ids):
    cluster_ids = tuple(cluster_ids)
    sample_ids = tuple(sample_ids)
    cluster_lookup = _label_lookup(cluster_ids, 'cluster', fixed_matrix.source)
    sample_lookup = _label_lookup(sample_ids, 'sample', fixed_matrix.source)
    resolved_entries = []

    for entry in fixed_matrix.entries:
        if entry.matrix in ('A', 'B'):
            if entry.row not in cluster_lookup:
                raise _format_error(
                    fixed_matrix.source,
                    entry.line_number,
                    f'unknown SEG cluster identifier {entry.row!r}',
                )
            row_index = cluster_lookup[entry.row]
            column_index = _clone_index(
                entry.column,
                fixed_matrix.clones,
                fixed_matrix.source,
                entry.line_number,
            )
        else:
            row_index = _clone_index(
                entry.row,
                fixed_matrix.clones,
                fixed_matrix.source,
                entry.line_number,
            )
            if entry.column not in sample_lookup:
                raise _format_error(
                    fixed_matrix.source,
                    entry.line_number,
                    f'unknown SEG sample identifier {entry.column!r}',
                )
            column_index = sample_lookup[entry.column]

        resolved_entries.append(
            ResolvedFixedMatrixEntry(
                matrix=entry.matrix,
                row_index=row_index,
                column_index=column_index,
                value=entry.value,
                source_entry=entry,
            )
        )

    return ResolvedFixedMatrix(
        clones=fixed_matrix.clones,
        entries=tuple(resolved_entries),
        cluster_ids=cluster_ids,
        sample_ids=sample_ids,
        source=fixed_matrix.source,
    )


def _resolved_values(fixed_matrix):
    return {(entry.matrix, entry.row_index, entry.column_index): entry for entry in fixed_matrix.entries}


def _copy_numbers_by_cluster(copy_numbers):
    return {str(cluster_id): tuple(state) for cluster_id, state in (copy_numbers or {}).items()}


def _weight_by_cluster(weights):
    if weights is None:
        return None
    if hasattr(weights, 'items'):
        return {str(cluster_id): float(weight) for cluster_id, weight in weights.items()}
    return None


def _entry_error(fixed_matrix, entry, message):
    return _format_error(fixed_matrix.source, entry.source_entry.line_number, message)


def validate_fixed_matrix(
    fixed_matrix,
    cn_max,
    mu,
    copy_numbers=None,
    purities=None,
    sample_clonal=None,
    ampdel=True,
    base=None,
    cnstates=None,
    weights=None,
    tol=0.001,
    max_homdel_size=0.005,
    ampdel_exempt=None,
):
    values = _resolved_values(fixed_matrix)
    copy_numbers = _copy_numbers_by_cluster(copy_numbers)
    ampdel_exempt = {str(cluster_id) for cluster_id in (ampdel_exempt or ())}
    if base is None:
        base = min(2, len(copy_numbers))

    for entry in fixed_matrix.entries:
        if entry.matrix not in ('A', 'B'):
            continue

        cluster_label = str(fixed_matrix.cluster_ids[entry.row_index])
        upper_bound = cn_max if cn_max is not None and cn_max > 0 else None
        if cluster_label in copy_numbers:
            clonal_state = copy_numbers[cluster_label]
            clonal_total = sum(clonal_state)
            upper_bound = clonal_total if upper_bound is None else max(upper_bound, clonal_total)
            if entry.column_index > 0:
                allele_index = 0 if entry.matrix == 'A' else 1
                if entry.value != clonal_state[allele_index]:
                    raise _entry_error(
                        fixed_matrix,
                        entry,
                        f'{entry.matrix}[{cluster_label}, {entry.source_entry.column}] conflicts with '
                        f'clonal copy number {clonal_state[allele_index]}',
                    )

        if upper_bound is not None and entry.value > upper_bound:
            raise _entry_error(
                fixed_matrix,
                entry,
                f'{entry.matrix}[{cluster_label}, {entry.source_entry.column}]={entry.value} exceeds '
                f'copy-number bound {upper_bound}',
            )

    for cluster_index, cluster_id in enumerate(fixed_matrix.cluster_ids):
        cluster_label = str(cluster_id)
        upper_bound = cn_max if cn_max is not None and cn_max > 0 else None
        if cluster_label in copy_numbers:
            clonal_total = sum(copy_numbers[cluster_label])
            upper_bound = clonal_total if upper_bound is None else max(upper_bound, clonal_total)

        for clone_index in range(fixed_matrix.clones):
            entry_a = values.get(('A', cluster_index, clone_index))
            entry_b = values.get(('B', cluster_index, clone_index))
            if entry_a is not None and entry_b is not None and upper_bound is not None:
                total = entry_a.value + entry_b.value
                if total > upper_bound:
                    raise _entry_error(
                        fixed_matrix,
                        entry_b,
                        f'fixed A+B total {total} for cluster {cluster_label!r}, clone {clone_index} '
                        f'exceeds copy-number bound {upper_bound}',
                    )

    if mu is not None and mu > 0:
        for entry in fixed_matrix.entries_for('U'):
            if entry.row_index > 0 and tol < entry.value < mu - tol:
                raise _entry_error(
                    fixed_matrix,
                    entry,
                    f'fixed tumor proportion {entry.value} is between 0 and minprop {mu}',
                )

    for sample_index, sample_id in enumerate(fixed_matrix.sample_ids):
        sample_entries = [
            values[('U', clone_index, sample_index)]
            for clone_index in range(fixed_matrix.clones)
            if ('U', clone_index, sample_index) in values
        ]
        if not sample_entries:
            continue

        fixed_sum = sum(entry.value for entry in sample_entries)
        error_entry = max(sample_entries, key=lambda entry: entry.source_entry.line_number)
        if fixed_sum > 1 + tol:
            raise _entry_error(
                fixed_matrix,
                error_entry,
                f'fixed U entries for sample {sample_id!r} sum to {fixed_sum}, which exceeds 1',
            )

        free_clones = [
            clone_index for clone_index in range(fixed_matrix.clones) if ('U', clone_index, sample_index) not in values
        ]
        residual = 1 - fixed_sum
        if not free_clones and abs(residual) > tol:
            raise _entry_error(
                fixed_matrix,
                error_entry,
                f'fully fixed U entries for sample {sample_id!r} must sum to 1, found {fixed_sum}',
            )
        if free_clones and 0 not in free_clones and mu is not None and mu > 0 and tol < residual < mu - tol:
            raise _entry_error(
                fixed_matrix,
                error_entry,
                f'remaining U proportion {residual} for sample {sample_id!r} is below minprop {mu}',
            )

    if purities is not None:
        if len(purities) != len(fixed_matrix.sample_ids):
            raise FixedMatrixFormatError(
                f'{fixed_matrix.source}: purity count {len(purities)} does not match '
                f'sample count {len(fixed_matrix.sample_ids)}'
            )
        for sample_index, purity in enumerate(purities):
            entry = values.get(('U', 0, sample_index))
            expected = 1 - purity
            if entry is not None and abs(entry.value - expected) > tol:
                raise _entry_error(
                    fixed_matrix,
                    entry,
                    f'fixed normal proportion {entry.value} conflicts with purity {purity}; expected {expected}',
                )

    if sample_clonal is not None:
        if purities is None or len(purities) != len(fixed_matrix.sample_ids):
            raise FixedMatrixFormatError(
                f'{fixed_matrix.source}: sample-clonal validation requires one purity per sample'
            )

        states_by_sample = {
            sample_index: {
                entry.cluster_index: (entry.cn_a, entry.cn_b)
                for entry in sample_clonal.entries_for_sample(sample_index)
            }
            for sample_index in range(len(fixed_matrix.sample_ids))
        }

        def clone_compatible(clone_index, sample_index):
            for cluster_index, state in states_by_sample[sample_index].items():
                entry_a = values.get(('A', cluster_index, clone_index))
                entry_b = values.get(('B', cluster_index, clone_index))
                if entry_a is not None and entry_a.value != state[0]:
                    return False
                if entry_b is not None and entry_b.value != state[1]:
                    return False
            return True

        fixed_positive_samples = {clone_index: [] for clone_index in range(1, fixed_matrix.clones)}
        for sample_index, sample_id in enumerate(fixed_matrix.sample_ids):
            tumor_entries = [
                values[('U', clone_index, sample_index)]
                for clone_index in range(1, fixed_matrix.clones)
                if ('U', clone_index, sample_index) in values
            ]
            fixed_tumor_sum = sum(entry.value for entry in tumor_entries)
            purity = purities[sample_index]
            if fixed_tumor_sum > purity + tol:
                error_entry = max(tumor_entries, key=lambda entry: entry.source_entry.line_number)
                raise _entry_error(
                    fixed_matrix,
                    error_entry,
                    f'fixed tumor proportions for sample {sample_id!r} sum to {fixed_tumor_sum}, '
                    f'which exceeds sample-clonal purity {purity}',
                )

            for entry in tumor_entries:
                if entry.value <= tol:
                    continue
                if not clone_compatible(entry.row_index, sample_index):
                    raise _entry_error(
                        fixed_matrix,
                        entry,
                        f'fixed-positive {entry.source_entry.row} is incompatible with sample-clonal states '
                        f'for sample {sample_id!r}',
                    )
                fixed_positive_samples[entry.row_index].append(sample_index)

            fixed_clones = {entry.row_index for entry in tumor_entries}
            free_compatible = [
                clone_index
                for clone_index in range(1, fixed_matrix.clones)
                if clone_index not in fixed_clones and clone_compatible(clone_index, sample_index)
            ]
            residual = purity - fixed_tumor_sum
            if residual > tol and not free_compatible:
                error_entry = max(tumor_entries, key=lambda entry: entry.source_entry.line_number, default=None)
                if error_entry is None:
                    raise FixedMatrixFormatError(
                        f'{fixed_matrix.source}: no tumor clone is eligible for remaining sample-clonal purity '
                        f'{residual} in sample {sample_id!r}'
                    )
                raise _entry_error(
                    fixed_matrix,
                    error_entry,
                    f'no tumor clone is eligible for remaining sample-clonal purity {residual} '
                    f'in sample {sample_id!r}',
                )
            if residual > tol and mu is not None and mu > 0 and residual < mu - tol:
                error_entry = max(tumor_entries, key=lambda entry: entry.source_entry.line_number, default=None)
                message = f'remaining sample-clonal purity {residual} for sample {sample_id!r} is below minprop {mu}'
                if error_entry is None:
                    raise FixedMatrixFormatError(f'{fixed_matrix.source}: {message}')
                raise _entry_error(fixed_matrix, error_entry, message)

        for clone_index, sample_indices in fixed_positive_samples.items():
            for left_index, left_sample in enumerate(sample_indices):
                for right_sample in sample_indices[left_index + 1 :]:
                    shared_clusters = set(states_by_sample[left_sample]).intersection(states_by_sample[right_sample])
                    if any(
                        states_by_sample[left_sample][cluster] != states_by_sample[right_sample][cluster]
                        for cluster in shared_clusters
                    ):
                        entry = values[('U', clone_index, right_sample)]
                        raise _entry_error(
                            fixed_matrix,
                            entry,
                            f'fixed-positive {entry.source_entry.row} is assigned to samples with incompatible '
                            f'sample-clonal states',
                        )

    if ampdel:
        for cluster_index, cluster_id in enumerate(fixed_matrix.cluster_ids):
            if str(cluster_id) in copy_numbers or str(cluster_id) in ampdel_exempt:
                continue
            for matrix in ('A', 'B'):
                entries = [
                    values[(matrix, cluster_index, clone_index)]
                    for clone_index in range(1, fixed_matrix.clones)
                    if (matrix, cluster_index, clone_index) in values
                ]
                if any(entry.value < base for entry in entries) and any(entry.value > base for entry in entries):
                    error_entry = max(entries, key=lambda entry: entry.source_entry.line_number)
                    raise _entry_error(
                        fixed_matrix,
                        error_entry,
                        f'fixed {matrix} entries for cluster {str(cluster_id)!r} force both amplification '
                        f'and deletion relative to base {base}',
                    )

    weight_map = _weight_by_cluster(weights)
    if weight_map:
        total_weight = sum(weight_map.values())
        if total_weight > 0:
            for cluster_index, cluster_id in enumerate(fixed_matrix.cluster_ids):
                cluster_label = str(cluster_id)
                cluster_weight = weight_map.get(cluster_label, 0)
                if cluster_weight <= 0 or cluster_weight / total_weight < max_homdel_size:
                    continue
                for clone_index in range(1, fixed_matrix.clones):
                    entry_a = values.get(('A', cluster_index, clone_index))
                    entry_b = values.get(('B', cluster_index, clone_index))
                    if entry_a is not None and entry_b is not None and entry_a.value + entry_b.value == 0:
                        raise _entry_error(
                            fixed_matrix,
                            entry_b,
                            f'large cluster {cluster_label!r} cannot have fixed tumor copy number (0, 0)',
                        )

    if cnstates is not None and cnstates > 0:
        for cluster_index, cluster_id in enumerate(fixed_matrix.cluster_ids):
            complete_states = []
            for clone_index in range(1, fixed_matrix.clones):
                entry_a = values.get(('A', cluster_index, clone_index))
                entry_b = values.get(('B', cluster_index, clone_index))
                if entry_a is not None and entry_b is not None:
                    complete_states.append((entry_a.value, entry_b.value, entry_b))
            distinct_states = {(state_a, state_b) for state_a, state_b, _ in complete_states}
            if len(distinct_states) > cnstates:
                error_entry = max(
                    (entry for _, _, entry in complete_states),
                    key=lambda entry: entry.source_entry.line_number,
                )
                raise _entry_error(
                    fixed_matrix,
                    error_entry,
                    f'cluster {str(cluster_id)!r} fixes {len(distinct_states)} distinct tumor states, '
                    f'exceeding cnstates {cnstates}',
                )

    return fixed_matrix


def verify_fixed_matrix_solution(fixed_matrix, cA, cB, u, tol=0.001):
    if fixed_matrix is None:
        return

    matrices = {'A': cA, 'B': cB, 'U': u}
    for entry in fixed_matrix.entries:
        actual = matrices[entry.matrix][entry.row_index][entry.column_index]
        if entry.matrix in ('A', 'B'):
            matches = int(round(actual)) == entry.value
        else:
            matches = abs(actual - entry.value) <= tol
        if not matches:
            source_entry = entry.source_entry
            raise RuntimeError(
                f'{fixed_matrix.source}:{source_entry.line_number}: solver returned {actual} for '
                f'{entry.matrix}[{source_entry.row}, {source_entry.column}], expected {entry.value}'
            )
