import math
import os
from dataclasses import dataclass
from typing import Tuple


FORMAT_MARKER = '# HATCHET_SAMPLE_CLONAL 1'
HEADER = 'sample\tcluster\tcn_a\tcn_b'
PURITY_TOLERANCE = 0.05


class SampleClonalFormatError(ValueError):
    pass


@dataclass(frozen=True)
class SampleClonalEntry:
    sample: str
    cluster: str
    cn_a: int
    cn_b: int
    line_number: int

    @property
    def total_copy_number(self):
        return self.cn_a + self.cn_b


@dataclass(frozen=True)
class SampleClonal:
    entries: Tuple[SampleClonalEntry, ...]
    source: str

    def entries_for_sample(self, sample):
        return tuple(entry for entry in self.entries if entry.sample == sample)


@dataclass(frozen=True)
class ResolvedSampleClonalEntry:
    sample_index: int
    cluster_index: int
    source_entry: SampleClonalEntry

    @property
    def cn_a(self):
        return self.source_entry.cn_a

    @property
    def cn_b(self):
        return self.source_entry.cn_b

    @property
    def total_copy_number(self):
        return self.source_entry.total_copy_number


@dataclass(frozen=True)
class ResolvedSampleClonal:
    entries: Tuple[ResolvedSampleClonalEntry, ...]
    cluster_ids: tuple
    sample_ids: tuple
    source: str

    def entries_for_sample(self, sample_index):
        return tuple(entry for entry in self.entries if entry.sample_index == sample_index)


@dataclass(frozen=True)
class SampleClonalCalibration:
    sample_ids: tuple
    gamma: tuple
    purity: tuple

    def gamma_by_sample(self):
        return dict(zip(self.sample_ids, self.gamma))

    def purity_by_sample(self):
        return dict(zip(self.sample_ids, self.purity))


def _format_error(source, line_number, message):
    return SampleClonalFormatError(f'{source}:{line_number}: {message}')


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


def _make_entry(fields, source, line_number):
    sample, cluster, raw_cn_a, raw_cn_b = fields
    if not sample:
        raise _format_error(source, line_number, 'sample identifier must not be empty')
    if not cluster:
        raise _format_error(source, line_number, 'cluster identifier must not be empty')

    return SampleClonalEntry(
        sample=sample,
        cluster=cluster,
        cn_a=_parse_copy_number(raw_cn_a.strip(), source, line_number),
        cn_b=_parse_copy_number(raw_cn_b.strip(), source, line_number),
        line_number=line_number,
    )


def _validate_entries(entries, source):
    if not entries:
        raise _format_error(source, 1, 'the specification must contain at least one entry')

    coordinates = set()
    by_sample = {}
    for entry in entries:
        coordinate = entry.sample, entry.cluster
        if coordinate in coordinates:
            raise _format_error(
                source,
                entry.line_number,
                f'duplicate sample-clonal entry for sample {entry.sample!r}, cluster {entry.cluster!r}',
            )
        coordinates.add(coordinate)
        by_sample.setdefault(entry.sample, []).append(entry)

    for sample, sample_entries in by_sample.items():
        if len(sample_entries) < 2:
            raise _format_error(
                source,
                sample_entries[-1].line_number,
                f'sample {sample!r} must have at least two entries, found {len(sample_entries)}',
            )
        # only the first two entries calibrate the sample; any further entries are only enforced
        if len({entry.total_copy_number for entry in sample_entries[:2]}) != 2:
            raise _format_error(
                source,
                sample_entries[1].line_number,
                f'sample {sample!r} must have two distinct total copy-number states',
            )


def parse_sample_clonal_file(path):
    source = os.fspath(path)
    try:
        with open(source, encoding='utf-8-sig') as sample_clonal_file:
            first_line = sample_clonal_file.readline()
            if first_line.rstrip('\r\n') != FORMAT_MARKER:
                raise _format_error(source, 1, f'expected format marker {FORMAT_MARKER!r}')

            header_seen = False
            entries = []
            for line_number, raw_line in enumerate(sample_clonal_file, start=2):
                line = raw_line.rstrip('\r\n')
                if not line.strip():
                    continue
                if line.startswith('#'):
                    if line.startswith('# HATCHET_SAMPLE_CLONAL'):
                        raise _format_error(source, line_number, 'the format marker may only appear on line 1')
                    continue
                if not header_seen:
                    if line != HEADER:
                        raise _format_error(source, line_number, f'expected tab-separated header {HEADER!r}')
                    header_seen = True
                    continue
                if line == HEADER:
                    raise _format_error(source, line_number, 'duplicate header')

                fields = line.split('\t')
                if len(fields) != 4:
                    raise _format_error(source, line_number, f'expected 4 tab-separated fields, found {len(fields)}')
                entries.append(_make_entry(fields, source, line_number))
    except UnicodeDecodeError as exc:
        raise SampleClonalFormatError(f'{source}: sample-clonal files must use UTF-8 encoding') from exc

    if not header_seen:
        raise _format_error(source, 1, 'missing sample-clonal header')
    _validate_entries(entries, source)
    return SampleClonal(entries=tuple(entries), source=source)


def _validate_inline_identifier(identifier, kind, source, line_number):
    if not identifier:
        raise _format_error(source, line_number, f'{kind} identifier must not be empty')
    if any(delimiter in identifier for delimiter in '=;,:'):
        raise _format_error(
            source,
            line_number,
            f'{kind} identifier {identifier!r} contains a reserved inline delimiter; use a TSV file',
        )


def parse_sample_clonal_inline(value, source='--sample-clonal'):
    entries = []
    for sample_number, raw_sample in enumerate(value.split(';'), start=1):
        sample_part = raw_sample.strip()
        if not sample_part:
            raise _format_error(source, sample_number, 'empty sample specification')
        if sample_part.count('=') != 1:
            raise _format_error(source, sample_number, 'expected SAMPLE=CLUSTER:A:B,CLUSTER:A:B')

        sample, raw_states = (part.strip() for part in sample_part.split('=', 1))
        _validate_inline_identifier(sample, 'sample', source, sample_number)
        states = raw_states.split(',')
        if len(states) < 2:
            raise _format_error(source, sample_number, f'sample {sample!r} must have at least two entries')

        for raw_state in states:
            state_fields = tuple(field.strip() for field in raw_state.strip().split(':'))
            if len(state_fields) != 3:
                raise _format_error(source, sample_number, 'expected each state as CLUSTER:A:B')
            _validate_inline_identifier(state_fields[0], 'cluster', source, sample_number)
            entries.append(_make_entry((sample,) + state_fields, source, sample_number))

    _validate_entries(entries, source)
    return SampleClonal(entries=tuple(entries), source=source)


def _label_lookup(labels, label_type, source):
    lookup = {}
    for index, label in enumerate(labels):
        text = str(label)
        if text in lookup:
            raise SampleClonalFormatError(
                f'{source}: input {label_type} identifiers are ambiguous when represented as text: {text!r}'
            )
        lookup[text] = index
    return lookup


def resolve_sample_clonal(sample_clonal, cluster_ids, sample_ids):
    cluster_ids = tuple(cluster_ids)
    sample_ids = tuple(sample_ids)
    cluster_lookup = _label_lookup(cluster_ids, 'cluster', sample_clonal.source)
    sample_lookup = _label_lookup(sample_ids, 'sample', sample_clonal.source)
    resolved_entries = []
    supplied_samples = set()

    for entry in sample_clonal.entries:
        if entry.sample not in sample_lookup:
            raise _format_error(
                sample_clonal.source,
                entry.line_number,
                f'unknown SEG sample identifier {entry.sample!r}',
            )
        if entry.cluster not in cluster_lookup:
            raise _format_error(
                sample_clonal.source,
                entry.line_number,
                f'unknown SEG cluster identifier {entry.cluster!r}',
            )
        supplied_samples.add(entry.sample)
        resolved_entries.append(
            ResolvedSampleClonalEntry(
                sample_index=sample_lookup[entry.sample],
                cluster_index=cluster_lookup[entry.cluster],
                source_entry=entry,
            )
        )

    missing_samples = [str(sample) for sample in sample_ids if str(sample) not in supplied_samples]
    if missing_samples:
        raise SampleClonalFormatError(
            f'{sample_clonal.source}: missing sample-clonal entries for SEG samples {missing_samples}'
        )

    return ResolvedSampleClonal(
        entries=tuple(resolved_entries),
        cluster_ids=cluster_ids,
        sample_ids=sample_ids,
        source=sample_clonal.source,
    )


def reindex_sample_clonal(sample_clonal, cluster_ids, sample_ids):
    if isinstance(sample_clonal, ResolvedSampleClonal):
        sample_clonal = SampleClonal(
            entries=tuple(entry.source_entry for entry in sample_clonal.entries),
            source=sample_clonal.source,
        )
    return resolve_sample_clonal(sample_clonal, cluster_ids, sample_ids)


def write_sample_clonal_file(sample_clonal, path):
    entries = tuple(getattr(entry, 'source_entry', entry) for entry in sample_clonal.entries)
    with open(path, 'w', encoding='utf-8', newline='') as output:
        output.write(FORMAT_MARKER + '\n')
        output.write(HEADER + '\n')
        for entry in entries:
            output.write(f'{entry.sample}\t{entry.cluster}\t{entry.cn_a}\t{entry.cn_b}\n')


def derive_sample_clonal_calibration(resolved, rdr, purity_tolerance=PURITY_TOLERANCE):
    gamma = []
    purity = []

    for sample_index, sample_id in enumerate(resolved.sample_ids):
        # entries after the first two are enforced during solving but do not calibrate the sample
        first, second = resolved.entries_for_sample(sample_index)[:2]
        rdr_first = float(rdr.iloc[first.cluster_index, sample_index])
        rdr_second = float(rdr.iloc[second.cluster_index, sample_index])
        total_first = first.total_copy_number
        total_second = second.total_copy_number

        if not math.isfinite(rdr_first) or not math.isfinite(rdr_second):
            raise SampleClonalFormatError(f'{resolved.source}: RDR values for sample {str(sample_id)!r} must be finite')

        denominator = rdr_first * (total_second - 2) + rdr_second * (2 - total_first)
        if math.isclose(denominator, 0.0, rel_tol=0.0, abs_tol=1e-12):
            raise SampleClonalFormatError(
                f'{resolved.source}: sample {str(sample_id)!r} has a singular sample-clonal calibration system'
            )

        sample_gamma = 2.0 * (total_second - total_first) / denominator
        non_diploid_entry, non_diploid_rdr = (first, rdr_first) if total_first != 2 else (second, rdr_second)
        sample_purity = (sample_gamma * non_diploid_rdr - 2.0) / (non_diploid_entry.total_copy_number - 2)

        if 1 <= sample_purity <= 1 + purity_tolerance:
            sample_purity = 1.0
        elif -purity_tolerance <= sample_purity <= 0:
            sample_purity = 0.0

        if not math.isfinite(sample_gamma) or sample_gamma <= 0:
            raise SampleClonalFormatError(
                f'{resolved.source}: sample {str(sample_id)!r} produced invalid scaling factor {sample_gamma}'
            )
        if not math.isfinite(sample_purity) or not 0 <= sample_purity <= 1:
            raise SampleClonalFormatError(
                f'{resolved.source}: sample {str(sample_id)!r} produced invalid purity {sample_purity}'
            )

        gamma.append(sample_gamma)
        purity.append(sample_purity)

    return SampleClonalCalibration(
        sample_ids=resolved.sample_ids,
        gamma=tuple(gamma),
        purity=tuple(purity),
    )


def validate_sample_clonal_constraints(
    resolved,
    purities,
    cn_max,
    base,
    ampdel=True,
    cnstates=None,
    weights=None,
    tol=0.001,
    max_homdel_size=0.005,
    ampdel_exempt=None,
):
    positive_entries = [entry for entry in resolved.entries if purities[entry.sample_index] > tol]
    ampdel_exempt = {str(cluster_id) for cluster_id in (ampdel_exempt or ())}
    weight_by_cluster = None
    if weights is not None:
        if hasattr(weights, 'items'):
            weight_by_cluster = {str(cluster): float(weight) for cluster, weight in weights.items()}
        else:
            weight_by_cluster = {
                str(cluster): float(weights[index]) for index, cluster in enumerate(resolved.cluster_ids)
            }
    total_weight = sum(weight_by_cluster.values()) if weight_by_cluster else 0.0

    for entry in positive_entries:
        if entry.total_copy_number > cn_max:
            raise _format_error(
                resolved.source,
                entry.source_entry.line_number,
                f'sample-clonal total copy number {entry.total_copy_number} exceeds copy-number bound {cn_max}',
            )
        cluster_id = str(resolved.cluster_ids[entry.cluster_index])
        cluster_weight = weight_by_cluster.get(cluster_id, 0.0) if weight_by_cluster else 0.0
        if (
            entry.total_copy_number == 0
            and total_weight > 0
            and cluster_weight > 0
            and cluster_weight / total_weight >= max_homdel_size
        ):
            raise _format_error(
                resolved.source,
                entry.source_entry.line_number,
                f'large cluster {cluster_id!r} cannot have required tumor copy number (0, 0)',
            )

    entries_by_cluster = {}
    for entry in positive_entries:
        entries_by_cluster.setdefault(entry.cluster_index, []).append(entry)

    if ampdel:
        for cluster_index, entries in entries_by_cluster.items():
            if str(resolved.cluster_ids[cluster_index]) in ampdel_exempt:
                continue
            for allele, attribute in (('A', 'cn_a'), ('B', 'cn_b')):
                values = [getattr(entry, attribute) for entry in entries]
                if any(value < base for value in values) and any(value > base for value in values):
                    error_entry = max(entries, key=lambda entry: entry.source_entry.line_number)
                    raise _format_error(
                        resolved.source,
                        error_entry.source_entry.line_number,
                        f'sample-clonal states force both amplification and deletion of allele {allele} '
                        f'for cluster {str(resolved.cluster_ids[cluster_index])!r} relative to base {base}',
                    )

    if cnstates is not None and cnstates > 0:
        for cluster_index, entries in entries_by_cluster.items():
            states = {(entry.cn_a, entry.cn_b) for entry in entries}
            if len(states) > cnstates:
                error_entry = max(entries, key=lambda entry: entry.source_entry.line_number)
                raise _format_error(
                    resolved.source,
                    error_entry.source_entry.line_number,
                    f'sample-clonal states for cluster {str(resolved.cluster_ids[cluster_index])!r} '
                    f'require at least {len(states)} copy-number states, exceeding cnstates {cnstates}',
                )


def verify_sample_clonal_solution(resolved, purities, c_a, c_b, u, tol=0.001):
    for sample_index, sample_id in enumerate(resolved.sample_ids):
        expected_normal = 1 - purities[sample_index]
        if abs(u[0][sample_index] - expected_normal) > tol:
            raise RuntimeError(
                f'Sample-clonal solution has normal proportion {u[0][sample_index]} for sample '
                f'{str(sample_id)!r}; expected {expected_normal}'
            )
        for clone_index in range(1, len(u)):
            if u[clone_index][sample_index] <= tol:
                continue
            for entry in resolved.entries_for_sample(sample_index):
                if c_a[entry.cluster_index][clone_index] != entry.cn_a:
                    raise RuntimeError(
                        f'Sample-clonal solution violates A[{str(resolved.cluster_ids[entry.cluster_index])!r}, '
                        f'clone{clone_index}] for sample {str(sample_id)!r}'
                    )
                if c_b[entry.cluster_index][clone_index] != entry.cn_b:
                    raise RuntimeError(
                        f'Sample-clonal solution violates B[{str(resolved.cluster_ids[entry.cluster_index])!r}, '
                        f'clone{clone_index}] for sample {str(sample_id)!r}'
                    )
