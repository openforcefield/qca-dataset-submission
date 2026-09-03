#!/usr/bin/env python3
"""Enumerate feasible charge/multiplicity variants for a metal complex.

Given a metal center and one known (total charge, multiplicity) pair -- but *not*
the oxidation state -- enumerate the other (charge, multiplicity) combinations
that are physically reasonable for the same ligand set.

Assumptions
-----------
**The ligands are closed-shell, and every electron gained or lost is the metal's.**
Two consequences follow. The metal alone sets the total multiplicity, so a
multiplicity that no oxidation state of that metal can produce is rejected as a
bad record. And the ligand set never changes, so a change in total charge shifts
the oxidation state by the same amount.

Redox-non-innocent ligands (dithiolenes, quinones, NO, reduced bipyridines) break
this by taking the electrons themselves, leaving the metal oxidation state
unchanged. Such a complex still gets the right charge and multiplicity but the
wrong oxidation state, so use the variants and distrust `Variant.assignments`.

Ligand charge itself is unconstrained and never checked. Whether an input is
accepted depends only on `metal` and `multiplicity`, never on `charge`, so a
charge carried by the ligands cannot cause a spurious rejection.

The oxidation state cannot be recovered from total charge and multiplicity alone,
so every oxidation state whose spin manifold contains the input multiplicity is
carried as a hypothesis and their results are unioned. Each hypothesis implies a
ligand-set charge ``q_ligands = q_input - ox``, and reaching a target total charge
means the metal absorbs the difference::

    ox_target = q_target - q_ligands = ox + (q_target - q_input)

Multiplicities allowed at ``ox_target`` come from the d-electron table, which
encodes how many d orbitals are free to host unpaired spins: d6 admits 1/3/5, d8
only 1/3, d10 only 1. A target charge equal to the input charge therefore still
yields spin-state variants, whenever the d manifold has room. Results are then
filtered by a draft per-element accessible-oxidation-state table scoped to this
dataset's chemistry. Elements absent from that table are left unrestricted.

Because charge changes land on the metal, multiplicity parity is tied to charge
parity: an odd change in
total charge always flips multiplicity parity. Parity holds regardless of where
the charge sits, so `parity_consistent` checks a record against its own electron
count without reference to either assumption.

Oxidation states requiring more than `MAX_D_ELECTRONS` d electrons are dropped as
unphysical, so an element's states may be non-contiguous and need not start at 0:
Cu starts at 1+ and Zn at 2+.

Examples
--------
Variants of a neutral low-spin Fe complex, at total charge -1, 0 and +1:

>>> variants = enumerate_variants("Fe", 0, 1)
>>> [(v.charge, v.multiplicity) for v in variants if not v.is_input]
[(-1, 2), (-1, 4), (-1, 6), (0, 3), (0, 5), (1, 2), (1, 4), (1, 6)]

The oxidation states supporting a given variant, here the neutral quintet:

>>> [v.assignments for v in variants if (v.charge, v.multiplicity) == (0, 5)]
[[Assignment(oxidation_state=2, configuration='d6'), Assignment(oxidation_state=4, configuration='d4')]]

Zn is singlet-only, so a neutral Zn complex has no variants but itself:

>>> variants = enumerate_variants("Zn", 0, 1)
>>> [(v.charge, v.multiplicity, v.is_input) for v in variants]
[(0, 1, True)]

The draft oxidation-state table treats Pd triplets at total charge +1 as Pd(II),
not Pd(IV) or Pd(VI), so reducing by two charges does not keep the triplet:

>>> variants = enumerate_variants("Pd", 1, 3)
>>> [(v.charge, v.multiplicity) for v in variants]
[(-1, 1), (1, 1), (1, 3)]

Lifting the oxidation-state filter admits states the d manifold allows but the
dataset's chemistry likely does not:

>>> variants = enumerate_variants("Ni", 0, 1, restrict_to_accessible=False)
>>> [(v.charge, v.multiplicity) for v in variants]
[(-1, 2), (-1, 4), (-1, 6), (0, 1), (0, 3), (0, 5), (1, 2), (1, 4), (1, 6)]
"""

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

# Total charges permitted for an output complex, per dataset scope.
DEFAULT_CHARGES = (-1, 0, 1)

# A d shell holds at most 10 electrons. The report tabulates d11/d12 as
# bookkeeping shorthand for nd10 (n+1)s^k, which are not real d configurations;
# `_flatten_spin` drops them so they are never offered as oxidation-state
# hypotheses. They must stay in the tables below regardless, because list
# position encodes the oxidation state -- deleting a row shifts every state
# after it.
MAX_D_ELECTRONS = 10

# ---------------------------------------------------------------------------
# Reference data, transcribed from .agents/reports/atom_substitutions.md.
#
# Per group: (d-electron count, allowed spin multiplicities) for oxidation
# states 0, 1+, 2+, ... in order. POSITION IS THE OXIDATION STATE: entry i is
# the i+ state, so rows may be extended but never reordered or shortened.
# ---------------------------------------------------------------------------

_SPIN_GROUPS: Dict[Tuple[str, ...], List[Tuple[int, List[int]]]] = {
    ("Sc", "Y", "La"): [(3, [2, 4]), (2, [1, 3]), (1, [2]), (0, [1])],
    ("Ti", "Zr", "Hf"): [(4, [1, 3, 5]), (3, [2, 4]), (2, [1, 3]), (1, [2]), (0, [1])],
    ("V", "Nb", "Ta"): [(5, [2, 4, 6]), (4, [1, 3, 5]), (3, [2, 4]), (2, [1, 3]),
                        (1, [2]), (0, [1])],
    ("Cr", "Mo", "W"): [(6, [1, 3, 5]), (5, [2, 4, 6]), (4, [1, 3, 5]), (3, [2, 4]),
                        (2, [1, 3]), (1, [2]), (0, [1])],
    ("Mn", "Tc", "Re"): [(7, [2, 4]), (6, [1, 3, 5]), (5, [2, 4, 6]), (4, [1, 3, 5]),
                         (3, [2, 4]), (2, [1, 3]), (1, [2])],
    ("Fe", "Ru", "Os"): [(8, [1, 3]), (7, [2, 4]), (6, [1, 3, 5]), (5, [2, 4, 6]),
                         (4, [1, 3, 5]), (3, [2, 4]), (2, [1, 3])],
    ("Co", "Rh", "Ir"): [(9, [2]), (8, [1, 3]), (7, [2, 4]), (6, [1, 3, 5]),
                         (5, [2, 4, 6]), (4, [1, 3, 5]), (3, [2, 4])],
    ("Ni", "Pd", "Pt"): [(10, [1]), (9, [2]), (8, [1, 3]), (7, [2, 4]), (6, [1, 3, 5]),
                         (5, [2, 4, 6]), (4, [1, 3, 5])],
    ("Cu", "Ag", "Au"): [(11, []), (10, [1]), (9, [2]), (8, [1, 3]), (7, [2, 4]),
                         (6, [1, 3, 5]), (5, [2, 4, 6])],
    ("Zn", "Cd", "Hg"): [(12, []), (11, []), (10, [1]), (9, [2]), (8, [1, 3]),
                         (7, [2, 4]), (6, [1, 3, 5])],
}

# Main-group centers have no d manifold; spin follows s/p occupancy alone.
# Not in the report's d-count table, which covers transition metals only.
_MAIN_GROUP_SPIN: Dict[Tuple[str, ...], List[Tuple[int, List[int]]]] = {
    ("Li", "Na", "K"): [(1, [2]), (0, [1])],
    ("Be", "Mg", "Ca"): [(2, [1]), (1, [2]), (0, [1])],
}

# Draft accessible oxidation states per element, for this dataset's chemistry.
# Elements absent here are left unrestricted.
_ACCESSIBLE_OXIDATION_GROUPS: Dict[Tuple[str, ...], Set[int]] = {
    ("Li",): {0, 1},
    ("Mg",): {0, 2},
    ("Sc", "Y"): {3},
    ("Ti", "Zr"): {2, 3, 4},
    ("V", "Nb"): {2, 3, 4, 5},
    ("Cr", "Mo"): {0, 1, 2, 3, 4, 5, 6},
    ("Mn", "Tc"): {2, 3, 4, 5, 6},
    ("Fe", "Ru"): {2, 3, 4},
    ("Co", "Rh"): {0, 1, 2, 3},
    ("Ni", "Pd"): {0, 1, 2, 3},
    ("Cu", "Ag"): {1, 2, 3},
    ("Zn", "Cd"): {2},
}


def _flatten_spin() -> Dict[str, Dict[int, Tuple[int, List[int]]]]:
    """Expand the grouped spin tables to one entry per element.

    List position in the source tables is the oxidation state. Entries whose d
    count exceeds `MAX_D_ELECTRONS` are the report's nd10 (n+1)s^k shorthand
    rather than real d configurations, and are dropped here so no caller can
    propose them; because the result is keyed by oxidation state rather than
    position, dropping them shifts nothing. Main-group counts are s occupancy
    and are never filtered.

    Returns
    -------
    dict
        ``{element: {oxidation_state: (electron_count, [multiplicities])}}``,
        where ``electron_count`` is the d count for transition metals and the s
        count for main-group centers. Oxidation states may be non-contiguous:
        Cu starts at 1+ and Zn at 2+, the lower states having been filtered.
    """
    out: Dict[str, Dict[int, Tuple[int, List[int]]]] = {}
    for source, is_d in ((_SPIN_GROUPS, True), (_MAIN_GROUP_SPIN, False)):
        for elements, states in source.items():
            for element in elements:
                out[element] = {
                    ox: (count, list(mults))
                    for ox, (count, mults) in enumerate(states)
                    if not is_d or count <= MAX_D_ELECTRONS
                }
    return out


def _flatten_accessible_oxidation_states() -> Dict[str, Set[int]]:
    """Expand the grouped accessible-oxidation-state table to one entry per element.

    Returns
    -------
    dict
        ``{element: {accessible oxidation states}}``. Elements absent from the
        table are absent here, and are treated as unrestricted.
    """
    return {element: set(states)
            for elements, states in _ACCESSIBLE_OXIDATION_GROUPS.items()
            for element in elements}


SPIN_STATES = _flatten_spin()
ACCESSIBLE_OXIDATION_STATES = _flatten_accessible_oxidation_states()

# s-block centers have no d manifold; their counts are s occupancy, not d.
MAIN_GROUP = {el for els in _MAIN_GROUP_SPIN for el in els}


def shell(metal: str) -> str:
    """Return the orbital label for an element's tabulated electron count.

    Parameters
    ----------
    metal : str
        Element symbol, capitalized.

    Returns
    -------
    str
        ``"s"`` for main-group centers, ``"d"`` for transition metals. Use it to
        label `Variant.assignments` counts, which are s occupancy for the former
        and d occupancy for the latter.

    Examples
    --------
    >>> shell("Fe"), shell("Li")
    ('d', 's')
    """
    return "s" if metal in MAIN_GROUP else "d"


def parity_consistent(neutral_electrons: int, charge: int, multiplicity: int) -> bool:
    """Test a (charge, multiplicity) pair against a molecule's electron count.

    Unpaired electrons cannot change the parity of a total electron count, so a
    complex holding ``N = neutral_electrons - charge`` electrons must have
    ``multiplicity - 1`` unpaired electrons of the same parity as ``N``. An
    even-electron species is therefore odd in multiplicity, and vice versa.

    Parameters
    ----------
    neutral_electrons : int
        Electron count of the *neutral* species, i.e. the sum of atomic numbers
        over all atoms. Independent of `charge`, which is applied here.
    charge : int
        Total charge of the complex.
    multiplicity : int
        Spin multiplicity (unpaired electrons + 1).

    Returns
    -------
    bool
        True if the combination is possible.

    Examples
    --------
    Neutral benzene, 42 electrons, is even and so must be a singlet:

    >>> parity_consistent(42, 0, 1), parity_consistent(42, 0, 2)
    (True, False)

    Removing one electron makes it odd, and so a doublet:

    >>> parity_consistent(42, 1, 2), parity_consistent(42, 1, 1)
    (True, False)
    """
    return (multiplicity - 1 - (neutral_electrons - charge)) % 2 == 0


def electron_count(atomic_numbers: Iterable[int], charge: int = 0) -> int:
    """Return the electron count of a species from its atomic numbers.

    Parameters
    ----------
    atomic_numbers : iterable of int
        Atomic number of every atom in the molecule.
    charge : int, optional
        Total charge. Default 0, giving the neutral count that
        `parity_consistent` expects.

    Returns
    -------
    int
        ``sum(atomic_numbers) - charge``.

    Examples
    --------
    >>> electron_count([26, 8, 8])           # FeO2, neutral
    42
    >>> electron_count([26, 8, 8], charge=2)
    40
    """
    return sum(atomic_numbers) - charge


@dataclass(frozen=True, order=True)
class Assignment:
    """One oxidation-state reading that supports a Variant.

    Attributes
    ----------
    oxidation_state : int
        Oxidation state of the metal center for one consistent reading.
    configuration : str
        Shell occupancy label for that reading, such as ``"d6"`` or ``"s1"``.
        Transition metals use d counts; main-group centers use s counts.
    """

    oxidation_state: int
    configuration: str


@dataclass
class Variant:
    """One feasible (total charge, multiplicity) combination for the complex.

    A combination reachable under several oxidation-state hypotheses is reported
    once, with every supporting hypothesis recorded in `assignments`.

    Attributes
    ----------
    charge : int
        Total charge of the variant complex, taken from the requested charge
        window. The ligand set is assumed unchanged from the input, so any
        difference from the input charge is absorbed by the metal.
    multiplicity : int
        Spin multiplicity (unpaired electrons + 1). Its parity is forced by
        `charge`: an odd change in total charge changes the electron count by
        one and so flips multiplicity parity.
    is_input : bool
        True for the single variant equal to the input charge and multiplicity.
        Filter on this to keep only new combinations.
    assignments : list of Assignment
        Sorted supporting readings of the metal that make this combination
        reachable. `Assignment.configuration` is a shell occupancy label such as
        ``"d6"`` or ``"s1"``. Never empty, and never contains a d count above
        `MAX_D_ELECTRONS`. More than one entry means the oxidation state is
        genuinely ambiguous, not that the variant is duplicated.
    """
    charge: int
    multiplicity: int
    is_input: bool
    assignments: List[Assignment] = field(default_factory=list)


def enumerate_variants(
    metal: str,
    charge: int,
    multiplicity: int,
    charges: Sequence[int] = DEFAULT_CHARGES,
    max_multiplicity: int = 6,
    restrict_to_accessible: bool = True,
    neutral_electrons: Optional[int] = None,
) -> List[Variant]:
    """Enumerate feasible (charge, multiplicity) variants of a metal complex.

    Parameters
    ----------
    metal : str
        Element symbol of the metal center, capitalized (e.g. ``"Fe"``). Must be
        a key of `SPIN_STATES`.
    charge : int
        Total charge of the input complex, not the metal oxidation state.
    multiplicity : int
        Spin multiplicity of the input complex (unpaired electrons + 1).
    charges : sequence of int, optional
        Total charges to enumerate variants at. These are absolute charges, not
        offsets from `charge`, so the input charge need not appear among them.
        Default `DEFAULT_CHARGES`, ``(-1, 0, 1)``.
    max_multiplicity : int, optional
        Upper bound on reported multiplicity, applied whatever
        `restrict_to_accessible` is set to. Default 6.
    restrict_to_accessible : bool, optional
        If True (default), restrict both input hypotheses and target states to
        this element's entry in `ACCESSIBLE_OXIDATION_STATES`. If False, use every
        oxidation state tabulated in `SPIN_STATES`. Elements absent from the
        oxidation-state table are unrestricted either way.
    neutral_electrons : int, optional
        Electron count of the neutral molecule (sum of atomic numbers), used to
        check `charge` and `multiplicity` against `parity_consistent`. Only the
        input pair is tested: every enumerated variant shifts charge and
        multiplicity together, so all variants share the input's parity status
        and one check settles them. Default None, skipping the check.

    Returns
    -------
    list of Variant
        Feasible combinations, sorted by charge then multiplicity, one entry per
        distinct pair. Includes the input combination, flagged by
        `Variant.is_input`, whenever `charge` is in `charges`. May be empty.

    Raises
    ------
    KeyError
        If `metal` has no entry in `SPIN_STATES`.
    ValueError
        If `multiplicity` is not compatible with this element's accessible
        oxidation states while `restrict_to_accessible` is True, is not tabulated
        for any oxidation state of `metal`, or fails the `neutral_electrons`
        parity check.

    Examples
    --------
    >>> variants = enumerate_variants("Pd", 0, 1)
    >>> [(v.charge, v.multiplicity) for v in variants]
    [(0, 1), (0, 3)]

    Multiplicity 6 requires a half-filled d5 shell, which pins the oxidation
    state and narrows the result:

    >>> variants = enumerate_variants("Fe", 1, 6)
    >>> [(v.charge, v.multiplicity) for v in variants]
    [(-1, 2), (-1, 4), (0, 1), (0, 3), (0, 5), (1, 2), (1, 4), (1, 6)]

    A neutral 100-electron complex cannot be a doublet, and is rejected before
    any variant is built:

    >>> enumerate_variants("Fe", 0, 2, neutral_electrons=100)
    Traceback (most recent call last):
        ...
    ValueError: charge 0 with multiplicity 2 is parity inconsistent: a 100-electron species must have odd multiplicity

    Oxidising it by one electron makes the doublet the consistent choice:

    >>> variants = enumerate_variants("Fe", 1, 2, neutral_electrons=100)
    >>> [(v.charge, v.multiplicity) for v in variants]
    [(-1, 2), (-1, 4), (-1, 6), (0, 1), (0, 3), (0, 5), (1, 2), (1, 4), (1, 6)]
    """
    if metal not in SPIN_STATES:
        raise KeyError(f"no spin-state data for element {metal!r}")

    if neutral_electrons is not None and not parity_consistent(
            neutral_electrons, charge, multiplicity):
        total = neutral_electrons - charge
        raise ValueError(
            f"charge {charge} with multiplicity {multiplicity} is parity "
            f"inconsistent: a {total}-electron species must have "
            f"{'odd' if total % 2 == 0 else 'even'} multiplicity"
        )

    table = SPIN_STATES[metal]
    allowed_oxidation_states = (
        ACCESSIBLE_OXIDATION_STATES.get(metal)
        if restrict_to_accessible else None
    )

    def permitted(mult: int) -> bool:
        return mult <= max_multiplicity

    def oxidation_state_permitted(oxidation_state: int) -> bool:
        return (allowed_oxidation_states is None
                or oxidation_state in allowed_oxidation_states)

    matching_oxidation_states = [
        ox for ox in sorted(table)
        if multiplicity in table[ox][1] and oxidation_state_permitted(ox)
    ]
    if not matching_oxidation_states:
        if allowed_oxidation_states is not None:
            raise ValueError(
                f"multiplicity {multiplicity} is not compatible with the "
                f"accessible oxidation states for {metal} "
                f"({', '.join(str(ox) for ox in sorted(allowed_oxidation_states))}); "
                f"pass restrict_to_accessible=False to override"
            )
        raise ValueError(
            f"multiplicity {multiplicity} is not tabulated for any oxidation "
            f"state of {metal}"
        )

    shell_label = shell(metal)

    # (charge, multiplicity) -> supporting oxidation-state readings
    found: Dict[Tuple[int, int], Set[Assignment]] = {}
    for ox in matching_oxidation_states:
        for target_charge in charges:
            target_ox = ox + (target_charge - charge)
            if target_ox not in table or not oxidation_state_permitted(target_ox):
                continue
            electron_count, mults = table[target_ox]
            for mult in mults:
                if permitted(mult):
                    found.setdefault((target_charge, mult), set()).add(
                        Assignment(
                            oxidation_state=target_ox,
                            configuration=f"{shell_label}{electron_count}",
                        )
                    )

    variants = [
        Variant(
            charge=q,
            multiplicity=mult,
            is_input=(q == charge and mult == multiplicity),
            assignments=sorted(readings),
        )
        for (q, mult), readings in sorted(found.items())
    ]
    return variants

