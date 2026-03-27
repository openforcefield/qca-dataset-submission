import sys
from unittest import mock

import pytest

# avoid requiring PyGitHub to run tests
sys.modules["github"] = mock.MagicMock()


@pytest.mark.parametrize(
    "input_tag,want",
    [
        ("compute-openff_mw-300-600", ([300.0, 600.0], "compute-openff")),
        (
            "compute-arbitrary-compute-tag_mw-300-600-900",
            ([300.0, 600.0, 900.0], "compute-arbitrary-compute-tag"),
        ),
        (
            "compute-arbitrary-compute-tag-without-mw-suffix",
            ([], "compute-arbitrary-compute-tag-without-mw-suffix"),
        ),
    ],
)
def test_parse_tags(input_tag, want):
    """Test that ``parse_tags`` parses arbitrary compute tags ending with
    ``_mw[-###]+`` (matching the ``lifecycle.SPLIT_TAG`` regex). This function
    is not intended to be called with non-matching tags (as in the third test
    entry), but it should also handle this gracefully, returning an empty
    sequence of bins and the original tag.
    """
    from lifecycle import parse_tags

    got = parse_tags(input_tag)
    assert got == want


@pytest.mark.parametrize(
    "dsname,dstype,bins,want",
    [
        (
            "OpenFF Sulfur Hessian Training Coverage Supplement v1.0",
            "singlepoint",
            [100.0, 200.0, 300.0],
            {0: 6, 1: 439, 2: 234, 3: 220},
        ),
        (
            "OpenFF Sulfur Optimization Training Coverage Supplement v1.0",
            "Optimization",
            [100.0, 200.0, 300.0],
            {0: 6, 1: 439, 2: 234, 3: 220},
        ),
        (
            "OpenFF Alkane Torsion Drives v1.0",
            "torsiondrive",
            [50.0, 75.0, 100.0],
            {0: 2, 1: 30, 2: 92, 3: 68},
        ),
    ],
)
def test_partition_records(dsname, dstype, bins, want):
    from qcportal import PortalClient

    from lifecycle import partition_records

    client = PortalClient("https://api.qcarchive.molssi.org:443/")

    ds = client.get_dataset(dstype, dsname)

    records = partition_records(ds, bins, include_complete=True)

    lens = [len(v) for v in records.values()]
    assert sum(lens) == len(list(ds.iterate_entries()))

    got = {bin_id: len(recs) for bin_id, recs in records.items()}
    assert got == want


@pytest.mark.parametrize(
    "dsname,dstype,tag,want",
    [
        (
            "OpenFF Sulfur Hessian Training Coverage Supplement v1.0",
            "singlepoint",
            "compute-openff_mw-100-200-300",
            [
                (6, "compute-openff-100"),
                (220, "compute-openff-large"),
                (234, "compute-openff-300"),
                (439, "compute-openff-200"),
            ],
        ),
    ],
)
def test_update_compute_tags(dsname, dstype, tag, want):
    from qcportal import PortalClient

    from lifecycle import update_compute_tags

    class DummyClient:
        def __init__(self):
            self.calls = list()

        def (self, record_ids, new_tag):
            self.calls.append((len(record_ids), new_tag))

    client = PortalClient("https://api.qcarchive.molssi.org:443/")
    ds = client.get_dataset(dstype, dsname)

    dummy = DummyClient()
    update_compute_tags(dummy, ds, list(), tag, include_complete=True)

    # calls depends on dict iteration order, so sort the output by number of
    # record_ids in each bin
    assert list(sorted(dummy.calls, key=lambda x: x[0])) == want
