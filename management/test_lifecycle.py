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
            "compute-pr393_mw-300-600-900",
            ([300.0, 600.0, 900.0], "compute-pr393"),
        ),
        ("compute-pr393", ([], "compute-pr393")),
    ],
)
def test_parse_tags(input_tag, want):
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
            [6, 220, 234, 439],
        ),
        (
            "OpenFF Sulfur Optimization Training Coverage Supplement v1.0",
            "Optimization",
            [100.0, 200.0, 300.0],
            [6, 220, 234, 439],
        ),
        (
            "OpenFF Alkane Torsion Drives v1.0",
            "torsiondrive",
            [50.0, 75.0, 100.0],
            [2, 30, 68, 92],
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

    assert list(sorted(lens)) == want


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

        def modify_records(self, record_ids, new_tag):
            self.calls.append((len(record_ids), new_tag))

    client = PortalClient("https://api.qcarchive.molssi.org:443/")
    ds = client.get_dataset(dstype, dsname)

    dummy = DummyClient()
    update_compute_tags(dummy, ds, list(), tag, include_complete=True)

    # calls depends on dict iteration order, so sort the output by number of
    # record_ids in each bin
    assert list(sorted(dummy.calls, key=lambda x: x[0])) == want
