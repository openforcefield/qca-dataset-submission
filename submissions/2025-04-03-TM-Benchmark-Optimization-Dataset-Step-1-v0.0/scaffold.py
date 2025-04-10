"""Functions to export and import elect information from QCFractal datasets

This module will allow users to write and read a json file that can be reformed into an initial dataset
**before** submission. This means that records are not present.
"""

import json
import bz2

from qcportal.serialization import encode_to_json


def to_json(ds, filename="scaffold.json", indent=4, compress=False):
    """Export a QCFractal dataset to a json file.

    Can be imported with :func:`fom_json` to make a new dataset.

    Args:
        ds (qcportal.*Dataset): QCFractal dataset
        filename (str, optional): Filename/path to store output json file. Defaults to "scaffold.json".
        indent (int, optional): Level of indent for the output json file. Defaults to 4.
        compress (bool, optional): If True, will compress to bz2. Defaults to False.
    """

    inputs = [
        "dataset_type",
        "name",
        "description",
        "tagline",
        "tags",
        "group",
        "provenance",
        "visibility",
        "default_tag",
        "default_priority",
        "metadata",
        "extras",
        "owner_group",
    ]  # Inputs for client.add_dataset(
    metadata = {key: value for key, value in ds.dict().items() if key in inputs}
    d = {
        "metadata": metadata,
        "entries": {entry.name: entry for entry in ds.iterate_entries()},
        "specifications": ds.specifications,
    }
    d_serializable = encode_to_json(d)
    
    if compress:
        with bz2.open(filename+".bz2", "wt", encoding="utf-8") as f:
            json.dump(d_serializable, f, ensure_ascii=False, indent=indent)
    else:
        with open(filename, "w") as f:
            json.dump(d_serializable, f, indent=indent)