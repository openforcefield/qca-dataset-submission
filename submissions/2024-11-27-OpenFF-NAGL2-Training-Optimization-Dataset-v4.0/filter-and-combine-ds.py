import json
import logging
import typing

import click

from openff.qcsubmit.results.filters import SinglepointRecordFilter

# suppress stereochemistry warnings
logging.getLogger("openff").setLevel(logging.ERROR)

if typing.TYPE_CHECKING:
    from qcportal.optimization.record_models import OptimizationRecord
    from openff.qcsubmit.results import TorsionDriveResultCollection, OptimizationResultCollection
    from openff.toolkit import ForceField, Molecule

@click.group()
def cli():
    pass

def download_and_filter_opt_data(
    opt_datasets: typing.List[str],
    verbose: bool = False,
) -> "OptimizationResultCollection":
    """Download and filter optimization datasets."""

    from qcportal import PortalClient
    from qcportal.record_models import RecordStatusEnum
    from openff.qcsubmit.results import OptimizationResultCollection
    from openff.qcsubmit.results.filters import (
        ConnectivityFilter,
        RecordStatusFilter,
        UnperceivableStereoFilter,
        ElementFilter,
        ConformerRMSDFilter,
    )

    # download dataset(s) from QCArchive
    client = PortalClient(address="https://api.qcarchive.molssi.org:443/")
    dataset = OptimizationResultCollection.from_server(
        client=client,
        datasets=opt_datasets,
        spec_name="default",
    )
    if verbose:
        print(f"Number of entries before filtering: {dataset.n_results}")
    

    # Filter out other unsuitable entries
    dataset = dataset.filter(
        RecordStatusFilter(status=RecordStatusEnum.complete),
        ConnectivityFilter(tolerance=1.2),
        UnperceivableStereoFilter(),
        #ConformerRMSDFilter(max_conformers=max_opt_conformers)

    )
    return dataset

@cli.command("download-opt")
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    help="The path to write the dataset to. Should be a JSON",
)
@click.option(
    "--confrmsd",
    "confrmsd_path",
    required=True,
    type=click.Path(exists=False, dir_okay=False, file_okay=True),
    help="The path to write the dataset to. Should be a JSON",
)
@click.option(
    "--core-opt-dataset",
    "core_opt_datasets",
    multiple=True,
    required=True,
    type=str,
    help=(
        "The name of an optimization dataset to download. "
        "These will have iodine molecules filtered out."
    ),
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    show_default=True,
    help="Whether to print out additional information.",
)
def download_opt_data(
    output_path: str,
    confrmsd_path: str,
    core_opt_datasets: typing.List[str],
    verbose: bool = True,
):
    """Download and filter optimization datasets.

    \f
    Parameters
    ----------
    core_opt_datasets
        The core optimization datasets to download.
    """
    from openff.qcsubmit.results import OptimizationResultCollection
    from openff.toolkit import ForceField
    from openff.qcsubmit.results.filters import ConformerRMSDFilter
    

    # suppress stereochemistry warnings
    logging.getLogger("openff").setLevel(logging.ERROR)

    # download and filter core dataset(s)
    core_dataset = download_and_filter_opt_data(
        core_opt_datasets, verbose=verbose
    )
    if verbose:
        print(f"Number of filtered core entries: {core_dataset.n_results}")

    with open(output_path, "w") as file:
        file.write(core_dataset.json(indent=2))
    if verbose:
        print(f"Saved to {output_path}")

    # apply conformer RMSD filter
    filtered_for_confrmsd = core_dataset.filter(ConformerRMSDFilter(rmsd_tolerance=0.01))
    if verbose:
        print(f"Number of entries after conformer RMSD check: {filtered_for_confrmsd.n_results}")

    with open(confrmsd_path, "w") as file:
        file.write(filtered_for_confrmsd.json(indent=2))
    if verbose:
        print(f"Saved to {confrmsd_path}")


if __name__ == "__main__":
    cli()
