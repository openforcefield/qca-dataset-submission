import os
from pprint import pformat
import traceback
from itertools import chain
from collections import defaultdict, Counter
from datetime import datetime

QCFRACTAL_URL = "https://api.qcarchive.molssi.org:443/"


DATASET_TYPES = {
    'dataset': 'singlepoint',
    'optimizationdataset': 'optimization',
    'torsiondrivedataset': 'torsiondrive'
}

def get_version_info():
    """
    Get the version info for the packages used to validate the submission.
    """
    import importlib
    import pandas as pd

    report = {}
    # list the core packages here
    packages = ["openff.qcsubmit", "openff.toolkit", "basis_set_exchange", "qcelemental"]
    for package in packages:
        module = importlib.import_module(package)
        report[package] = pd.Series({"version": module.__version__})

    # now try openeye else use rdkit
    try:
        import openeye

        report["openeye"] = pd.Series({"version": openeye.__version__})
    except ImportError:
        import rdkit

        report["rdkit"] = pd.Series({"version": rdkit.__version__})

    return pd.DataFrame(report).transpose()


def create_dataset(dataset_data):
    from openff.qcsubmit.datasets import BasicDataset, OptimizationDataset, TorsiondriveDataset

    datasets = {
        "dataset": BasicDataset,
        "optimizationdataset": OptimizationDataset,
        "torsiondrivedataset": TorsiondriveDataset,
    }

    if "type" in dataset_data:
        dataset_type = dataset_data["type"].lower()
    elif "dataset_type" in dataset_data:
        dataset_type = dataset_data["dataset_type"].lower()

    dataset_class = datasets.get(dataset_type, None)
    if dataset_class is not None:
        return dataset_class.parse_obj(dataset_data)
    else:
        raise RuntimeError(f"The dataset type {dataset_type} is not supported.")


class SubmittableBase:
    def __init__(self, submittable, submission,
                 priority=1, computetag='openff'):
        """Create new Submittable instance linking a submission dataset to its PR.

        Parameters
        ----------
        submittable : path-like
            Path to submission file.
        submission : Submission
            Submission instance corresponding to the dataset submission.
        repo : str
            Github repo where datasets are tracked.
        priority : int
            Priority to use for the dataset if set by method calls;
            one of 0, 1, or 2, in increasing-priority order.
        computetag : str
            Compute tag to use for the dataset if set by method calls;
            tasks with a given compute tag will only be computed by managers
            configured to service that tag.

        """
        self.submittable = submittable
        self.submission = submission
        self.item = submission.item
        self.priority = priority
        self.computetag = computetag


    def _parse_spec(self):
        spec = self._load_submittable()
        dataset_name = spec["dataset_name"]
        if "type" in spec:
            dataset_type = DATASET_TYPES[spec["type"].lower()]
        elif "dataset_type" in spec:
            dataset_type = DATASET_TYPES[spec["dataset_type"].lower()]
        dataset_specs = spec.get("qc_specifications", None)
        return dataset_name, dataset_type, dataset_specs

    def _load_submittable(self):
        from openff.qcsubmit.serializers import deserialize

        spec = deserialize(self.submittable)
        return spec

    def _get_qca_client(self):
        import qcportal as ptl

        client = ptl.PortalClient(
            address=QCFRACTAL_URL,
            username=os.environ["QCA_USER"],
            password=os.environ["QCA_KEY"]
        )
        return client

    def _get_meta(self):
        import pandas as pd

        datehr = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
        dataset_name, dataset_type, dataset_specs = self._parse_spec()

        meta = {
            "**Dataset Name**": dataset_name,
            "**Dataset Type**": dataset_type,
            "**UTC Datetime**": datehr,
        }

        return pd.DataFrame(pd.Series(meta, name=""))

    def _version_info_report(self):
        version = get_version_info()

        comment = f"""
        <details>
        <summary><b>QCSubmit</b> version information(<i>click to expand</i>)</summary>
        <!-- have to be followed by an empty line! -->

        {version.to_markdown()}
        </details>
        """

        return comment

    def execute_queued_submit(self, max_retries=3):
        """Submit, perhaps with some retry logic.

        """
        client = self._get_qca_client()

        # load dataset into QCSubmit class
        ds = self._load_submittable()
        dataset_qcs = create_dataset(ds)

        try:
            # submit to QCArchive
            output = self.submit(dataset_qcs, client)
            self._queued_submit_report(output, success=True)
        except:
            self._queued_submit_report(traceback.format_exc(), success=False)
            return {"new_state": "Queued for Submission"}
        else:
            return {"new_state": "Error Cycling"}

    def _queued_submit_report(self, output, success):
        success_text = "**SUCCESS**" if success else "**FAILED**"

        comment = f"""
        ## Lifecycle - QCSubmit Submission Report : {success_text}

        {self._get_meta().to_markdown()}
        
        Response from public QCArchive:

        ```
        {output}
        ```

        ----------
        {self._version_info_report()}

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.item.add_issue_comment(comment)

    def execute_errorcycle(self,
                           reset_errors=False,
                           set_priority=False,
                           set_computetag=False):
        """Obtain complete, incomplete, error stats for submittable and report.

        For suspected random errors, we perform restarts.

        If submittable complete, recommend state "Archived/Complete".

        """
        client = self._get_qca_client()

        dataset_name, dataset_type, dataset_specs = self._parse_spec()
        ds = client.get_dataset(dataset_type, dataset_name)

        if dataset_type == "torsiondrive":
            complete = self._errorcycle_torsiondrive(
                    ds, client, dataset_specs,
                    reset_errors=reset_errors, set_priority=set_priority,
                    set_computetag=set_computetag)

        elif dataset_type == "optimization":
            complete = self._errorcycle_dataset(
                    ds, client, dataset_specs,
                    self._errorcycle_optimization_report,
                    reset_errors=reset_errors, set_priority=set_priority,
                    set_computetag=set_computetag)

        elif dataset_type == "singlepoint":
            complete = self._errorcycle_dataset(
                    ds, client, dataset_specs,
                    self._errorcycle_dataset_report,
                    reset_errors=reset_errors, set_priority=set_priority,
                    set_computetag=set_computetag)

        if complete:
            return {"new_state": "Archived/Complete"}
        else:
            return {"new_state": "Error Cycling"}

    def comment_archived_complete(self):
        comment = f"""
        ## Lifecycle - Archived/Complete

        {self._get_meta().to_markdown()}

        **Dataset Complete!**

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.item.add_issue_comment(comment)

    @staticmethod
    def count_unique_error_messages(errors_in, pretty_print=False):
        errors = defaultdict(set)
    
        for id, error in errors_in.items():
            errors["\n".join([error[i] for i in ['error_type', 'error_message']])].add(id)
    
        errors = dict(errors)
    
        content = ""
        if pretty_print:
            for count, key, value in sorted([(len(value), key, value) for key, value in errors.items()], reverse=True):
                content += '-------------------------------------\n'
                content += f"count : {count}\n"
                content += '\n'
                content += f'{key}\n'
                content += '\n'
                content += 'ids : \n'
                content += f'{pformat(value, width=80, compact=True)}\n'
                content += '-------------------------------------\n'
            return content
        else:
            return errors

    def _errorcycle_torsiondrive(self, ds, client, dataset_specs,
            reset_errors=False, set_priority=False, set_computetag=False):
        import pandas as pd
        from qcportal.record_models import RecordStatusEnum

        if dataset_specs is None:
            dataset_specs = ds.specification_names

        df_status = self._errorcycle_get_status(ds, dataset_specs)

        if reset_errors:
            erred_rec_ids = []
            erred_opts = {}
            status_counts = {}
            for ds_spec in dataset_specs:
                recs = ds.iterate_records(
                        specification_names=[ds_spec], 
                        #status='error'
                        )

                # build up optimization statuses and errors, if present
                erred_opts[ds_spec] = []
                status_counts[ds_spec] = Counter({status.value.upper(): 0 for status in list(RecordStatusEnum)})
                for entry, spec, rec in recs:
                    if rec.status == 'error':
                        erred_rec_ids.append(rec.id)
                    for opt in chain.from_iterable(rec.optimizations.values()):
                        status_counts[ds_spec][opt.status.value.upper()] += 1

                        if opt.status == 'error':
                            erred_opts[ds_spec].append((opt.id, opt.error))

            # create status counts dataframe
            df_opt_status = pd.DataFrame(status_counts).transpose()
            df_opt_status = df_opt_status[['COMPLETE', 'RUNNING', 'WAITING', 'ERROR', 'CANCELLED', 'INVALID', 'DELETED']]
            df_opt_status.index.name = 'specification'

            # aggregate all errors to get single set of counts for error messages
            errors = {}
            for ds_spec in erred_opts:
                errors.update({r[0]: r[1] for r in erred_opts[ds_spec]})

            error_counts = self.count_unique_error_messages(errors, pretty_print=True)

            self._errorcycle_torsiondrive_report(df_status, df_opt_status, error_counts)

        if df_status[["WAITING", "RUNNING", "ERROR"]].sum().sum() == 0:
            complete = True
        else:
            if reset_errors:
                client.reset_records(erred_rec_ids)
            if set_priority:
                ds.modify_records(specification_names=list(dataset_specs),
                                  new_priority=self.priority)
            if set_computetag:
                ds.modify_records(specification_names=list(dataset_specs),
                                  new_tag=self.computetag)
            complete = False

        return complete

    def _errorcycle_torsiondrive_report(self, df_tdr, df_tdr_opt, opt_error_counts):

        if len(opt_error_counts) > 60000:
            opt_error_counts = opt_error_counts[:60000]
            opt_error_counts += "\n--- Too many errors; truncated here ---\n"

        comment = f"""
        ## Lifecycle - Error Cycling Report

        {self._get_meta().to_markdown()}

        All errored tasks and services will be restarted.
        Errored states prior to restart reported below.

        ### `TorsionDriveRecord` current status

        {df_tdr.to_markdown()}

        ### `OptimizationRecord` current status

        {df_tdr_opt.to_markdown()}

        #### `OptimizationRecord` Error Tracebacks:

        <details>
        <summary><b>Tracebacks</b> (<i>click to expand</i>)</summary>
        <!-- have to be followed by an empty line! -->

        ```
        {opt_error_counts}
        ```
        </details>

        ----------
        {self._version_info_report()}

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.item.add_issue_comment(comment)

    def _errorcycle_get_status(self, ds, dataset_specs):
        import pandas as pd
        from qcportal.record_models import RecordStatusEnum

        if dataset_specs is None:
            dataset_specs = ds.specification_names
        
        status = ds.status()
        status_ = {key: {status.value.upper(): counts.get(status, 0)
                         for status in list(RecordStatusEnum)}
                   for key, counts in status.items() if key in dataset_specs.keys()}

        df = pd.DataFrame(status_).transpose()
        df = df[['COMPLETE', 'RUNNING', 'WAITING', 'ERROR', 'CANCELLED', 'INVALID', 'DELETED']]
        df.index.name = 'specification'

        return df

    def _errorcycle_optimization_report(self, df_status, opt_error_counts):

        if len(opt_error_counts) > 60000:
            opt_error_counts = opt_error_counts[:60000]
            opt_error_counts += "\n--- Too many errors; truncated here ---\n"

        comment = f"""
        ## Lifecycle - Error Cycling Report

        {self._get_meta().to_markdown()}

        All errored tasks will be restarted.
        Errored states prior to restart reported below.

        ### `OptimizationRecord` current status

        {df_status.to_markdown()}

        #### `OptimizationRecord` Error Tracebacks:

        <details>
        <summary><b>Tracebacks</b> (<i>click to expand</i>)</summary>
        <!-- have to be followed by an empty line! -->

        ```
        {opt_error_counts}
        ```
        </details>

        ----------
        {self._version_info_report()}

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.item.add_issue_comment(comment)

    def _errorcycle_dataset(self, ds, client, dataset_specs, report_method,
            reset_errors=False, set_priority=False, set_computetag=False):

        if dataset_specs is None:
            dataset_specs = ds.specification_names

        df_status = self._errorcycle_get_status(ds, dataset_specs)

        if reset_errors:
            erred_recs = ds.iterate_records(
                    specification_names=list(dataset_specs), 
                    status='error')

            errors = {r.id: r.error for entry, spec, r in erred_recs}
            error_counts = self.count_unique_error_messages(errors, pretty_print=True)

            report_method(df_status, error_counts)

        if df_status[["WAITING", "RUNNING", "ERROR"]].sum().sum() == 0:
            complete = True
        else:
            if reset_errors:
                client.reset_records(list(errors))
            if set_priority:
                ds.modify_records(specification_names=list(dataset_specs),
                                  new_priority=self.priority)
            if set_computetag:
                ds.modify_records(specification_names=list(dataset_specs),
                                  new_tag=self.computetag)
            complete = False

        return complete

    def _errorcycle_dataset_report(self, df_res, res_error_counts):

        if len(res_error_counts) > 60000:
            res_error_counts = res_error_counts[:60000]
            res_error_counts += "\n--- Too many errors; truncated here ---\n"

        comment = f"""
        ## Lifecycle - Error Cycling Report

        {self._get_meta().to_markdown()}

        All errored tasks will be restarted.
        Errored states prior to restart reported below.

        ### `ResultRecord` current status

        {df_res.to_markdown()}

        #### `ResultRecord` Error Tracebacks:

        <details>
        <summary><b>Tracebacks</b> (<i>click to expand</i>)</summary>
        <!-- have to be followed by an empty line! -->

        ```
        {res_error_counts}
        ```
        </details>

        ----------
        {self._version_info_report()}

        """

        # postprocess due to raw spacing above
        comment = "\n".join([substr.strip() for substr in comment.split("\n")])

        # submit comment
        self.item.add_issue_comment(comment)

    def submit(self, dataset_qcs, client):
        return dataset_qcs.submit(client=client, ignore_errors=True)


class DataSet(SubmittableBase):
    """A dataset submitted to QCArchive.

    A dataset has a lifecycle with well-defined states.
    The state of a dataset is the state of its submission PR.

    """
    ...


class Compute(SubmittableBase):
    """Supplemental compute submitted to QCArchive.

    """
    ...
