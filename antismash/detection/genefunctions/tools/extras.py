# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Classifies genes based on a set of unrelated profiles
"""

from dataclasses import dataclass
from typing import Any, Iterable, Optional, Self

from antismash.common import hmmer, json, path
from antismash.common.secmet.features import CDSFeature
from antismash.common.secmet.qualifiers.gene_functions import GeneFunction
from antismash.config import ConfigType

from .core import HMMFunctionResults as Results, HMMHit, Tool, scan_profiles_for_functions

TOOL_NAME = "extras"
DATABASE = path.get_full_path(__file__, "data", "extras.hmm")
ENTRY_DATA = path.get_full_path(__file__, "data", "extras.json")


@dataclass
class Entry:
    identifier: str
    description: str
    cutoff: float
    function: GeneFunction
    source: str

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        data["function"] = GeneFunction.from_string(data["function"])
        return cls(**data)


def check_prereqs(options: ConfigType) -> list[str]:
    """ Check all prerequisites are satisfied
    """
    failure_messages = []

    for binary_name in ['hmmscan', 'hmmpress']:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r}")
    failure_messages.extend(prepare_data(logging_only=True))
    return failure_messages


def classify(cds_features: Iterable[CDSFeature], _options: ConfigType) -> Results:
    """ Finds possible classifications for the provided CDS features.

        Arguments:
            cds_features: a list of CDSFeatures to classify
            options: the current run options

        Returns:
            a results instance with best hits and classification for each CDS
    """
    metadata = _load_metadata()
    hits: dict[str, HMMHit] = scan_profiles_for_functions(cds_features, DATABASE, hmmscan_opts=["-E", "1E-16"])
    filtered = {}
    cds_name_to_function = {}
    for cds_name, hit in hits.items():
        entry = metadata[hit.hit_id]
        if hit.bitscore < entry.cutoff:
            continue
        cds_name_to_function[cds_name] = entry.function
        hit.description = entry.description
        filtered[cds_name] = hit

    return Results(tool=TOOL_NAME, best_hits=filtered, function_mapping=cds_name_to_function)


def _load_metadata(metadata: str = ENTRY_DATA) -> dict[str, Entry]:
    with open(metadata, encoding="utf-8") as handle:
        raw = json.load(handle)
    entries = [Entry.from_json(entry) for entry in raw["entries"]]
    return {entry.identifier: entry for entry in entries}


def prepare_data(logging_only: bool = False) -> list[str]:
    """ Ensures tool data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failures = []
    failures.extend(hmmer.ensure_database_pressed(DATABASE, return_not_raise=logging_only))
    return failures


def regenerate_results(data: Optional[dict[str, Any]]) -> Optional[Results]:
    if not data:
        return None
    return Results.from_json(data)


TOOL = Tool(
    name=TOOL_NAME,
    check_prereqs=check_prereqs,
    classify=classify,
    prepare_data=prepare_data,
)
