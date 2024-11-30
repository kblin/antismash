# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A module to find and mark gene functions in records with a variety of tools
"""

from collections import defaultdict
from dataclasses import dataclass, field
from itertools import chain
import logging
import os
from typing import Any, Dict, Iterator, List, Optional, Self

from markupsafe import Markup

from antismash.common import hmmer, module_results, path
from antismash.common.pfamdb import check_db
from antismash.common.secmet import Record
from antismash.common.secmet.qualifiers.gene_functions import ECGroup
from antismash.config import ConfigType, get_config
from antismash.config.args import ModuleArgs
from antismash.detection import DetectionStage

from antismash.common.html_renderer import HTMLSections, FileTemplate
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer
from antismash.detection.genefunctions.tools.core import HMMHit
from antismash.detection.genefunctions.tools.halogenases.data_structures import HalogenaseHit

from .tools import (
    FunctionResults,
    Hit,
    Tool,
    extras,
    halogenases,
    mite,
    resistance,
    smcogs,
)

from .tools.extras import Results as ExtrasResults
from .tools.halogenases import HalogenaseResults
from .tools.mite import MiteHit, Results as MiteResults
from .tools.smcogs import Results as SmcogsResults
from .tools.resistance import Results as ResistanceResults

NAME = "genefunctions"
SHORT_DESCRIPTION = "Gene function annotations"
DETECTION_STAGE = DetectionStage.PER_AREA

TOOLS = [
    extras.TOOL,
    halogenases.TOOL,
    mite.TOOL,
    smcogs.TOOL,
    resistance.TOOL,
]


@dataclass
class ToolResults:
    """ A container for all tools generating gene function results """
    # members are named by tool to simplify downstream typing, avoiding a great
    # deal of extra type casting
    extras: Optional[ExtrasResults] = None
    halogenases: Optional[HalogenaseResults] = None
    mite: Optional[MiteResults] = None
    resist: Optional[ResistanceResults] = None
    smcogs: Optional[SmcogsResults] = None

    def __iter__(self) -> Iterator[FunctionResults[Any]]:
        if self.extras:
            yield self.extras
        if self.mite:
            yield self.mite
        if self.resist:
            yield self.resist
        if self.smcogs:
            yield self.smcogs
        if self.halogenases:
            yield self.halogenases

    def add_tool_results(self, tool: Tool, results: FunctionResults[Hit]) -> None:
        """ Adds the given tool's results to the instance.

            Arguments:
                tool: the Tool from which the results were generated
                results: the results to be added
        """
        if not hasattr(self, tool.name.lower()):
            raise NotImplementedError(f"Results storage has no handling for {tool.name}")
        if getattr(self, tool.name.lower()) is not None:
            raise ValueError(f"Gene function results already exist for tool: {results.tool}")
        if tool is smcogs.TOOL:
            assert isinstance(results, smcogs.Results)
            self.smcogs = results
            return
        if tool is resistance.TOOL:
            assert isinstance(results, resistance.Results)
            self.resist = results
            return
        if tool is extras.TOOL:
            assert isinstance(results, ExtrasResults)
            self.extras = results
            return
        if tool is mite.TOOL:
            assert isinstance(results, mite.Results)
            self.mite = results
            return
        if tool is halogenases.TOOL:
            assert isinstance(results, halogenases.HalogenaseResults)
            self.halogenases = results
            return

        raise NotImplementedError(f"Results storage has no handling for {tool.name}")

    @classmethod
    def from_json(cls, data: dict[str, dict[str, Any]]) -> Self:
        """ Reconstructs an instance from the given data """
        return cls(
            extras=extras.regenerate_results(data.get(extras.TOOL.name.lower())),
            mite=mite.regenerate_results(data.get(mite.TOOL.name.lower())),
            smcogs=smcogs.regenerate_results(data.get(smcogs.TOOL.name.lower())),
            resist=resistance.regenerate_results(data.get(resistance.TOOL.name.lower())),
            halogenases=halogenases.regenerate_results(data.get(halogenases.TOOL.name.lower())),
        )


class AllFunctionResults(module_results.DetectionResults):
    """ A collection of results for a variety of gene function detections """
    schema_version = 2

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self.tool_results = ToolResults()

    def add_tool_results(self, tool: Tool, results: FunctionResults[Hit]) -> None:
        """ Add results for a tool, tool name must be unique and user-friendly """
        self.tool_results.add_tool_results(tool, results)

    def add_to_record(self, record: Record) -> None:
        for result in self.tool_results:
            result.add_to_record(record)

    @staticmethod
    def from_json(json: Dict[str, Any], record: Record) -> Optional["AllFunctionResults"]:
        if json.get("schema_version") != AllFunctionResults.schema_version:
            logging.debug("Schema version mismatch, discarding %s results", NAME)
            return None
        if record.id != json.get("record_id"):
            logging.debug("Record ID mismatch, discarding %s results", NAME)
            return None
        all_results = AllFunctionResults(record.id)
        all_results.tool_results = ToolResults.from_json(json["tools"])
        return all_results

    def to_json(self) -> Dict[str, Any]:
        return {"schema_version": self.schema_version,
                "record_id": self.record_id,
                "tools": self.tool_results,
                }


def get_arguments() -> ModuleArgs:
    """ Build and return arguments. No extra options beyond a switch to enable """
    args = ModuleArgs('Additional analysis', 'genefunctions', enabled_by_default=True)
    for tool in TOOLS:
        tool.add_arguments(args)
    return args


def check_options(options: ConfigType) -> list[str]:
    """ Checks options for conflicts or issues.
    """
    return list(chain.from_iterable(tool.check_options(options) for tool in TOOLS))


def prepare_data(logging_only: bool = False) -> List[str]:
    """ Ensures packaged data is fully prepared

        Arguments:
            logging_only: whether to return error messages instead of raising exceptions

        Returns:
            a list of error messages (only if logging_only is True)
    """
    failures: list[str] = []
    failures.extend(chain.from_iterable(tool.prepare_data(logging_only) for tool in TOOLS))
    return failures


def check_prereqs(options: ConfigType) -> List[str]:
    """Check for prerequisites
    """
    failure_messages = []

    for binary_name in ['hmmscan', 'hmmpress']:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r}")

    database = os.path.join(options.database_dir, 'resfam', 'Resfams.hmm')
    if "mounted_at_runtime" not in database and path.locate_file(database) is None:
        failure_messages.append(f"Failed to locate Resfam database in {database!r}")

    # normally a module would check the data is prepared here,
    # but only tools currently have data, so prepare_data here would only duplicate those checks

    generator = (tool.check_prereqs(options) for tool in TOOLS if tool.check_options is not None)
    failure_messages.extend(chain.from_iterable(generator))
    return failure_messages


def is_enabled(options: ConfigType) -> bool:
    """ Should the module be run with these options """
    return not options.minimal or options.genefunctions_enabled


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[AllFunctionResults]:
    """ Regenerate the previous results from JSON format. """
    if not previous:
        return None
    return AllFunctionResults.from_json(previous, record)


def run_on_record(record: Record, results: AllFunctionResults, options: ConfigType) -> AllFunctionResults:
    """ Run the detection, unless the previous results apply to the given record """
    if isinstance(results, AllFunctionResults) and results.record_id == record.id:
        return results
    results = AllFunctionResults(record.id)
    cds_features = record.get_cds_features_within_regions()
    for tool in TOOLS:
        results.add_tool_results(tool, tool.classify(cds_features, options))
    return results


@dataclass(order=True)
class TailoringEntry:
    """ An entry in the Tailoring Enzymes table """
    name: str
    smcogs: HMMHit | None = field(default=None, kw_only=True)
    mite: MiteHit | None = field(default=None, kw_only=True)
    halogenase: HalogenaseHit | None = field(default=None, kw_only=True)
    extra: HMMHit | None = field(default=None, kw_only=True)
    url: str | None = field(default=None, kw_only=True)

    def add_tool_results(self, results: FunctionResults[Any]) -> None:
        """ Add the hit info for a tool """
        if isinstance(results, SmcogsResults):
            self.smcogs = results.best_hits.get(self.name)
        elif isinstance(results, MiteResults):
            self.mite = results.best_hits.get(self.name)
            self.url = results.url
        elif isinstance(results, HalogenaseResults):
            self.halogenase = results.best_hits.get(self.name)
        elif isinstance(results, ExtrasResults):
            self.extra = results.best_hits.get(self.name)


    @property
    def description(self) -> str:
        """ Get the most specific description of the entry """
        if not any((self.smcogs, self.mite, self.halogenase, self.extra)):
            raise ValueError("Need at least one tailoring descriptor set")

        if self.halogenase:
            return self.halogenase.description
        if self.mite:
            return self.mite.description
        if self.extra:
            return self.extra.description

        assert self.smcogs  # shut up, mypy
        return self.smcogs.description

    def __hash__(self) -> int:
        return hash(self.name)

    def __iter__(self) -> Iterator[Markup]:
        if self.halogenase:
            yield Markup(f"<dt>Halogenase</dt><dd>{self.halogenase.get_html_fragment(hide_id=True)}</dd>")
        if self.mite:
            metadata = {"url": self.url}
            yield Markup(f"<dt>MITE</dt><dd>{self.mite.get_html_fragment(metadata, hide_id=True)}</dd>")
        if self.extra:
            yield Markup(f"<dt>Extra</dt><dd>{self.extra.get_html_fragment(hide_id=True)}</dd>")
        if self.smcogs:
            yield Markup(f"<dt>smCoG</dt><dd>{self.smcogs.get_html_fragment(hide_id=True)}</dd>")


def generate_html(region_layer: RegionLayer, results: AllFunctionResults,
                  record_layer: RecordLayer, options_layer: OptionsLayer) -> HTMLSections:
    """ Generate the details panel HTML with results from the tailoring module """

    name_to_entry: dict[str, TailoringEntry] = {}
    entries_by_group: dict[ECGroup, set[TailoringEntry]] = defaultdict(set)


    for tool in results.tool_results:
        for cds, mappings in tool.group_mapping.items():
            if cds not in name_to_entry:
                name_to_entry[cds] = TailoringEntry(cds)
            entry = name_to_entry[cds]
            entry.add_tool_results(tool)
            for mapping in mappings:
                entries_by_group[mapping].add(entry)

    sorted_entries = {group: sorted(entries) for group, entries in entries_by_group.items()}

    # description of EC numbers and
    html = HTMLSections("tailoring")
    template = FileTemplate(path.get_full_path(__file__, "templates", "tailoring.html"))
    section = template.render(record = record_layer, region=region_layer,
                              results = results,
                              entries = sorted_entries,
                              tooltip = "placeholder tooltip")
    html.add_detail_section("Tailoring", section, "gene-function-details")
    return html

def will_handle(products: list[str], _product_categories: set[str]) -> bool:
    """ Returns true if one or more relevant products or product categories are present """
    return True
