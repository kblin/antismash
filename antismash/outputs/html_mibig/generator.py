# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Responsible for creating the single web page results """

import json
import string
import os
import re
from typing import Any, Dict, List, Tuple, Union, cast

from antismash.common import html_renderer, path, module_results
from antismash.common.html_renderer import FileTemplate, HTMLSections
from antismash.common.json import JSONOrf
from antismash.common.layers import RecordLayer, OptionsLayer
from antismash.common.secmet import Record
from antismash.custom_typing import AntismashModule
from antismash.detection import nrps_pks_domains
from antismash.config import ConfigType
from antismash.outputs.html import js
from antismash.outputs.html.generator import (
    find_plugins_for_cluster,
    generate_searchgtr_htmls,
    write_regions_js,
)


def build_json_data(records: List[Record], results: List[Dict[str, module_results.ModuleResults]],
                    options: ConfigType) -> Tuple[List[Dict[str, Any]], List[Dict[str, Union[str, List[JSONOrf]]]]]:
    """ Builds JSON versions of records and domains for use in drawing SVGs with
        javascript.

        Arguments:
            records: a list of Records to convert
            results: a dictionary mapping record id to a list of ModuleResults to convert
            options: antiSMASH options

        Returns:
            a tuple of
                a list of JSON-friendly dicts representing records
                a list of JSON-friendly dicts representing domains
    """

    from antismash import get_all_modules  # TODO break circular dependency
    js_records = js.convert_records(records, results, options)

    js_domains = []

    for i, record in enumerate(records):
        json_record = js_records[i]
        # replace antismash cds_detail with mibig's one
        cds_annotations = results[i]["antismash.detection.mibig"].data["cluster"]["genes"]["annotations"]
        update_cds_description(json_record, cds_annotations)

        json_record['seq_id'] = "".join(char for char in json_record['seq_id'] if char in string.printable)
        for region, json_region in zip(record.get_regions(), json_record['regions']):
            handlers = find_plugins_for_cluster(get_all_modules(), json_region)
            if nrps_pks_domains not in handlers and nrps_pks_domains.domain_drawing.has_domain_details(region):
                handlers.append(cast(AntismashModule, nrps_pks_domains))
            for handler in handlers:
                # if there's no results for the module, don't let it try
                if handler.__name__ not in results[i]:
                    continue
                if "generate_js_domains" in dir(handler):
                    domains_by_region = handler.generate_js_domains(region, record)
                    if domains_by_region:
                        js_domains.append(domains_by_region)

    return js_records, js_domains


def generate_html_sections(records: List[RecordLayer], results: Dict[str, Dict[str, module_results.ModuleResults]],
                           options: ConfigType) -> Dict[str, Dict[int, List[HTMLSections]]]:
    """ Generates a mapping of record->region->HTMLSections for each record, region and module

        Arguments:
            records: a list of RecordLayers to pass through to the modules
            results: a dictionary mapping record name to
                        a dictionary mapping each module name to its results object
            options: the current antiSMASH config

        Returns:
            a dictionary mapping record id to
                a dictionary mapping region number to
                    a list of HTMLSections, one for each module
    """
    details = {}
    for record in records:
        record_details = {}
        record_result = results[record.id]
        for region in record.regions:
            # work around mibig module not creating protoclusters with the expected types
            assert len(region.subregions) == 1 and region.subregions[0].tool == "mibig"
            if nrps_pks_domains.domain_drawing.has_domain_details(region.region_feature):
                region.handlers.append(cast(AntismashModule, nrps_pks_domains))

            sections = []
            for handler in region.handlers:
                if handler.will_handle(region.products) or handler is nrps_pks_domains:
                    handler_results = record_result.get(handler.__name__)
                    if handler_results is None:
                        continue
                    sections.append(handler.generate_html(region, handler_results, record, options))
            record_details[region.get_region_number()] = sections
        details[record.id] = record_details
    return details


def generate_webpage(records: List[Record], results: List[Dict[str, module_results.ModuleResults]],
                     options: ConfigType) -> None:
    """ Generates and writes the HTML itself """

    generate_searchgtr_htmls(records, options)
    json_records, js_domains = build_json_data(records, results, options)
    write_regions_js(json_records, options.output_dir, js_domains)

    with open(os.path.join(options.output_dir, 'index.html'), 'w') as result_file:
        template = FileTemplate(path.get_full_path(__file__, "templates", "overview.html"))

        options_layer = OptionsLayer(options)
        record_layers_with_regions = []
        record_layers_without_regions = []
        results_by_record_id = {}  # type: Dict[str, Dict[str, module_results.ModuleResults]]
        for record, record_results in zip(records, results):
            if record.get_regions():
                record_layers_with_regions.append(RecordLayer(record, None, options_layer))
            else:
                record_layers_without_regions.append(RecordLayer(record, None, options_layer))
            results_by_record_id[record.id] = record_results

        regions_written = sum(len(record.get_regions()) for record in records)
        job_id = os.path.basename(options.output_dir)

        mibig_id = os.path.splitext(os.path.basename(options.mibig_json))[0]
        annotation_filename = "{}.json".format(mibig_id)
        page_title = mibig_id

        html_sections = generate_html_sections(record_layers_with_regions, results_by_record_id, options)

        svg_tooltip = ("Shows the layout of the region, marking coding sequences and areas of interest. "
                       "Clicking a gene will select it and show any relevant details. "
                       "Clicking an area feature (e.g. a candidate cluster) will select all coding "
                       "sequences within that area. Double clicking an area feature will zoom to that area. "
                       "Multiple genes and area features can be selected by clicking them while holding the Ctrl key."
                       )

        aux = template.render(records=record_layers_with_regions, options=options_layer,
                              version=options.version, extra_data=js_domains,
                              regions_written=regions_written, sections=html_sections,
                              results_by_record_id=results_by_record_id,
                              config=options, job_id=job_id, page_title=page_title,
                              records_without_regions=record_layers_without_regions,
                              svg_tooltip=svg_tooltip,
                              annotation_filename=annotation_filename)
        result_file.write(aux)


def update_cds_description(js_record, annotations):
    id_to_annotations = {}
    template = html_renderer.FileTemplate(path.get_full_path(__file__, "templates", "cds_detail.html"))
    for annotation in annotations:
        if annotation["id"]:
            id_to_annotations[annotation["id"]] = annotation
        if annotation["name"]:
            id_to_annotations[annotation["name"]] = annotation
    for reg_idx, js_region in enumerate(js_record["regions"]):
        for cds_idx, js_cds in enumerate(js_region["orfs"]):
            desc_html = js_cds["description"]
            re_locus = re.search("Locus tag: (.+)<br>", desc_html)
            locus_tag = re_locus.group(1)
            if locus_tag in id_to_annotations:
                annotation = id_to_annotations[locus_tag]
            else:
                re_protid = re.search("Protein ID: (.+)<br>", desc_html)
                protein_id = re_protid.group(1)
                if protein_id in id_to_annotations:
                    annotation = id_to_annotations[protein_id]
                else:
                    re_gid = re.search("Gene: (.+)<br>", desc_html)
                    gene_id = re_gid.group(1)
                    if gene_id in id_to_annotations:
                        annotation = id_to_annotations[gene_id]
                    else:
                        annotation = None
            additional_annots = "<br>added<br>"
            js_cds["description"] = re.sub("(\(total: (.+) nt\)<br>)", r"\1{}".format(template.render(annotation=annotation)), str(js_cds["description"]))