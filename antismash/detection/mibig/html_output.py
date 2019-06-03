# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for the MIBiG sideloader """

from typing import List

from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.html_renderer import FileTemplate, HTMLSections, docs_link
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer


def will_handle(products: List[str]) -> bool:
    """ Returns true if one or more relevant products are present """
    return True


def generate_html(region_layer: RegionLayer, results: ModuleResults,
                  _record_layer: RecordLayer, _options_layer: OptionsLayer) -> List[HTMLSections]:

    data = results.data

    html_compound = HTMLSections("mibig-compounds")
    html_compound.add_detail_section("Compounds", FileTemplate(path.get_full_path(__file__, "templates", "compounds.html")).render(compounds=results.data["cluster"]["compounds"]))

    html_general = HTMLSections("mibig-general")
    html_general.add_detail_section("General", FileTemplate(path.get_full_path(__file__, "templates", "general.html")).render(data=results.data))

    return [html_compound, html_general]
