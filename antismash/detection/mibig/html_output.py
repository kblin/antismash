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
    all_htmls = []

    html = HTMLSections("mibig-general")
    html.add_detail_section("General", FileTemplate(path.get_full_path(__file__, "templates", "general.html")).render(data=results.data))
    all_htmls.append(html)

    html = HTMLSections("mibig-compounds")
    html.add_detail_section("Compounds", FileTemplate(path.get_full_path(__file__, "templates", "compounds.html")).render(compounds=results.data["cluster"]["compounds"]))
    all_htmls.append(html)

    html = HTMLSections("mibig-genes")
    genes = []
    annots = [annot for annot in data["cluster"].get("genes", {}).get("annotations", [])]
    for cds_feature in _record_layer.get_cds_features_within_location(region_layer.location):
        gene = {
            "locus_tag": cds_feature.locus_tag,
            "protein_id": cds_feature.protein_id,
            "gene": cds_feature.gene,
            "start": cds_feature.location.start + 1,
            "end": cds_feature.location.end,
            "strand": cds_feature.location.strand,
            "product": cds_feature.product,
            "aa_seq": cds_feature.translation,
            "nt_seq": cds_feature.extract(_record_layer.seq)
        }
        gene["functions"] = []
        for function in cds_feature.gene_functions:
            function_text = str(function.function)
            if function.tool == "rule-based-clusters":
                function_text += " ({})".format(function.description)
            elif function.tool == "smcogs":
                function_text += " ({})".format(function.description.split(" (")[0])
            gene["functions"].append(function_text)
        annot_idx = -1
        for i, annot in enumerate(annots):
            if annot["id"] == cds_feature.locus_tag or annot["id"] == cds_feature.protein_id or annot["name"] == cds_feature.gene:
                annot_idx = i
                break
        if (annot_idx >= 0):
            annot = annots.pop(annot_idx)
            for function in annot.get("functions", []):
                function_text = function["category"]
                if len(annot.get('tailoring', [])) > 0:
                    function_text += " ({}) ".format(", ".join(annot['tailoring']))
                for evidence in function['evidence']:
                    function_text += "<span class='mibig-gf-evidence-{}' title='{}'>{}</span>".format(evidence[0], evidence, evidence[0])
                gene["functions"].append(function_text)
            if "mut_pheno" in gene and len(gene.get("mut_pheno", "")) > 0:
                function_text = "Mutation phenotype: {}".format(gene["mut_pheno"])
                gene["functions"].append(function_text)
            if "product" in annot:
                gene["product"] = annot["product"]
        genes.append(gene)
    for annot in annots:
        gene = {
            "locus_tag": annot.get("id", "None"),
            "protein_id": "None",
            "gene": annot.get("name", "None"),
            "product": annot.get("product", "")
        }
        gene["functions"] = []
        for function in annot.get("functions", []):
            function_text = function["category"]
            if len(annot.get('tailoring', [])) > 0:
                function_text += " ({}) ".format(", ".join(annot['tailoring']))
            for evidence in function['evidence']:
                function_text += "<span class='mibig-gf-evidence-{}' title='{}'>{}</span>".format(evidence[0], evidence, evidence[0])
            gene["functions"].append(function_text)
        genes.append(gene)
    html.add_detail_section("Genes", FileTemplate(path.get_full_path(__file__, "templates", "genes.html")).render(genes=genes))
    all_htmls.append(html)

    if "polyketide" in data["cluster"]:
        html = HTMLSections("mibig-polyketide")
        html.add_detail_section("Polyketide", FileTemplate(path.get_full_path(__file__, "polyketide", "details.html")).render(pk=results.data["cluster"]["polyketide"]))
        all_htmls.append(html)

    if "nrp" in data["cluster"]:
        html = HTMLSections("mibig-nrp")
        html.add_detail_section("NRP", FileTemplate(path.get_full_path(__file__, "nrp", "details.html")).render(nrp=results.data["cluster"]["nrp"]))
        all_htmls.append(html)

    if "ripp" in data["cluster"]:
        html = HTMLSections("mibig-ripp")
        html.add_detail_section("RiPP", FileTemplate(path.get_full_path(__file__, "ripp", "details.html")).render(ripp=results.data["cluster"]["ripp"]))
        all_htmls.append(html)

    if "saccharide" in data["cluster"]:
        html = HTMLSections("mibig-saccharide")
        html.add_detail_section("Saccharide", FileTemplate(path.get_full_path(__file__, "saccharide", "details.html")).render(sac=results.data["cluster"]["saccharide"]))
        all_htmls.append(html)

    if "terpene" in data["cluster"]:
        html = HTMLSections("mibig-terpene")
        html.add_detail_section("Terpene", FileTemplate(path.get_full_path(__file__, "terpene", "details.html")).render(trp=results.data["cluster"]["terpene"]))
        all_htmls.append(html)

    html = HTMLSections("mibig-logs")
    html.add_detail_section("History", FileTemplate(path.get_full_path(__file__, "templates", "logs.html")).render(logs=sorted(results.data["changelog"], key=lambda log: log["version"])))
    all_htmls.append(html)

    return all_htmls
