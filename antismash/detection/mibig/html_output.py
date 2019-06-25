# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for the MIBiG sideloader """

from typing import List

from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.html_renderer import FileTemplate, HTMLSections
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from .mibig import MibigAnnotations


def will_handle(_products: List[str]) -> bool:
    """ Returns true if one or more relevant products are present """
    return True


def generate_html(region_layer: RegionLayer, results: ModuleResults,
                  record_layer: RecordLayer, _options_layer: OptionsLayer) -> HTMLSections:
    assert isinstance(results, MibigAnnotations)
    data = results.data

    html = HTMLSections("mibig-general")
    taxonomy_text = " > ".join(["{} ({})".format(taxobj["name"], taxobj["rank"]) for taxobj in data["cluster"]["taxonomy"]])
    publications_links = []
    for pub in data["cluster"]["publications"]:
        [category, index] = pub.split(":")
        if category == "pubmed":
            url = "https://www.ncbi.nlm.nih.gov/pubmed/{}".format(index)
        elif category == "patent":
            url = "https://patents.google.com/patent/{}".format(index)
        elif category == "doi":
            url = "https://dx.doi.org/{}".format(index)
        elif category == "url":
            url = index
        publications_links.append(url)
    html.add_detail_section("General", FileTemplate(path.get_full_path(__file__, "templates", "general.html")).render(data=results.data, taxonomy_text=taxonomy_text, publications_links=publications_links))

    for compound in results.data["cluster"]["compounds"]:
        compound["keys"] = [key for key in compound.keys() if (key not in ["compound", "chem_struct"]) and (not isinstance(compound[key], list) or len(compound[key]) > 0)]
    html.add_detail_section("Compounds", FileTemplate(path.get_full_path(__file__, "templates", "compounds.html")).render(compounds=results.data["cluster"]["compounds"]),
                            class_name="mibig-compounds")

    genes = []
    annots = [annot for annot in data["cluster"].get("genes", {}).get("annotations", [])]
    for cds_feature in record_layer.get_cds_features_within_location(region_layer.location):
        gene = {
            "locus_tag": cds_feature.locus_tag,
            "protein_id": cds_feature.protein_id,
            "gene": cds_feature.gene,
            "start": cds_feature.location.start + 1,
            "end": cds_feature.location.end,
            "strand": cds_feature.location.strand,
            "product": cds_feature.product,
            "aa_seq": cds_feature.translation,
            "nt_seq": cds_feature.extract(record_layer.seq)
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
        if annot_idx >= 0:
            annot = annots.pop(annot_idx)
            for function in annot.get("functions", []):
                function_text = function["category"]
                if "tailoring" in annot:
                    function_text += " ({}) ".format(", ".join(annot['tailoring']))
                else:
                    function_text += " "
                for evidence in function['evidence']:
                    function_text += "<span class='mibig-gf-evidence-{}' title='{}'>{}</span>".format(evidence[0], evidence, evidence[0])
                gene["functions"].append(function_text)
            if "mut_pheno" in gene:
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
            if "tailoring" in annot:
                function_text += " ({}) ".format(", ".join(annot['tailoring']))
            else:
                function_text += " "
            for evidence in function['evidence']:
                function_text += "<span class='mibig-gf-evidence-{}' title='{}'>{}</span>".format(evidence[0], evidence, evidence[0])
            gene["functions"].append(function_text)
        genes.append(gene)
    html.add_detail_section("Genes", FileTemplate(path.get_full_path(__file__, "templates", "genes.html")).render(genes=genes),
                            class_name="mibig-genes")

    if "polyketide" in data["cluster"]:
        html.add_detail_section("Polyketide", FileTemplate(path.get_full_path(__file__, "polyketide", "details.html")).render(pk=results.data["cluster"]["polyketide"]),
                                class_name="mibig-polyketide")

    if "nrp" in data["cluster"]:
        html.add_detail_section("NRP", FileTemplate(path.get_full_path(__file__, "nrp", "details.html")).render(nrp=results.data["cluster"]["nrp"]),
                                class_name="mibig-nrp")

    if "ripp" in data["cluster"]:
        html.add_detail_section("RiPP", FileTemplate(path.get_full_path(__file__, "ripp", "details.html")).render(ripp=results.data["cluster"]["ripp"]),
                                class_name="mibig-ripp")

    if "saccharide" in data["cluster"]:
        html.add_detail_section("Saccharide", FileTemplate(path.get_full_path(__file__, "saccharide", "details.html")).render(sac=results.data["cluster"]["saccharide"]),
                                class_name="mibig-saccharide")

    if "terpene" in data["cluster"]:
        html.add_detail_section("Terpene", FileTemplate(path.get_full_path(__file__, "terpene", "details.html")).render(trp=results.data["cluster"]["terpene"]),
                                class_name="mibig-terpene")

    html.add_detail_section("History", FileTemplate(path.get_full_path(__file__, "templates", "logs.html")).render(logs=sorted(results.data["changelog"], key=lambda log: log["version"])),
                                class_name="mibig-logs")

    return html
