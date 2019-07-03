# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" MiBIG specific sideloading """

import json
from os import path
import sys
from typing import Any, Dict, List, Optional
import logging

from Bio import Entrez
from xml.dom import minidom
import time

from antismash.common.module_results import DetectionResults
from antismash.common.secmet import CDSFeature, SubRegion, Record
from antismash.common.secmet.locations import FeatureLocation, CompoundLocation


class MibigAnnotations(DetectionResults):
    def __init__(self, record_id: str, area: SubRegion, data: Dict[str, Any], cache_file: str) -> None:
        super().__init__(record_id)
        self.data = data # holds the original annotation json data
        # save calculated loci (relative to record), not annotated ones
        loci = data["cluster"]["loci"]
        assert (loci["accession"] == record_id) or (loci["accession"].split("MIBIG:")[-1] == record_id)
        self.area = area
        # fetch/update cached information (for taxonomy, etc.)
        cached = load_cached_information(data, cache_file)
        # save extra information from cache
        self.taxonomy = cached["taxonomy"][data["cluster"]["ncbi_tax_id"]]
        
    def get_predicted_subregions(self) -> List[SubRegion]:
        return [self.area]

    def to_json(self) -> Dict[str, Any]:
        # save only information critical for deciding reusability
        loci = self.data["cluster"]["loci"]
        return {
            "record_id": loci["accession"],
            "coords": (loci.get("start_coord", -1), loci.get("end_coord", -1)),
            "gene_annotations": self.data["cluster"].get("genes", {}).get("annotations", []),
            "extra_genes": self.data["cluster"].get("genes", {}).get("extra_genes", [])
        }

    @staticmethod
    def from_json(prev: Dict[str, Any], record: Record, annotations_file: str, cache_file: str) -> Optional["MibigAnnotations"]:
        with open(annotations_file) as handle:
            data = json.load(handle)
        
        # compare old vs new annotation, decide if we can 'reuse'
        can_reuse = True
        loci = data["cluster"]["loci"]
        gene_annotations = data["cluster"].get("genes", {}).get("annotations", [])
        extra_genes = data["cluster"].get("genes", {}).get("extra_genes", [])
        if loci["accession"] != prev["record_id"]:
            logging.debug("Previous result's record_id is not the same as the new one")
            can_reuse = False
        elif loci.get("start_coord", -1) != prev["coords"][0] or loci.get("end_coord", -1) != prev["coords"][1]:
            logging.debug("Previous result's start/end coordinate is not the same as the new one")
            can_reuse = False
        elif len(gene_annotations) != len(prev["gene_annotations"]):
            # lame implementation, fix when have time
            logging.debug("Updated gene annotations")
            can_reuse = False
        elif len(extra_genes) != len(prev["extra_genes"]):
            # lame implementation, fix when have time
            logging.debug("Updated gene annotations")
            can_reuse = False

        # if we can't reuse, stop running antismash, because CDS annotations won't be correct
        if can_reuse:
            product = ", ".join(data["cluster"]["biosyn_class"])
            loci_region = FeatureLocation(
                loci.get("start_coord", 1) - 1,
                loci.get("end_coord", len(record.seq))
            )
            area = SubRegion(loci_region, tool="mibig", label=product)
            return MibigAnnotations(record.id, area, data, cache_file)
        else:
            logging.error("Can't reuse MIBiG annotation, please turn off --reuse-results. Exiting..")
            sys.exit(1)


def mibig_loader(annotations_file: str, cache_file: str, record: Record) -> MibigAnnotations:
    """This method will be called only when not reusing data"""
    with open(annotations_file) as handle:
        data = json.load(handle)

    product = ", ".join(data["cluster"]["biosyn_class"])
    loci = data["cluster"]["loci"]
    loci_region = FeatureLocation(
        loci.get("start_coord", 1) - 1,
        loci.get("end_coord", len(record.seq))
    )
    area = SubRegion(loci_region, tool="mibig", label=product)

    # re-annotate CDSes in the Record
    if "genes" in data["cluster"]:
        # extra genes
        for gene in data["cluster"]["genes"].get("extra_genes", []):
            if "id" in gene and "location" in gene:
                # todo: check if exist in gbk
                exons = [FeatureLocation(exon["start"] - 1, exon["end"], strand=gene["location"]["strand"]) for exon in gene["location"]["exons"]]
                location = CompoundLocation(exons) if len(exons) > 1 else exons[0]
                translation = gene.get("translation", record.get_aa_translation_from_location(location))
                cds_feature = CDSFeature(location=location, locus_tag=gene["id"], translation=translation)
                record.add_cds_feature(cds_feature)
        # re-annotation
        for cds_feature in record.get_cds_features_within_location(area.location):
            locus_tag = cds_feature.locus_tag
            protein_id = cds_feature.protein_id
            name = cds_feature.gene
            for annot in data["cluster"]["genes"].get("annotations", []):
                if locus_tag and annot["id"] == locus_tag:
                    pass
                elif protein_id and annot["id"] == protein_id:
                    pass
                elif name and annot.get("name", None) == name:
                    pass
                else:
                    continue
                if "product" in annot:
                    cds_feature.product = annot["product"]

    return MibigAnnotations(record.id, area, data, cache_file)


def load_cached_information(annotations, cache_json_path, update=True):
    """"""
    if len(cache_json_path) > 0 and (path.exists(cache_json_path)):
        with open(cache_json_path) as handle:
            cached = json.load(handle)
    else:
        cached = {}

    ncbi_email = "mibig@secondarymetabolites.org"

    assert isinstance(cached, dict)

    # fetch taxonomy information
    if "taxonomy" not in cached:
        cached["taxonomy"] = {}
    ncbi_tax_id = annotations["cluster"]["ncbi_tax_id"]
    if ncbi_tax_id not in cached["taxonomy"]:
        cached["taxonomy"][ncbi_tax_id] = get_ncbi_taxonomy(ncbi_tax_id, ncbi_email)
    
    # fetch BibTex for publications
    # ....

    if update:
        # update cache file
        save_cached_information(cached, cache_json_path)

    return cached


def save_cached_information(cached, cache_json_path):
    """"""
    with open(cache_json_path, "w") as handle:
        handle.write(json.dumps(cached, indent=4, separators=(',', ': '), sort_keys=True))

    
def get_ncbi_taxonomy(tax_id, email):
    """fetch taxonomy information from ncbi_tax_id"""
    taxonomy = []    
    Entrez.email = email
    num_try = 1
    while num_try < 6:
        try:
            logging.debug("Fetching taxonomy information from NCBI for tax_id:{}...".format(tax_id))
            dom = minidom.parse(Entrez.efetch(db="taxonomy", id=tax_id))
            for dom_taxon in dom.getElementsByTagName('Taxon'):
                taxid = dom_taxon.getElementsByTagName("TaxId")[0].firstChild.nodeValue
                name = dom_taxon.getElementsByTagName("ScientificName")[0].firstChild.nodeValue
                rank = dom_taxon.getElementsByTagName("Rank")[0].firstChild.nodeValue
                taxonomy.append({"name": name, "taxid": taxid, "rank": rank})
            break
        except:
            pass
        num_try += 1
        time.sleep(5)
    if len(taxonomy) > 1: # shuffle species to the end of the list
        taxonomy.append(taxonomy.pop(0))
        
    return taxonomy