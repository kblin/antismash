# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
    Categorisation of halogenases based on family and function.
    It utilizes pHMMs with thresholds and the presence or absence of
    characteristic/signature residues to identify the substrate or
    the number of halogenations.

    Signature residues for the pyrrole substrates were determined by
    using information from an experimental study (https://doi.org/10.1073/pnas.1519695113).

    Signature residues for all the other substrates were extracted manually,
    by looking at the positions that have the same value.
"""


from collections import defaultdict
import re
from typing import Iterable, Optional

from antismash.common import (
    fasta,
    json,
    subprocessing,
    utils,
)
from antismash.common.secmet import (
    CDSFeature,
)
from antismash.common.secmet.qualifiers import GeneFunction
from antismash.config import ConfigType


from ...core import scan_profiles_for_functions
from ..data_structures import (
    Conventionality,
    Group,
    HalogenaseHit,
    HalogenaseResults,
    Match,
    MotifDetails,
    Profile,
    get_file_path,
)

GENERAL_FDH_PROFILES = [
    Profile(name="conventional", profile_name="all_conventional_FDH",
            description="Flavin-dependent halogenase",
            cofactor="flavin", family="flavin-dependent",
            cutoffs=(100,), filename=get_file_path("all_conventional_FDH.hmm")),
    Profile(name="unconventional", profile_name="unconventional_FDH",
            description="Unconventional flavin-dependent halogenase",
            cofactor="flavin", family="flavin-dependent",
            cutoffs=(100,), filename=get_file_path("unconventional_FDH.hmm")),
]

ALL_FDH_PROFILES = get_file_path("FDH.hmm")

# all conventional
GENERAL_FDH_MOTIFS = {
    "W.W.I.": list(range(206, 212)),
    "F.*P.*S.G": list(range(280, 305)),
}


def _load_subgroups(metadata_file: str) -> dict[str, Group]:
    with open(metadata_file, encoding="utf-8") as handle:
        raw = json.load(handle)
    return {key: Group.from_json(val) for key, val in raw.items()}


SUBGROUPS: dict[str, Group] = _load_subgroups(get_file_path("metadata.json"))


def _get_substrate_specific_profiles() -> list[Profile]:
    """ Collects unique substrate-specific pHMM profiles from the substrate-specific submodules"""
    profiles = []
    for submodule in SUBGROUPS.values():
        for profile in submodule.profiles.values():
            profiles.append(profile)
    return profiles


def retrieve_fdh_signature_residues(translation: str, hmm_result: HalogenaseHit,
                                    motifs: Iterable[MotifDetails],
                                    ) -> dict[str, str]:
    """ Extracts residues for each of the given motifs from an HMM hit

        Arguments:
            sequence: protein sequence
            hmm_result: instance of HmmResult class,
                        which contains information about the hit in a pHMM
            motifs: the motifs for which to extract signatures

        Returns:
            signature residues which were retrieved from a certain pHMM
    """
    signature_residues: dict[str, str] = {}
    for motif in motifs:
        if not motif.positions:  # then it's always present, just 'empty'
            signature_residues[motif.name] = ""
            continue
        residues = extract_residues(translation, motif.positions, hmm_result)
        if residues:
            signature_residues[motif.name] = residues
    return signature_residues


def extract_residues(sequence: str, positions: Iterable[int],
                     hmm_result: HalogenaseHit,
                     max_evalue: float = 0.1) -> Optional[str]:
    """ Get the signature residues from the pHMM for the searched protein sequence

        Arguments:
            sequence: protein sequence
            positions: list of position numbers in the pHMM
            hmm_result: properties of the halogenase pHMM hit
            max_evalue: maximum e-value

        Returns:
            residues that are present in the given positions
    """
    if not positions:
        raise ValueError("cannot extract residues without positions")
    assert hmm_result.profile
    alignment = subprocessing.hmmpfam.get_alignment_against_profile(sequence, hmm_result.profile.filename,
                                                                    hmm_result.reference_id, max_evalue=max_evalue)
    if not alignment:
        return None
    return utils.extract_from_alignment(alignment, positions)


def search_conserved_motif(cds: CDSFeature, motif_positions: list[int],
                           hmm_result: HalogenaseHit, motif_pattern: str,
                           ) -> str:
    """ Finds motifs of the given pattern within the CDS translation.

        Arguments:
            cds: gene/CDS and its properties
            motif_positions: positions of the conserved motifs in the pHMM
            hmm_result: details of the hit (e.g. bitscore, name of the profile, etc.)
            motif_pattern: pattern of the motif in regex

        Returns:
            the conserved motifs, if present, otherwise an empty string
    """
    signature_residues = extract_residues(cds.translation, motif_positions, hmm_result)

    if not signature_residues:
        return ""

    motif = re.search(motif_pattern, signature_residues)
    if not motif:
        return ""

    return motif[0]


def run_halogenase_phmms(cluster_fasta: str, profiles: list[Profile],
                         ) -> dict[str, list[HalogenaseHit]]:
    """ Check if protein sequences hit any pHMM

        Arguments:
            cluster_fasta: string of protein sequences in a fasta format
            profiles: a list of HMM profiles to run

        Returns:
            a dictionary mapping query names to a list of hits for that query, one per profile
    """
    halogenase_hmms_by_id: dict[str, list[HalogenaseHit]] = defaultdict(list)
    for profile in profiles:
        run_results = subprocessing.hmmsearch.run_hmmsearch(profile.filename, cluster_fasta)
        for run_result in run_results:
            for hsp in run_result.hsps:
                if hsp.bitscore > profile.cutoffs[-1]:
                    cds_name = hsp.hit_id
                    hit = HalogenaseHit(
                        query_id=cds_name,
                        reference_id=hsp.query_id,
                        bitscore=hsp.bitscore,
                        evalue=hsp.evalue,
                        profile=profile,
                    )
                    halogenase_hmms_by_id[cds_name].append(hit)
    return halogenase_hmms_by_id


def categorise_on_substrate_level(cds: CDSFeature, hmm_results: list[HalogenaseHit],
                                  ) -> list[Match]:
    """ Finds matches for each hit for a specific CDS within specific subgroups of
        halogenases.

        Arguments:
            cds: the CDS feature to categorise
            hmm_results: a list of HMM hits for the CDS

        Returns:
            a list of matches
    """
    if not hmm_results:
        return []

    matches = []

    for hit in hmm_results:
        for subgroup in SUBGROUPS.values():
            for profile in subgroup.get_matching_profiles(hit):
                residues = retrieve_fdh_signature_residues(cds.translation, hit, profile.motifs)
                matches.extend(profile.get_matches_from_hit(residues, hit))
                if matches:
                    # substrates are always within the conventional set at time of writing
                    hit.conventionality = Conventionality.CONVENTIONAL
                    break
    return matches


def set_conventionality(cds: CDSFeature, hits: list[HalogenaseHit]) -> Optional[HalogenaseHit]:
    """ Sets the conventionality of each given hit

        Arguments:
            cds: the CDS feature the hits belong to
            hits: a list of hits for which to set conventionality

        Returns:
            the singular best conventionality hit, if any
    """
    if not hits:
        return None
    # handle the case where multiple hits exist, whether in fragments or multiple domains
    best_hit = hits[0]
    for hit in hits[1:]:
        if hit.bitscore > best_hit.bitscore:
            best_hit = hit

    assert best_hit.profile

    for motif_pattern, positions in GENERAL_FDH_MOTIFS.items():
        conserved_motif = search_conserved_motif(cds, positions, best_hit, motif_pattern)
        if conserved_motif:
            best_hit.conventionality_residues[motif_pattern] = conserved_motif

    if best_hit.profile.name == "conventional" and not best_hit.conventionality_residues:
        best_hit.conventionality = Conventionality.AMBIGUOUS
        return best_hit

    best_hit.conventionality = Conventionality[best_hit.profile.name.upper()]
    return best_hit


def classify(cds_features: Iterable[CDSFeature], _options: ConfigType,
             ) -> HalogenaseResults:
    """ Categorisation of enzyme, categorises any halogenase in a cds in regions

        Arguments:
            cds_features: the CDS features to classify

        Returns:
            a results object containing all classifications and additional details
    """
    cds_by_name = {cds.get_name(): cds for cds in cds_features}

    hits = scan_profiles_for_functions(cds_features, ALL_FDH_PROFILES)
    if not hits:
        return HalogenaseResults()

    relevant = [cds_by_name[hit] for hit in hits]

    specific_fasta = fasta.get_fasta_from_features(relevant)
    general_hmm_hits = run_halogenase_phmms(specific_fasta, GENERAL_FDH_PROFILES)

    substrate_hmm_hits = {}
    if general_hmm_hits:
        specific_profiles = _get_substrate_specific_profiles()
        substrate_hmm_hits = run_halogenase_phmms(specific_fasta, specific_profiles)
    best_hits: dict[str, HalogenaseHit] = {}
    all_hits: dict[str, list[HalogenaseHit]] = defaultdict(list)

    for cds_name, conventionality_hits in general_hmm_hits.items():
        if not conventionality_hits:
            continue
        cds = cds_by_name[cds_name]

        cds_hits = all_hits[cds_name]
        general = set_conventionality(cds, conventionality_hits)
        if not general:
            continue

        cds_hits.append(general)
        best_hits[cds_name] = general

        if not general.is_conventional():
            continue

        substrate_hits = substrate_hmm_hits[cds_name]

        if not categorise_on_substrate_level(cds, substrate_hits):
            continue

        substrate_hits = sorted(substrate_hits, key=lambda x: x.confidence, reverse=True)
        cds_hits.extend(substrate_hits)
        # don't set a substrate-specific best hit if it's not the only one with the same confidence
        if len(substrate_hits) == 1 or substrate_hits[0].confidence > substrate_hits[1].confidence:
            best_hits[cds_name] = substrate_hits[0]

    functions = {name: GeneFunction.ADDITIONAL for name in best_hits}
    return HalogenaseResults(best_hits=best_hits, function_mapping=functions, all_hits=all_hits)
