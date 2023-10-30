# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import os
import unittest
from unittest.mock import Mock, patch

from antismash.common import json, subprocessing, utils
from antismash.common.secmet.test.helpers import DummyCDS
from antismash.common.test.helpers import FakeHSPHit, FakeHit
from antismash.detection.genefunctions.tools import core
from antismash.detection.genefunctions.tools.halogenases.data_structures import (
    Conventionality,
    HalogenaseHit,
    Match,
    MotifDetails,
    Profile,
)
from antismash.detection.genefunctions.tools.halogenases.flavin_dependent import substrate_analysis
from antismash.detection.genefunctions.tools.halogenases.flavin_dependent.substrate_analysis import (
    categorise_on_substrate_level,
    classify,
    extract_residues,
    run_halogenase_phmms,
    set_conventionality,
)


def create_dummy_profile(**kwargs):
    options = {
        "name": "name",
        "description": "desc",
        "profile_name": "prof_name",
        "filename": os.devnull,
        "cutoffs": (5, 10),
    }
    options.update(kwargs)
    return Profile(**options)


def create_dummy_hit(**kwargs):
    options = {
        "query_id": "query",
        "reference_id": "reference",
        "bitscore": 50.,
        "evalue": 1e-20,
        "conventionality_residues": {},
        "profile": create_dummy_profile(),
    }
    options.update(**kwargs)
    return HalogenaseHit(**options)


class BlockHmmer(unittest.TestCase):
    def setUp(self):
        self.hmmpfam = patch.object(subprocessing, "run_hmmpfam2",
                                    side_effect=ValueError("hmmpfam2 should not run"))
        self.hmmpfam.start()
        self.hmmscan = patch.object(subprocessing.hmmscan, "run_hmmscan",
                                    side_effect=ValueError("hmmscan should not run"))
        self.hmmscan.start()
        self.hmmsearch = patch.object(subprocessing.hmmsearch, "run_hmmsearch",
                                      side_effect=ValueError("hmmsearch should not run"))
        self.hmmsearch.start()

    def tearDown(self):
        self.hmmpfam.stop()
        self.hmmscan.stop()
        self.hmmsearch.stop()


class TestComponents(BlockHmmer):
    def test_categorisation_with_no_hits(self):
        with patch.object(substrate_analysis, "retrieve_fdh_signature_residues", return_value={}):
            result = categorise_on_substrate_level(DummyCDS(), [])
        assert not result

    def test_conversion(self):
        fdh = create_dummy_hit()
        fdh.potential_matches.extend([
            Match(
                "prof_A",
                "flavin",
                "flavin-dependent",
                confidence=0.5,
                consensus_residues="MAGIC",
                substrate="some_sub",
                target_positions=(3, 4),
                number_of_decorations="deco",
            ),
            Match(
                "prof_B",
                "flavin",
                "flavin-dependent",
                confidence=0.1,
                consensus_residues="BEANS",
                substrate="other_sub",
                target_positions=(1,),
            ),
        ])
        rebuilt = HalogenaseHit.from_json(json.loads(json.dumps(fdh.to_json())))
        assert isinstance(rebuilt, HalogenaseHit)
        assert json.loads(json.dumps(fdh.to_json())) == json.loads(json.dumps(rebuilt.to_json()))

    @patch.object(subprocessing.hmmsearch, "run_hmmsearch",
                  return_value=[FakeHit(1, 2, 1000, "foo")])
    def test_run_halogenase_phmms(self, _patched):
        signature = create_dummy_profile(cutoffs=(300,))
        for value in _patched.return_value:
            value.hsps = [FakeHSPHit("foo", "query", bitscore=250)]

        negative_test_halogenase_hmms_by_id = run_halogenase_phmms("", [signature])
        assert not negative_test_halogenase_hmms_by_id

        for value in _patched.return_value:
            value.hsps = [FakeHSPHit("foo", "query_name", bitscore=1000)]

        results = run_halogenase_phmms("", [signature])["query_name"]
        assert results
        for hit in results:
            assert isinstance(hit, HalogenaseHit)

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="PREFIXWIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_conventional_non_specific(self, _patched_extract):
        name = "CtoA"
        hit = HalogenaseHit(
            query_id=name, bitscore=200, evalue=1e-20, reference_id="all_conventional_FDH",
            profile=create_dummy_profile(name="conventional"),
        )
        assert not hit.conventionality_residues
        result = substrate_analysis.set_conventionality(DummyCDS(), [hit])
        assert result is hit
        assert hit.conventionality_residues == {"W.W.I.": "WIWVIR"}
        specific_matches = substrate_analysis.categorise_on_substrate_level(
            DummyCDS(locus_tag=name), [hit],
        )
        assert not hit.potential_matches
        assert not specific_matches


class TestGetBest(BlockHmmer):
    def setUp(self):
        super().setUp()
        self.high_confidence_match = Match(
            "prof_A",
            "flavin",
            "flavin-dependent",
            confidence=0.9,
            consensus_residues="MAGIC",
            substrate="some_sub",
            target_positions=(3, 4),
            number_of_decorations="deco",
        )
        self.low_confidence_match = Match(
            "prof_B",
            "flavin",
            "flavin-dependent",
            confidence=0.1,
            consensus_residues="BEANS",
            substrate="other_sub",
            target_positions=(1,),
        )
        self.matches = [self.high_confidence_match, self.low_confidence_match]

    def test_no_matches(self):
        assert not create_dummy_hit().get_best_matches()

    def test_winner(self):
        best = create_dummy_hit(
            potential_matches=self.matches,
            conventionality=Conventionality.CONVENTIONAL,
        ).get_best_matches()
        assert len(best) == 1
        assert best[0] is self.high_confidence_match

    def test_multiple_equal(self):
        matches = [self.high_confidence_match, self.high_confidence_match]
        fdh = create_dummy_hit(potential_matches=matches)
        assert fdh.potential_matches == matches

    def test_single(self):
        matches = [self.low_confidence_match]
        best = create_dummy_hit(potential_matches=matches).get_best_matches()
        assert best == matches


class TestExtract(BlockHmmer):
    def test_no_alignment(self):
        with patch.object(utils, "extract_from_alignment", side_affect=RuntimeError("should not have been called")):
            with patch.object(subprocessing.hmmpfam, "get_alignment_against_profile", return_value=None):
                assert extract_residues("MAGIC", [1, 5, 6], Mock()) is None

    def test_alignment_passed(self):
        positions = [1, 5, 6]
        alignment = Mock()
        with patch.object(utils, "extract_from_alignment", return_value="dummy") as patched_util:
            with patch.object(subprocessing.hmmpfam, "get_alignment_against_profile", return_value=alignment):
                result = extract_residues("MAGIC", positions, Mock())
                assert result == "dummy"
            patched_util.assert_called_once_with(alignment, positions)

    def test_invalid_positions(self):
        with self.assertRaisesRegex(ValueError, "without positions"):
            extract_residues("MAGIC", [], Mock())


class TestSpecificAnalysis(BlockHmmer):
    def specific_analysis_test(self, name, *, bitscore=200, profile_id="dummy_profile", profile=None):
        cds = DummyCDS(locus_tag=name)

        if profile is None:
            profile = create_dummy_profile(profile_name=profile_id)
        hit = HalogenaseHit(query_id=name, bitscore=bitscore, evalue=1e-10, reference_id=profile_id, profile=profile)
        scan_result = [Mock(id=name, hits=True)]
        search_result = [Mock(id=name, hsps=[Mock(hit_id=name, query_id=profile_id, bitscore=bitscore, evalue=1e-10)])]
        with patch.object(subprocessing.hmmsearch, "run_hmmsearch", return_value=search_result):
            with patch.object(subprocessing.hmmscan, "run_hmmscan", return_value=scan_result):
                with patch.object(substrate_analysis, "scan_profiles_for_functions", return_value={name: hit}):
                    return classify([cds], _options=None)

    @patch.object(substrate_analysis, "scan_profiles_for_functions", return_value={})
    def test_no_hits(self, _patched_hmmscan):
        results = classify([DummyCDS()], _options=None)
        assert not results.best_hits

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="VALAMIVALAMI")
    def test_unconventional(self, _patched_extract_residues):
        name = "VatD"
        # have both types hit and ensure that only the best comes back
        conventional = create_dummy_hit(query_id=name, reference_id="conventional_FDH", bitscore=150,
                                        profile=create_dummy_profile(name="conventional"))
        unconventional = create_dummy_hit(query_id=name, reference_id="unconventional_FDH", bitscore=200,
                                          profile=create_dummy_profile(name="unconventional"))
        with patch.object(substrate_analysis, "search_conserved_motif", return_value="mismatch"):
            hit = set_conventionality(DummyCDS(), [conventional, unconventional])
        assert hit.conventionality == Conventionality.UNCONVENTIONAL
        assert not hit.is_conventional()
        assert not hit.potential_matches

    @patch.object(substrate_analysis, "extract_residues",
                  return_value="WIWVIRYGMIGDAASVIDAYYSQGVSLALVT")
    def test_conventional(self, _patched_extract_residues):
        name = "CtoA"
        result = self.specific_analysis_test(name, bitscore=200, profile_id="all_general_FDH")
        assert len(result.best_hits) == 1
        assert result.best_hits[name].conventionality_residues == {"W.W.I.": "WIWVIR"}


class TestDataStructures(unittest.TestCase):
    def new_via_json(self, item):
        # ensure it's proper json by converting to a string and back
        new = type(item).from_json(json.loads(json.dumps(item)))
        assert json.loads(json.dumps(new)) == json.loads(json.dumps(item))
        return new

    def test_match_conversion(self):
        match = Match(
            profile="dummy_profile", cofactor="some_cofactor", family="a family",
            confidence=0.2, consensus_residues="MAGIC", substrate="a substrate",
            target_positions=(5, 7), number_of_decorations="two_or_three_maybe",
        )

        rebuilt = self.new_via_json(match)
        assert rebuilt == match

    def test_match_equality(self):
        values = {
            "profile": "dummy_profile",
            "cofactor": "some_cofactor",
            "family": "a family",
            "confidence": 0.2,
            "consensus_residues": "MAGIC",
        }

        match = Match(**values)
        assert Match(**values) == match

        # and with any modification, equality should fail
        for key, val in values.items():
            copy = dict(values)
            copy[key] = val * 2  # arbitrary, any difference will do
            assert Match(**copy) != match

    def test_motif_string_equality(self):
        motif = MotifDetails(name="motif name", positions=(1, 6), residues="MT",
                             substrate="something", decorations="dec")
        assert motif == "MT" == motif.residues

    def test_motif_class_equality(self):
        # only the residues are meaningful for comparison
        residues = "YES"
        first = MotifDetails(name="motif name", positions=(1, 2, 6), residues=residues,
                             substrate="something", decorations="dec")
        different = MotifDetails(name="motif name", positions=(1, 6), residues="NO",
                                 substrate="something", decorations="dec")
        assert first.residues != different.residues
        assert first != different
        same = MotifDetails(name="new name", positions=(5, 6, 12), residues=residues,
                            substrate="a", decorations="other")
        assert first.residues == same.residues
        assert first == same
