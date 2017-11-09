# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
    Runs minowa and kr_streochem over PKS domains
"""


import logging

from .results import PKSResults
from .substrates_pks import count_pks_genes, run_minowa_predictor_pks_at, run_minowa_predictor_pks_cal, run_kr_stereochemistry_predictions, extract_pks_genes

def run_pks_substr_spec_predictions(genes) -> PKSResults:
    pks_genes = extract_pks_genes(genes)
    counted = count_pks_genes(genes)
    if counted != len(pks_genes):
        logging.critical("mismatching PKS genes counted and PKS names extracted")
    results = PKSResults()
    if pks_genes:
        signature_results, minowa_at_results = run_minowa_predictor_pks_at(pks_genes)
        results.method_results["signature"] = signature_results
        results.method_results["minowa_at"] = minowa_at_results
    if counted:
        minowa_cal_results = run_minowa_predictor_pks_cal(genes)
        kr_activity, kr_stereo = run_kr_stereochemistry_predictions(genes)
        results.method_results["minowa_cal"] = minowa_cal_results
        results.method_results["kr_activity"] = kr_activity
        results.method_results["kr_stereochem"] = kr_stereo
    return results