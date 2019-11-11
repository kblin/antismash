#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Script to generate mibig html output from mibig json file."""

import json
from os import path
from sys import argv, exit
from subprocess import call
from shutil import rmtree
from datetime import datetime
from mibig.converters.read.top import Everything
from antismash.detection.mibig.mibig import load_cached_information

def write_log(text, file_path):
    with open(file_path, "a") as o:
        o.write("[{}] {}\n".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), text))

def _main():
    json_path = argv[1]
    gbk_folder = argv[2]
    cache_folder = argv[3]
    output_folder = argv[4]
    log_file_path = argv[5]
    if len(argv) > 6:
        use_source = argv[6] == "1"
    else:
        use_source = False
    antismash_path = path.abspath(path.dirname(__file__))

    with open(json_path, "r") as json_file:
        raw = json.load(json_file)
        data = Everything(raw)
        mibig_acc = data.cluster.mibig_accession
        gbk_acc = data.cluster.loci.accession
        gbk_path = path.join(gbk_folder, "{}.gbk".format(gbk_acc))
        cache_json_path = path.join(cache_folder, "{}.cache.json".format(mibig_acc))
        output_path = path.join(output_folder, mibig_acc)
        reusable_json_path = path.join(output_path, "{}.1.json".format(mibig_acc))

        # load/populate cache in advance, because it is needed to fetch taxonomy information
        cached = load_cached_information(data, cache_json_path, True)
        taxonomy = [tax_obj["name"] for tax_obj in cached["taxonomy"][data.cluster.ncbi_tax_id]]

        if "Bacteria" in taxonomy:
            taxon = "bacteria"
        elif "Fungi" in taxonomy:
            taxon = "fungi"
        elif "Viridiplantae":
            taxon = "plants"
        else:
            write_log("Unrecognizable taxons {} ({})".format(mibig_acc, ":".join(taxonomy)), log_file_path)
            exit(1)

        commands = [
            "python",
            path.join(antismash_path, "run_antismash.py"),
            "-v",
            "--cb-known",
            "--minlength",
            "1",
            "--taxon",
            taxon,
            "--mibig-mode",
            "--mibig-json",
            json_path,
            "--mibig-cache-json",
            cache_json_path,
            "--output-dir",
            output_path,
            gbk_path
        ]
        if not use_source:
            commands = ["antismash"] + commands[2:]
        print("Generating MIBiG output for {}".format(mibig_acc))
        reuse_success = False
        if path.exists(reusable_json_path):
            reuse_success = call(commands[:-1] + ["--reuse-results", reusable_json_path]) == 0
            if reuse_success:
                write_log("Successfully reused JSON file {}".format(reusable_json_path), log_file_path)
            else:
                write_log("Failed to reuse JSON file {}".format(reusable_json_path), log_file_path)
        if path.exists(output_path) and not reuse_success:
            # remove output path, proceed with caution!
            rmtree(output_path)
            write_log("Removed {}".format(output_path), log_file_path)
        if not reuse_success and call(commands) == 1:
            write_log("Failed to generate MIBiG page for {}".format(mibig_acc), log_file_path)
            exit(1)
        else:
            write_log("Successfully generated MIBiG page for {}".format(mibig_acc), log_file_path)
            print("Generating antiSMASH output for {}".format(mibig_acc))
            with open(path.join(output_path, "{}.1.json".format(mibig_acc)), "r") as result_json_txt:
                result_json = json.load(result_json_txt)
                assert len(result_json["records"]) == 1 and "antismash.detection.mibig" in result_json["records"][0]["modules"]
                region_acc = "{}.1.region001".format(mibig_acc)
            reusable_as5_json_path = path.join(output_path, "generated", "{}.json".format(region_acc))
            region_gbk_path = path.join(output_path, "{}.gbk".format(region_acc))
            if taxon in ["bacteria", "fungi"]:
                commands = [
                    "python",
                    path.join(antismash_path, "run_antismash.py"),
                    "-v",
                    "--minlength",
                    "1",
                    "--html-title",
                    "{} (generated)".format(mibig_acc),
                    "--cb-known",
                    "--taxon",
                    taxon,
                    "--output-dir",
                    path.join(output_path, "generated"),
                    region_gbk_path
                ]
                if not use_source:
                    commands = ["antismash"] + commands[2:]
                reuse_as5_success = False
                if path.exists(reusable_as5_json_path):
                    reuse_as5_success = call(commands[:-1] + ["--reuse-results", reusable_as5_json_path]) == 0
                if reuse_as5_success or call(commands) == 0:
                    write_log("Successfully generated antiSMASH5 page for {}".format(mibig_acc), log_file_path)
                else:
                    write_log("Failed to generate antiSMASH5 page for {}".format(mibig_acc), log_file_path)
                    exit(1)
            elif taxon == "plants":
                "Plant BGC is temporarily not supported"
                write_log("Plant BGCs {}".format(mibig_acc), log_file_path)

if __name__ == "__main__":
    _main()
