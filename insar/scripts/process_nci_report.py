#!/usr/bin/env python

import datetime
import os
from pathlib import Path
import pandas as pd
import structlog
import shutil
import json
from datetime import datetime, timedelta

# from insar.scripts.process_gamma import ARDWorkflow

# Note: While working with dates, we primarily keep things in the original DT_FMT_SHORT string
# formatting, we only convert to a date if we need to do date/time related ops.
# This just keeps the code simpler (no useless conversions to/from string)
DT_FMT_SHORT = "%Y%m%d"
DT_FMT_LONG = "%Y-%m-%d"


def query_out_dir(dir: Path):
    """
    Scans a gamma insar output directory for information about the job that was processed.
    """
    workflow_achieved = None

    # Load top level data
    all_scene_dates = []
    all_ifg_date_pairs = []

    with (dir / "lists" / "scenes.list").open("r") as scenes_file:
        all_scene_dates = [line for line in scenes_file]

    with (dir / "lists" / "ifgs.list").open("r") as ifgs_file:
        all_ifg_date_pairs = [line.split(",") for line in ifgs_file]

    with (dir / "metadata.json").open("r") as metadata_file:
        metadata = json.load(metadata_file)

    pols = metadata["polarizations"]

    # Sanity check
    assert len(all_scene_dates) == int(metadata["num_scene_dates"])

    # Iterate all SLC detect completed products
    completed_slcs = []
    completed_coreg_slcs = []
    completed_backscatter = []

    for date_str in all_scene_dates:
        date_dir = dir / "SLC" / date_str
        if not date_dir.exists():
            continue

        # Check to see if data outputs exist
        is_complete = True

        for pol in pols:
            is_complete = len(list(date_dir.glob(""))) > 0

        if is_complete:
            completed_coreg_slcs.append(date_str)

    # Iterate all IFGs to detect completed products
    completed_ifgs = []

    for master_date, slave_date in all_ifg_date_pairs:
        date_dir = dir / "INT" / f"{master_date}-{slave_date}"
        if not date_dir.exists():
            continue

        # Check to see if data outputs exist
        is_complete = False  # TODO

        if is_complete:
            completed_ifgs.append(date_dir.name.split("-"))

    # Determine missing products
    missing_slc_dates = set(all_scene_dates) - set(completed_coreg_slcs)
    missing_backscatter_dates = set(all_scene_dates) - set(completed_backscatter)
    missing_ifgs = set(all_ifg_date_pairs) - set(completed_ifgs)

    return {
        "all_scene_dates": all_scene_dates,
        "all_ifg_date_pairs": all_ifg_date_pairs,
        "completed_slcs": completed_coreg_slcs,
        "missing_slcs": missing_slc_dates,
        "completed_coregs": completed_coreg_slcs,
        "missing_coregs": missing_slc_dates,
        "completed_backscatter": completed_coreg_slcs,
        "missing_backscatter": missing_slc_dates,
        "completed_ifgs": completed_ifgs,
        "missing_ifgs": missing_ifgs,
        "metadata": metadata,
    }


def query_job_dir(dir: Path):
    """
    Scans a gamma insar job directory for information about the job that was processed.
    """
    status_log_path = dir / "status.log"
    insar_log_path = dir / "insar-log.jsonl"

    # Load metadata
    with (dir / "metadata.json").open("r") as metadata_file:
        metadata = json.load(metadata_file)

    # Parse all job stdouts
    job_runs = []

    for stdout_path in dir.glob("job*.bash.o*"):
        with stdout_path.open("r") as stdout_file:
            parsing_enabled = False

            for line in stdout_file:
                # Skip until we hit an NCI resource usage report
                if not parsing_enabled:
                    prefix = "resource usage on "
                    idx = line.lower().find(prefix)
                    if idx > -1:
                        idx += len(prefix)
                        timestamp = datetime.strptime(
                            line[idx:].strip(), "%Y-%m-%d %H:%M:%S:"
                        )
                        parsing_enabled = True
                        parsing_line = 1

                    continue

                # Example:
                #                 Resource Usage on 2021-05-11 13:42:00:
                #   Job Id:             22160355.gadi-pbs
                #   Project:            dg9
                #   Exit Status:        0
                #   Service Units:      77.79
                #   NCPUs Requested:    48                     NCPUs Used: 48
                #                                           CPU Time Used: 12:53:21
                #   Memory Requested:   192.0GB               Memory Used: 131.02GB
                #   Walltime requested: 02:00:00            Walltime Used: 00:48:37
                #   JobFS requested:    400.0GB                JobFS used: 0B
                if parsing_line == 1:
                    prefix, job_id = line.split(":")
                    assert prefix.strip().lower() == "job id"
                elif parsing_line == 2:
                    prefix, compute_project = line.split(":")
                    assert prefix.strip().lower() == "project"
                elif parsing_line == 3:
                    prefix, exit_status = line.split(":")
                    assert prefix.strip().lower() == "exit status"
                elif parsing_line == 4:
                    prefix, service_units = line.split(":")
                    assert prefix.strip().lower() == "service units"
                elif parsing_line == 8:
                    prefix = "walltime used:"
                    idx = line.lower().find(prefix)
                    assert idx > -1
                    idx += len(prefix)
                    h, m, s = [int(i) for i in line[idx:].split(":")]
                    walltime_used = timedelta(hours=h, minutes=m, seconds=s)

                parsing_line += 1

                if parsing_line > 9:
                    parsing_enabled = False

                    job_runs.append(
                        {
                            "timestamp": timestamp,
                            "job_id": job_id,
                            "compute_project": compute_project,
                            "exit_status": exit_status,
                            "service_units": float(service_units),
                            "walltime_used": walltime_used,
                        }
                    )

    # Determine total walltime and SU used
    total_walltime = timedelta()
    total_service_units = 0

    for run in job_runs:
        total_walltime += run["walltime_used"]
        total_service_units += run["service_units"]

    # Check Luigi structure to detect completed products
    tfs = metadata["track_frame_sensor"]
    num_completed_backscatter = len(list((dir / tfs).glob("*_coreg_logs.out")))
    num_completed_ifgs = len(list((dir / tfs).glob("*_ifg_*_status_logs.out")))

    # Parse status.log to detect run/fail/completed products
    started_backscatter = {}
    started_ifgs = {}
    failed_backscatter = {}
    failed_ifgs = {}
    completed_backscatter = {}
    completed_ifgs = {}

    # In the future we may need to distinguish between completed/failed products for different polarisations?
    with (dir / "status.log").open("r") as status_file:
        for line in status_file:
            entry = json.loads(line[line.index("{") :])
            timestamp = entry["timestamp"]

            if "Beginning SLC coregistration" in entry["event"]:
                scene_date = entry["slave_date"]
                entry = [
                    timestamp,
                    entry["polarization"],
                    entry["master_date"],
                    entry["slave_date"],
                ]

                if scene_date in completed_backscatter:
                    del completed_backscatter[scene_date]

                if scene_date in failed_backscatter:
                    del failed_backscatter[scene_date]

                started_backscatter[scene_date] = entry
            elif "SLC coregistration complete" in entry["event"]:
                scene_date = entry["slave_date"]
                completed_backscatter[scene_date] = [
                    timestamp,
                    entry["polarization"],
                    entry["master_date"],
                    entry["slave_date"],
                ]

                # TODO: Record time taken to compute?
            elif "SLC coregistration failed with exception" in entry["event"]:
                # TODO: This reall means coreg 'or' backscatter failed (unsure which until we separate them)
                scene_date = entry["slave_date"]
                error = entry.exception
                failed_backscatter[scene_date] = [
                    timestamp,
                    entry["polarization"],
                    entry["master_date"],
                    entry["slave_date"],
                    error,
                ]

            elif "Beginning interferogram processing" in entry["event"]:
                date_pair = [entry["master_date"], entry["slave_date"]]
                entry = [timestamp, entry["master_date"], entry["slave_date"]]

                # Remove failed/completed status if we resumed!
                if date_pair in completed_ifgs:
                    del completed_ifgs[date_pair]

                if date_pair in failed_ifgs:
                    del failed_ifgs[date_pair]

                started_ifgs[date_pair] = entry
            elif "Interferogram complete" in entry["event"]:
                date_pair = [entry["master_date"], entry["slave_date"]]
                entry = [timestamp, entry["master_date"], entry["slave_date"]]
                completed_ifgs[date_pair] = entry

                # TODO: Record time taken to compute?
            elif "Interferogram failed with exception" in entry["event"]:
                date_pair = [entry["master_date"], entry["slave_date"]]
                error = entry.exception

                failed_ifgs[date_pair] = [
                    timestamp,
                    entry["master_date"],
                    entry["slave_date"],
                    error,
                ]

    print("num_completed_backscatter:", num_completed_backscatter)
    print("completed_backscatter:", len(completed_backscatter))
    print("failed_backscatter:", len(failed_backscatter))
    print("num_completed_ifgs:", num_completed_ifgs)
    print("completed_ifgs:", len(completed_ifgs))
    print("failed_ifgs:", len(failed_ifgs))
    assert num_completed_backscatter == len(completed_backscatter) + len(
        failed_backscatter
    )
    assert num_completed_ifgs == len(completed_ifgs) + len(failed_ifgs)

    return {
        "job_runs": job_runs,
        "total_walltime": total_walltime,
        "total_service_units": total_service_units,
        "started_backscatter": started_backscatter,
        "started_ifgs": started_ifgs,
        "failed_backscatter": failed_backscatter,
        "failed_ifgs": failed_ifgs,
        "completed_backscatter": completed_backscatter,
        "completed_ifgs": completed_ifgs,
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate a simple summary report on the outputs/progress of a processing job."
    )

    parser.add_argument(
        "path",
        type=str,
        help="The path of the job to report on (or a base dir from which --job-id queries from)",
    )

    parser.add_argument(
        "--job-id",
        type=str,
        help="The name of the processing job to look for, relative to the path (if path is a base path)",
    )

    parser.add_argument(
        "--json",
        action="store_true",
        help="Changes the output from plain human readable text into JSON",
    )

    args = parser.parse_args()

    # Find/query directories
    target_dir = Path(args.path)
    if not target_dir.exists():
        print("Provided path does not exist:", target_dir)
        exit(1)

    job_dir = None
    out_dir = None

    if args.job_id:
        job_dirs = list(target_dir.glob(f"job*{args.job_id}"))
        assert len(job_dirs) == 0 or len(job_dirs) == 1

        if job_dirs:
            job_dir = job_dirs[0]

        if (target_dir / args.job_id).exists():
            out_dir = target_dir / args.job_id

    else:
        is_job_dir = (target_dir / "insar-log.jsonl").exists()
        if is_job_dir:
            job_dir = target_dir
        else:
            out_dir = target_dir

    if not job_dir and not out_dir:
        print("Failed to find any job logs or data in target directory")
        exit(1)

    # Load metadata (should exist in both work and job dir)
    with (target_dir / "metadata.json").open("r") as metadata_file:
        metadata = json.load(metadata_file)

    pols = metadata["polarizations"]

    # Try and fill any gaps if the original dirs are still valid
    if not job_dir and Path(metadata["original_job_dir"]).exists():
        job_dir = Path(metadata["original_job_dir"])

    if not out_dir and Path(metadata["original_work_dir"]).exists():
        job_dir = Path(metadata["original_work_dir"])

    # Query data from job/output directories
    job_query = query_job_dir(job_dir) if job_dir else None
    out_query = query_out_dir(out_dir) if out_dir else None

    # Sanity check job dir <-> work dir correlation and use both to get complete reporting coverage
    if job_query and out_query:
        num_total_slc_scenes = len(out_query["all_scene_dates"])
        num_total_ifg_scenes = len(out_query["all_ifg_date_pairs"])

        num_completed_backscatter = len(job_query["completed_backscatter"])
        num_completed_ifgs = len(job_query["completed_ifgs"])
        num_completed_backscatter_outdir = len(out_query["completed_backscatter"])
        num_completed_ifgs_outdir = len(out_query["completed_ifgs"])
        assert num_completed_backscatter == num_completed_backscatter_outdir
        assert num_completed_ifgs == num_completed_ifgs_outdir

        num_failed_backscatter = len(job_query["failed_backscatter"])
        num_failed_ifgs = len(job_query["failed_ifgs"])

        num_missing_backscatter = len(out_query["missing_ifgs"])
        num_missing_ifgs = len(out_query["missing_ifgs"])

        # TODO: Remove failed from missing!

        assert (
            num_completed_backscatter + num_failed_backscatter + num_missing_backscatter
            == num_total_slc_scenes
        )
        assert (
            num_completed_ifgs + num_failed_ifgs + num_missing_ifgs
            == num_total_ifg_scenes
        )

    # If we only have a job dir & data doesn't exist, we can't get an accurate list of expected outputs from
    # the files that alredy exist (Note: we 'could' run the DB query again...)
    elif job_query:
        # Crude estimate if we don't have an out dir to read scene lists from...
        num_total_slc_scenes = len(job_query["started_backscatter"])
        num_total_ifg_scenes = len(job_query["started_ifgs"])

        num_completed_backscatter = len(job_query["completed_backscatter"])
        num_completed_ifgs = len(job_query["completed_ifgs"])

        num_failed_backscatter = len(job_query["failed_backscatter"])
        num_failed_ifgs = len(job_query["failed_ifgs"])

        num_missing_backscatter = "?"
        num_missing_ifgs = "?"

    # If we only have a data dir, we don't have any failure info (as all that is logged in the job dir)
    else:
        num_total_slc_scenes = len(out_query["all_scene_dates"])
        num_total_ifg_scenes = len(out_query["all_ifg_date_pairs"])

        num_completed_backscatter = len(out_query["completed_backscatter"])
        num_completed_ifgs = len(out_query["completed_ifgs"])

        num_failed_backscatter = "?"
        num_failed_ifgs = "?"

        num_missing_backscatter = len(out_query["missing_ifgs"])
        num_missing_ifgs = len(out_query["missing_ifgs"])

    # For each failed scene, report errors
    if job_query:
        for date_pair, entry in job_query["failed_backscatter"]:
            print(date_pair, "backscatter failed with:")
            print(entry["error"])
            print("")

            # TODO: Extract gamma stdout if it was a gamma error?

        for date_pair, entry in job_query["failed_ifgs"]:
            print(date_pair, "IFG failed with:")
            print(entry["error"])
            print("")

            # TODO: Extract gamma stdout if it was a gamma error?

    # Report NCI processing details (resource usage / node info / total job success or fail / etc)
    if job_query:
        for run in job_query["job_runs"]:
            started = run["timestamp"]
            walltime = run["walltime_used"]
            su_used = run["service_units"]
            proj = run["compute_project"]
            print(f"[{started}] Job ran for {walltime} using {su_used} SU from {proj}")

        print("")
        print("Total walltime:", job_query["total_walltime"])
        print("Total SU used:", job_query["total_service_units"])

    # Report completion metrics (including estimated walltime/KSU to complete)
    if num_total_slc_scenes and not isinstance(num_total_slc_scenes, str):
        completed_backscatter_pct = round(
            (num_completed_backscatter / num_total_slc_scenes) * 100.0
        )
    else:
        completed_backscatter_pct = num_total_slc_scenes  # copy the str value

    if num_total_ifg_scenes and not isinstance(num_total_ifg_scenes, str):
        completed_ifg_pct = round((num_completed_ifgs / num_total_ifg_scenes) * 100.0)
    else:
        completed_ifg_pct = num_total_ifg_scenes  # copy the str value

    print(
        f"Completed {num_completed_backscatter}/{num_total_slc_scenes} ({completed_backscatter_pct}%) backscatter products!"
    )
    print(
        f"Completed {num_completed_ifgs}/{num_total_ifg_scenes} ({completed_ifg_pct}%) IFG products!"
    )
    print(f"Missing {num_missing_backscatter} backscatter products!")
    print(f"Missing {num_missing_ifgs} IFG products!")
    print(f"Failed {num_failed_backscatter} backscatter products!")
    print(f"Failed {num_failed_ifgs} IFG products!")
