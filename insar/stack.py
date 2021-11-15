import datetime
import itertools
import json
from pathlib import Path
import re
from click.core import Option
import pandas as pd
from pandas._libs import missing
import geopandas as gpd
from shapely.geometry import box

from typing import Optional, List, Tuple, Union

from insar.project import ProcConfig
from insar.constant import SCENE_DATE_FMT
from insar.sensors import identify_data_source, get_data_swath_info, S1_ID, RSAT2_ID, PALSAR_ID
from insar.generate_slc_inputs import query_slc_inputs, slc_inputs
from insar.logs import STATUS_LOGGER

from insar.workflow.luigi.utils import simplify_dates, one_day

def load_stack_config(stack_proc_path: Union[str, Path]) -> ProcConfig:
    """
    Loads a stack's .proc config file from a supplied stack path.

    :param stack_proc_path:
        A path either directly to a stack's config.proc file that resides in the stack's
        output directory, or alternatively this may be a path to either the stack's
        output or job directories (which will be used to find the stack's .proc file)
    :returns:
        The stack's loaded ProcConfig
    """

    stack_proc_path = Path(stack_proc_path)

    if not stack_proc_path.exists():
        raise ValueError("Specified path does not exist!")

    if stack_proc_path.is_dir():
        metadata = stack_proc_path / "metadata.json"
        if not metadata.exists():
            raise ValueError("Expected stack dir - metadata.json missing from specified directory!")

        with metadata.open("r") as file:
            metadata = json.load(file)

        og_config = Path(metadata["original_work_dir"]) / "config.proc"

        return load_stack_config(og_config)

    else:
        if stack_proc_path.suffix != ".proc":
            raise ValueError("Expected stack dir or .proc config file!")

        with stack_proc_path.open("r") as proc_file_obj:
            return ProcConfig.from_file(proc_file_obj)


def load_stack_scene_dates(proc_config: ProcConfig) -> List[List[str]]:
    list_dir = proc_config.output_path / proc_config.list_dir
    result = []

    append_idx = 0
    scene_file = Path(list_dir / "scenes.list")

    while scene_file.exists():
        result.append(scene_file.read_text().strip().splitlines())

        append_idx += 1
        scene_file = Path(list_dir / f"scenes{append_idx}.list")

    return result


def load_stack_scenes(proc_config: ProcConfig) -> List[Tuple[datetime.date, List[Path]]]:
    scene_dates = load_stack_scene_dates(proc_config)
    slc_dir = Path(proc_config.output_path) / proc_config.slc_dir
    result = []

    for date in itertools.chain(*scene_dates):
        scene_dir = slc_dir / date
        if not scene_dir.exists():
            continue

        metadata_files = list(scene_dir.glob("metadata_*.json"))
        metadata = json.loads(metadata_files[0].read_text())

        date_source_files = []

        for values in metadata.values():
            if "src_url" in values:
                print("src_url!", values["src_url"])
                date_source_files.append(Path(values["src_url"]))
            else:
                print("UH OH", date)

        result.append((date, date_source_files))

    return result


def resolve_stack_scene_additional_files(
    slc_inputs_df: pd.DataFrame,
    proc_config: ProcConfig,
    polarisations: List[str],
    include_source_files: List[Path],
    shape_file: Optional[Path] = None,
) -> pd.DataFrame:
    """
    This function inserts additional source files into a stack query
    by resolving their dates and determining burst availability before
    finally inserting them into the stack's pandas dataframe.
    """
    download_dir = Path(proc_config.output_path) / proc_config.raw_data_dir

    swath_info_by_date = {}

    # Gather swath information for all source files, grouping by date
    for data_path in include_source_files:
        _, _, scene_date = identify_data_source(data_path)

        if scene_date not in swath_info_by_date:
            swath_info_by_date[scene_date] = []

        for swath_data in get_data_swath_info(data_path, download_dir / scene_date):
            if swath_data["polarization"] not in polarisations:
                STATUS_LOGGER.info(
                    "Skipping source data which does not match stack polarisations",
                    source_file=data_path,
                    source_pol=swath_data["polarization"],
                    stack_pols=polarisations
                )
                continue

            swath_info_by_date[scene_date].append(swath_data)

    if shape_file is not None:
        shp_df = gpd.read_file(shape_file)

    # Add each date's swath data info to the stack dataframe
    for date, swath_data_list in swath_info_by_date.items():
        # If we have a shapefile, determine burst availability for each swath
        #
        # Note: all acquisitions swaths share the same burst availability in a date
        # - as it's expected if many acquisitions exist for a single date, they're
        # - concatenated together (as the stack is geometrically larger than a single
        # - satellite acquisition)
        if shape_file is not None:
            # Note: This only applies to S1, and only S1 will have swaths w/ these names
            for swath in ["IW1", "IW2", "IW3"]:
                shp_swath_df = shp_df[shp_df.swath == swath]
                contained_primary_bursts = []

                for swath_data in swath_data_list:
                    if swath_data["swath"] != swath:
                        continue

                    burst_numbers = swath_data["burst_number"]
                    burst_centroids = swath_data["burst_centroid"]
                    if not burst_centroids:
                        continue

                    contained_swath_bursts = []

                    # Remove the centroids from the data object (we don't want them recorded in dataframe)
                    del swath_data["burst_centroid"]

                    for idx, row in shp_swath_df.iterrows():
                        for burst_num, centroid in zip(burst_numbers, burst_centroids):
                            if row.geometry.contains(centroid):
                                contained_primary_bursts.append(row.burst_num)
                                contained_swath_bursts.append(burst_num)

                    # Rewrite "burst_number" to only include those included by the shapefile
                    swath_data["burst_number"] = contained_swath_bursts

                # TODO: de-duplicate shapefile burst availability logic w/ the same code
                # in generate_slc_inputs.py
                missing_bursts = set(shp_swath_df.burst_num.values) - set(contained_primary_bursts)
                missing_bursts = [int(i) for i in missing_bursts]

                for swath_data in swath_data_list:
                    if swath_data["swath"] != swath:
                        continue

                    swath_data["missing_primary_bursts"] = missing_bursts

        slc_inputs_df = slc_inputs_df.append(swath_data_list, ignore_index=True)

    return slc_inputs_df


def resolve_stack_scene_query(
    proc_config: ProcConfig,
    include_queries: List[Union[Path, datetime.date, str]],
    sensors: List[str],
    orbit: str,
    polarisations: List[str],
    sensor_filters: List[Optional[str]],
    shape_file: Optional[Path] = None,
) -> List[Tuple[datetime.date, List[Path]]]:
    """
    TODO: Documentation
    """

    include_dates = []
    include_source_files = []

    for query in include_queries:
        # Strings may be a YYYYMMDD date, or a file path
        if isinstance(query, str):
            if re.match("\d{8}", query):
                query = datetime.datetime.strptime(query, SCENE_DATE_FMT).date()
                include_dates.append((query,query))

            else:
                include_source_files.append(Path(query))

        elif isinstance(query, tuple):
            assert(len(query) == 2)
            assert(all(isinstance(i, datetime.date) for i in query))

            include_dates.append(query)

        elif isinstance(query, datetime.date):
            include_dates.append((query,query))

        elif isinstance(query, Path):
            include_source_files.append(query)

    resolved_source_files = []
    slc_inputs_df = pd.DataFrame()

    if shape_file:
        for sensor, sensor_filter in zip(sensors, sensor_filters):
            # If we have a shape file, query the DB for scenes in that extent
            # TBD: The database geospatial/temporal query is currently Sentinel-1 only
            # GH issue: https://github.com/GeoscienceAustralia/gamma_insar/issues/261
            if sensor == "S1":
                # get the relative orbit number, which is int value of the numeric part of the track name
                # Note: This is S1 specific...
                rel_orbit = int(re.findall(r"\d+", str(proc_config.track))[0])

                # Find the maximum extent of the queried dates
                min_date = include_dates[0][0]
                max_date = max([d[1] for d in include_dates])

                # Query SLCs that match our search criteria for the maximum span
                # of dates that covers all of our include dates.
                slc_query_results = query_slc_inputs(
                    str(proc_config.database_path),
                    str(shape_file),
                    min_date,
                    max_date,
                    orbit,
                    rel_orbit,
                    polarisations,
                    sensor_filter
                )

                if slc_query_results is None:
                    continue

                slc_inputs_df = slc_inputs_df.append(
                    [slc_inputs(slc_query_results[pol]) for pol in polarisations],
                    ignore_index=True
                )

                # Filter out dates we don't care about - as our search query is for
                # a single giant span of dates, but our include dates may be more fine
                # grained than the query supports.
                exclude_indices = []

                for index, row in slc_inputs_df.iterrows():
                    date = row["date"]

                    keep = any(date >= lhs or date <= rhs for lhs,rhs in include_dates)
                    if not keep:
                        exclude_indices.append(index)

                slc_inputs_df.drop(exclude_indices, inplace=True)

            else:
                raise NotImplementedError(f"{sensor} does not support geospatial queries")

    # Add any explicit source data files into the "inputs" data frame
    slc_inputs_df = resolve_stack_scene_additional_files(
        slc_inputs_df,
        proc_config,
        polarisations,
        include_source_files,
        shape_file
    )

    # Convert to simpler date -> [source paths]
    if not slc_inputs_df.empty:
        for date in slc_inputs_df["date"].unique():
            urls = slc_inputs_df[slc_inputs_df["date"] == date].url.unique()
            urls = [Path(i) for i in urls.astype(str)]

            resolved_source_files.append((date, urls))

    # Note: returning both simplified and pandas data, both are useful in
    # different contexts / not sure it justifies two functions though
    return resolved_source_files, slc_inputs_df
