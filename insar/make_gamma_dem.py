#!/usr/bin/env python3

from os.path import join as pjoin
from pathlib import Path
from typing import Optional, Union
import tempfile

import shapely
import shapely.wkt
from shapely.ops import cascaded_union
import rasterio
from spatialist import Vector
import structlog

from insar import constant as const
from insar.subprocess_utils import run_command
from insar.py_gamma_ga import GammaInterface, auto_logging_decorator, subprocess_wrapper

_LOG = structlog.get_logger("insar")


class GammaDemException(Exception):
    pass


pg = GammaInterface(
    subprocess_func=auto_logging_decorator(subprocess_wrapper, GammaDemException, _LOG)
)


def create_gamma_dem(
    gamma_dem_dir: Union[Path, str],
    dem_img: Union[Path, str],
    track_frame: str,
    shapefile: Union[Path, str],
    buffer_width: Optional[float] = 0.3,
    create_png: Optional[bool] = False,
) -> None:
    """
    Automatically creates a DEM and par file for use with GAMMA.
    :param gamma_dem_dir:
        A directory to where gamma dem and par file will be written
    :param dem_img:
        A DEM from where gamma dem will be extracted from
    :param track_frame:
        A track and frame name
    :param shapefile:
        A 'Path' to a shapefile
    :param buffer_width:
        Additional buffer to include in a subset of shapefile extent
    :param create_png:
        A flag to create preview of dem
    """
    vector_object = Vector(shapefile)

    def _get_bounds():
        if isinstance(vector_object, Vector):
            return (
                cascaded_union(
                    [
                        shapely.wkt.loads(extent)
                        for extent in vector_object.convert2wkt(set3D=False)
                    ]
                )
                .buffer(buffer_width, cap_style=2, join_style=2)
                .bounds
            )

    with tempfile.TemporaryDirectory() as tmpdir:
        min_lon, min_lat, max_lon, max_lat = _get_bounds()

        # subset dem_img for the area of interest as a geotiff format
        outfile = "{}_temp.tif".format(track_frame)
        command = [
            "gdal_translate",
            "-projwin",
            str(min_lon),
            str(max_lat),
            str(max_lon),
            str(min_lat),
            dem_img,
            outfile,
        ]
        run_command(command, tmpdir)

        # get nodata value from geotiff file to make into gamma compatible format
        # sets the nodata pixels to 0.0001 and update nodata value to None
        with rasterio.open(pjoin(tmpdir, outfile), "r") as src:
            data = src.read(1)
            mask = data == src.nodata
            data[mask] = 0.0001
            profile = src.profile
            profile.update(nodata=None)
            outfile_new = pjoin(tmpdir, "{}.tif".format(track_frame))

            with rasterio.open(outfile_new, "w", **profile) as dst:
                dst.write(data, 1)

            # dem_import function is a call to GAMMA software which creates a gamma compatible DEM
            input_dem_pathname = outfile_new
            dem_pathname = str(Path(gamma_dem_dir).joinpath(f"{track_frame}.dem"))
            par_pathname = f"{dem_pathname}.par"
            input_type = 0  # GeoTIFF
            priority = 1  # input DEM parameters have priority

            pg.dem_import(
                input_dem_pathname,
                dem_pathname,
                par_pathname,
                input_type,
                priority,
                const.NOT_PROVIDED,  # geoid or constant geoid height value
                const.NOT_PROVIDED,  # geoid_par, geoid DEM_par file
                const.NOT_PROVIDED,  # geoid_type, global geoid in EQA coordinates
                const.NOT_PROVIDED,  # lat_n_shift, latitude or Northing constant shift to apply
                const.NOT_PROVIDED,  # lon_e_shift, longitude or Easting constant shift
                const.NOT_PROVIDED,  # zflg, no_data values in input file are kept in output file
                const.NOT_PROVIDED,  # no_data, value defined in input metadata
            )

            # TODO: replace/add PIL option?
            # create a preview of dem file
            if create_png:
                dem_png_file = "{}_dem_preview.png".format(track_frame)
                command = [
                    "gdal_translate",
                    "-of",
                    "PNG",
                    "-outsize",
                    "10%",
                    "10%",
                    outfile_new,
                    dem_png_file,
                ]
                run_command(command, gamma_dem_dir)
