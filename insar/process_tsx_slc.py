import tempfile

from insar import constant
from insar.process_utils import convert

import structlog
from pathlib import Path
from dataclasses import dataclass

from insar.py_gamma_ga import GammaInterface, auto_logging_decorator, subprocess_wrapper
from insar.paths.slc import SlcPaths


# Customise Gamma shim to automatically handle basic error checking and logging
class ProcessSlcException(Exception):
    pass


_LOG = structlog.get_logger("insar")
pg = GammaInterface(
    subprocess_func=auto_logging_decorator(subprocess_wrapper, ProcessSlcException, _LOG)
)


# TODO: does there need to be any polarisation filtering?
#       Possibly - proc file specifies a polarisation to use...
#
# TODO: does this dataclass need to get paths from the main proc file for SLCs?

@dataclass
class TSXPaths:
    # output files
    slc: Path
    slc_par: Path
    sigma0_slc: Path
    sigma0_slc_par: Path
    sigma0_slc_bmp: Path
    sigma0_slc_png: Path

    @classmethod
    def create(cls, scene_date, dest_dir):
        dest_dir = Path(dest_dir) if isinstance(dest_dir, str) else dest_dir
        slc = (dest_dir / scene_date).with_suffix(".slc")
        slc_par = (dest_dir / scene_date).with_suffix(".slc.par")
        s0_slc = (dest_dir / scene_date).with_suffix(".sigma0.slc")
        s0_slc_par = (dest_dir / scene_date).with_suffix(".sigma0.slc.par")
        bmp = (dest_dir / scene_date).with_suffix(".sigma0.slc.bmp")
        png = (dest_dir / scene_date).with_suffix(".sigma0.slc.png")

        return TSXPaths(slc, slc_par, s0_slc, s0_slc_par, bmp, png)


def _verify_date_dir(dirs, product_path):
    if len(dirs) == 0:
        msg = f"No scene date directory found in {product_path}"
        raise ProcessSlcException(msg)
    elif len(dirs) > 1:
        msg = f"Multiple scene date dirs found in {product_path}:\n{dirs}"
        raise ProcessSlcException(msg)


def _verify_tsx_data_dirs(dirs, parent_path):
    if len(dirs) == 0:
        msg = f"No TSX data directory found in {parent_path}"
        raise ProcessSlcException(msg)
    elif len(dirs) > 1:
        msg = f"Multiple dirs found in {parent_path}:\n{dirs}"
        raise ProcessSlcException(msg)


def _verify_tsx_annotation_file(xmls, tsx_dir):
    if len(xmls) == 0:
        msg = f"No annotation XML file found in {tsx_dir}"
        raise ProcessSlcException(msg)
    elif len(xmls) > 1:
        msg = f"Multiple XML files found in {tsx_dir}, should only contain 1:\n{xmls}"
        raise ProcessSlcException(msg)


def _verify_cosar_file(cos_files, image_dir):
    if len(cos_files) == 0:
        msg = f"No COSAR file found in {image_dir}"
        raise ProcessSlcException(msg)
    elif len(cos_files) > 1:
        msg = f"Multiple COSAR files found in {image_dir}, should only contain 1:\n{cos_files}"
        raise ProcessSlcException(msg)


def process_tsx_slc(
    product_path: Path,
    output_dir: Path,
    slc_paths: SlcPaths = None,
):
    """
    Processes TSX/TDX data into GAMMA SLC data.

    :param product_path:
        Path to the TSX product directory. This is the 2ND subdir with all the files:
         e.g. "TDX1_SAR__SSC______SM_S_SRA_20170411T192821_20170411T192829"
    :param output_dir:
        Dir path to write outputs to.
    :param slc_paths:
        Optional SlcPaths if the default file locations need to be changed.
    """
    # if not product_path.exists():
    #     raise RuntimeError(f"The provided product path does not exist! product_path={product_path}")
    #
    # if not output_dir.is_dir():
    #     raise RuntimeError("output_dir is a file...")
    #
    # if not output_dir.exists():
    #     raise RuntimeError(f"The provided output dir path does not exist! output_dir={output_dir}")

    _LOG.info(f"process_tsx_slc(): product_path={product_path}")

    # scene_date = product_path.name
    # tsx_paths = tsx_paths if tsx_paths else TSXPaths.create(scene_date, output_dir)
    #
    # # find second scene date from root dir in the .tar.gz archive
    # scene_date_pattern = "[0-9]" * 8
    # date_dirs = list(product_path.glob(scene_date_pattern))
    # _verify_date_dir(date_dirs, product_path)
    # date_dir = date_dirs[0]  # is a Path
    #
    # # find long TSX dir under the "root" date dir
    # pattern = "T[SD]X[0-9]_SAR_*T[0-9][0-9][0-9][0-9][0-9][0-9]"
    # dirs = list(date_dir.glob(pattern))
    # _verify_tsx_data_dirs(dirs, date_dir)
    # tsx_dir = dirs[0]  # big ugly dir below the "root" date dir

    # find the annotation XML file (duplicates the big ugly name & adds .xml suffix
    pattern = "T[SD]X[0-9]_SAR_*T[0-9][0-9][0-9][0-9][0-9][0-9]"
    xmls = list(product_path.glob(pattern + ".xml"))
    _verify_tsx_annotation_file(xmls, product_path)
    xml_meta = xmls[0]

    # locate the image data file
    image_dir = product_path / "IMAGEDATA"
    cos_files = list(image_dir.glob("IMAGE_*.cos"))
    _verify_cosar_file(cos_files, image_dir)
    cosar = cos_files[0]

    with tempfile.tempdir() as td:
        base_name = slc_paths.slc.name
        gamma_slc = Path(base_name + ".slc")
        gamma_slc_par = Path(base_name + ".slc.par")

        # Read TSX data and produce SLC and parameter files in GAMMA format
        pg.par_TX_SLC(xml_meta,  # TSX product annotation file
                      cosar,  # COSAR SSC stripmap
                      gamma_slc_par,  # output param file
                      gamma_slc,  # output SLC data
        )

        # Apply stated calFactor from the xml file, scale according to sin(inc_angle) and
        # convert from scomplex to fcomplex. Output is sigma0
        # This does the file MOVE step inline compared to BASH version
        pg.radcal_SLC(gamma_slc,
                      gamma_slc_par,
                      slc_paths.slc,  # SLC output file, the sigma0
                      slc_paths.slc_par,  # SLC PAR output file, the sigma0 par
                      3,  # fcase: scomplex --> fcomplex
                      constant.NOT_PROVIDED,  # antenna gain file
                      0,  # rloss_flag
                      0,  # ant_flag
                      1,  # refarea_flag
                      0,  # sc_dB
                      constant.NOT_PROVIDED,  # K_dB
        )

        # Make quick-look png image of SLC
        par = pg.ParFile(slc_paths.slc_par.as_posix())
        width = par.get_value("range_samples", dtype=int, index=0)
        lines = par.get_value("azimuth_lines", dtype=int, index=0)
        bmp_path = slc_paths.slc.with_suffix("bmp")
        png_path = slc_paths.slc.with_suffix("png")

        pg.rasSLC(slc_paths.slc,
                  width,
                  1,
                  lines,
                  50,
                  20,
                  constant.NOT_PROVIDED,
                  constant.NOT_PROVIDED,
                  1,
                  0,
                  0,
                  bmp_path,
        )

        convert(bmp_path, png_path)
        bmp_path.unlink()
