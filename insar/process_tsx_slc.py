from insar import constant
from insar.process_utils import convert

import structlog
from pathlib import Path
from dataclasses import dataclass

from insar.py_gamma_ga import GammaInterface, auto_logging_decorator, subprocess_wrapper


# Customise Gamma shim to automatically handle basic error checking and logging
class ProcessSlcException(Exception):
    pass


_LOG = structlog.get_logger("insar")
pg = GammaInterface(
    subprocess_func=auto_logging_decorator(subprocess_wrapper, ProcessSlcException, _LOG)
)


# TODO: does there need to be any polarisation filtering? Possibly - proc file specifies a polarisation to use...

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


def _verify_tsx_data_dirs(dirs, product_path):
    if len(dirs) == 0:
        msg = f"No TSX data directory found in {product_path}"
        raise ProcessSlcException(msg)
    elif len(dirs) > 1:
        msg = f"Multiple dirs found in {product_path}:\n{dirs}"
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
    tsx_paths: TSXPaths = None,
):
    """
    Processes TSX/TDX data into GAMMA SLC data.

    :param product_path:
        Path to the TSX product directory. This is the dir with the 8 digit YYYYMMDD timestamp.
    :param output_dir:
        Dir path to write outputs to.
    :param tsx_paths:
        Optional TSXPaths if the default file locations need to be changed.
    """
    if not product_path.exists():
        raise RuntimeError("The provided product path does not exist!")

    if not output_dir.exists():
        raise RuntimeError("The provided output dir path does not exist!")

    scene_date = product_path.name
    tsx_paths = tsx_paths if tsx_paths else TSXPaths.create(scene_date, output_dir)

    # find the long TSX dir under the "root" date dir
    pattern = "T[SD]X[0-9]_SAR_*T[0-9][0-9][0-9][0-9][0-9][0-9]"
    dirs = list(product_path.glob(pattern))
    _verify_tsx_data_dirs(dirs, product_path)
    tsx_dir = dirs[0]  # the big ugly dir below the "root" date dir

    # find the annotation XML file (duplicates the big ugly name & adds .xml suffix
    xmls = list(tsx_dir.glob(pattern + ".xml"))
    _verify_tsx_annotation_file(xmls, tsx_dir)
    xml_meta = xmls[0]

    # locate the image data file
    image_dir = tsx_dir / "IMAGEDATA"
    cos_files = list(image_dir.glob("IMAGE_*.cos"))
    _verify_cosar_file(cos_files, image_dir)
    cosar = cos_files[0]

    # Read TSX data and produce SLC and parameter files in GAMMA format
    pg.par_TX_SLC(xml_meta,  # TSX product annotation file
                  cosar,  # COSAR SSC stripmap
                  tsx_paths.slc_par,  # output param file
                  tsx_paths.slc,  # output SLC data
    )

    # Apply stated calFactor from the xml file, scale according to sin(inc_angle) and
    # convert from scomplex to fcomplex. Output is sigma0
    pg.radcal_SLC(tsx_paths.slc,
                  tsx_paths.slc_par,
                  tsx_paths.sigma0_slc,  # SLC output file
                  tsx_paths.sigma0_slc_par,  # SLC PAR output file
                  3,  # fcase: scomplex --> fcomplex
                  constant.NOT_PROVIDED,  # antenna gain file
                  0,  # rloss_flag
                  0,  # ant_flag
                  1,  # refarea_flag
                  0,  # sc_dB
                  constant.NOT_PROVIDED,  # K_dB
    )

    # Make quick-look png image of SLC
    par = pg.ParFile(tsx_paths.sigma0_slc_par.as_posix())
    width = par.get_value("range_samples", dtype=int, index=0)
    lines = par.get_value("azimuth_lines", dtype=int, index=0)

    pg.rasSLC(tsx_paths.sigma0_slc,
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
              tsx_paths.sigma0_slc_bmp,
    )

    convert(tsx_paths.sigma0_slc_bmp, tsx_paths.sigma0_slc_png)
    tsx_paths.sigma0_slc_bmp.unlink()
