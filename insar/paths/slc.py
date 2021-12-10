from typing import Optional, List, Union
from pathlib import Path

from insar.project import ProcConfig
from insar.stack import load_stack_config

class SlcPaths:
    """
    Contains the path names of key SLC data and related files.

    All code should use this class when referring to pathnames relating to SLC data to
    avoid duplicating/repeating pathnames to avoid refactoring/renaming related errors.
    """

    dir: Path

    slc: Path
    slc_par: Path
    slc_tops_par: Path
    slc_tab: Path

    iw_slc: List[Path]
    iw_slc_par: List[Path]
    iw_slc_tops_par: List[Path]

    mli: Optional[Path]
    mli_par: Optional[Path]

    def __init__(
        self,
        stack_config: Union[ProcConfig, Path],
        date: str,
        pol: str,
        rlks: Optional[int] = None
    ):
        """
        TODO: Documentation
        """

        if not isinstance(stack_config, ProcConfig):
            stack_config = load_stack_config(stack_config)

        # TODO: Should be based on some kind of 'stack' paths?
        slc_dir = Path(stack_config.output_path) / stack_config.slc_dir / date
        self.dir = slc_dir

        self.slc = slc_dir / f"{date}_{pol}.slc"
        self.slc_par = slc_dir / f"{date}_{pol}.slc.par"
        self.slc_tops_par = slc_dir / f"{date}_{pol}.slc.TOPS_par"
        self.slc_tab = slc_dir / f"{date}_{pol}_tab"

        swaths = [1, 2, 3]
        self.iw_slc = [slc_dir / f"{date}_{pol}_IW{i}.slc" for i in swaths]
        self.iw_slc_par = [slc_dir / f"{date}_{pol}_IW{i}.slc.par" for i in swaths]
        self.iw_slc_tops_par = [slc_dir / f"{date}_{pol}_IW{i}.slc.TOPS_par" for i in swaths]

        if rlks is not None:
            self.mli = slc_dir / f"{date}_{pol}_{rlks}rlks.mli"
            self.mli_par = slc_dir / f"{date}_{pol}_{rlks}rlks.mli.par"
