import tempfile
import shutil
from pathlib import Path
from unittest import mock
import pytest
from PIL import Image

from tests.py_gamma_test_proxy import PyGammaTestProxy

import insar.coregister_slc
from insar.coregister_slc import CoregisterSlc
from insar.project import ProcConfig, IfgFileNames, DEMFileNames


def get_test_context():
    temp_dir = tempfile.TemporaryDirectory()
    data_dir = Path(temp_dir.name) / '20151127'

    pgp = PyGammaTestProxy()
    pgp = mock.Mock(spec=PyGammaTestProxy, wraps=pgp)

    # Make offset_fit return parseable stdout as required for coregister_slc to function
    def offset_fit_sideffect(offs: str, ccp: str, OFF_par: str, coffs: str, coffsets: str, thres, npoly, interact_flag):
        shutil.copyfile(data_dir / 'offset_fit.start', OFF_par)
        return 0, 'final model fit std. dev. (samples) range:   0.3699  azimuth:   0.1943', ''

    # raspwr needs to create a dummy bmp
    def raspwr_sideffect(pwr: str, width, start, nlines, pixavr, pixavaz, scale, exp, LR, rasf: str, data_type = None, hdrz = None):
        slave_gamma0_eqa = data_dir / rasf
        Image.new('RGB', size=(50, 50), color=(155, 0, 0)).save(slave_gamma0_eqa)
        return 0, '', ''

    pgp.raspwr.side_effect = raspwr_sideffect
    pgp.raspwr.return_value = 0, '', ''

    pgp.offset_fit.side_effect = offset_fit_sideffect
    pgp.offset_fit.return_value = (0, 'final model fit std. dev. (samples) range:   0.3699  azimuth:   0.1943', '')

    # Copy test data
    shutil.copytree(Path(__file__).parent.absolute() / 'data' / '20151127', data_dir)

    data = {
        'proc': None,  # FIXME: We need to get a valid gamma proc file
        'slc_master': data_dir / '20151127_VV.slc',
        'slc_slave': data_dir / '20151127_VV.slc',  # if slave/master are the same... everything should still run i assume? just useless outputs?
        'slave_mli': data_dir / '20151127_VV_8rlks.mli',
        'range_looks': 1,
        'azimuth_looks': 1,
        'ellip_pix_sigma0': data_dir / 'TODO_ellip_pix_sigma0',
        'dem_pix_gamma0': data_dir / '20151127_VV_8rlks_rdc_pix_gamma0',
        'r_dem_master_mli': data_dir / 'r20151127_VV_8rlks.mli',  # FIXME: Pretty sure this isn't right (but can't find any DEM runs with an mli)
        'rdc_dem': data_dir / '20151127_VV_8rlks_rdc.dem',
        'eqa_dem_par': data_dir / '20180127_VV_8rlks_eqa.dem.par',  # HACK: This is from a totally different scene
        'dem_lt_fine': data_dir / '20151127_VV_8rlks_eqa_to_rdc.lt',
    }

    return pgp, data, temp_dir


def test_valid_data(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_slc, 'pg', pgp)

    with temp_dir as temp_path:
        coreg = CoregisterSlc(
            *data.values(),
            Path(temp_path)
        )

        assert(str(coreg.out_dir) == temp_path)

        coreg.main()

        # Assert no failure status for any gamma call
        # assert()

        # Assert outputs exist
        #assert(Path(f"{coreg.master_slave_prefix}.ovr_results").exists())

        # Assert quick-look images exist
        slave_png = coreg.out_dir.joinpath(f"{coreg.slave_gamma0_eqa.name}.png")
        assert(slave_png.exists())

        # TODO: Assert commands/stats in pgp


def test_set_tab_files(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_slc, 'pg', pgp)

    with temp_dir as temp_path:
        coreg = CoregisterSlc(
            *data.values(),
            Path(temp_path)
        )

        coreg.set_tab_files()
        assert(coreg.slave_slc_tab.exists())
        assert(coreg.r_slave_slc_tab.exists())
        assert(coreg.master_slc_tab.exists())

        custom_dir = 'test123abc'
        (Path(temp_path) / custom_dir).mkdir()
        coreg.set_tab_files(Path(temp_path) / custom_dir)
        assert(coreg.slave_slc_tab.exists() and coreg.slave_slc_tab.parent.name == custom_dir)
        assert(coreg.r_slave_slc_tab.exists() and coreg.r_slave_slc_tab.parent.name == custom_dir)
        assert(coreg.master_slc_tab.exists() and coreg.master_slc_tab.parent.name == custom_dir)


def test_get_lookup(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_slc, 'pg', pgp)

    with temp_dir as temp_path:
        coreg = CoregisterSlc(
            *data.values(),
            Path(temp_path)
        )

        # Create dummy inputs that are expected
        coreg.r_dem_master_mli_par.touch()
        coreg.rdc_dem.touch()
        coreg.slave_mli_par.touch()

        # Run function
        coreg.get_lookup()

        # Ensure the output is produced
        assert(coreg.slave_lt.exists())


def test_resample_full(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_slc, 'pg', pgp)

    with temp_dir as temp_path:
        coreg = CoregisterSlc(
            *data.values(),
            Path(temp_path)
        )

        coreg.set_tab_files()
        coreg.resample_full()

        assert(Path(coreg.r_slave_slc_tab).exists())
        assert(Path(coreg.r_slave_slc).exists())
        assert(Path(coreg.r_slave_slc_par).exists())


def test_multi_look(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_slc, 'pg', pgp)

    with temp_dir as temp_path:
        coreg = CoregisterSlc(
            *data.values(),
            Path(temp_path)
        )

        coreg.multi_look()

        assert(Path(coreg.r_slave_mli).exists())
        assert(Path(coreg.r_slave_mli_par).exists())


def test_generate_normalised_backscatter(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_slc, 'pg', pgp)

    with temp_dir as temp_path:
        coreg = CoregisterSlc(
            *data.values(),
            Path(temp_path)
        )

        coreg.generate_normalised_backscatter()

        slave_gamma0 = coreg.out_dir / f"{coreg.slave_mli.stem}.gamma0"
        slave_gamma0_eqa = coreg.out_dir / f"{coreg.slave_mli.stem}_eqa.gamma0"
        slave_png = coreg.out_dir / f"{slave_gamma0_eqa.name}.png"

        assert(slave_gamma0.exists())
        assert(slave_gamma0_eqa.exists())
        assert(slave_png.exists())

        assert(slave_gamma0_eqa.with_suffix(".gamma0.tif").exists())

        assert(slave_gamma0_eqa.with_suffix(".sigma0").exists())
        assert(slave_gamma0_eqa.with_suffix(".sigma0.tif").exists())


# TODO: Test more specific corner cases (what are they?)
