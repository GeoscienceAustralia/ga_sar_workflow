import tempfile
import shutil
from pathlib import Path
from unittest import mock
import pytest
from PIL import Image

from tests.py_gamma_test_proxy import PyGammaTestProxy

import insar.coregister_dem
from insar.coregister_dem import CoregisterDem
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

    def SLC_copy_se(SLC_in: str, SLC_par_in: str, SLC_out: str, SLC_par_out: str, *args):
        #shutil.copyfile(SLC_in, SLC_out)
        shutil.copyfile(SLC_par_in, SLC_par_out)
        print(SLC_par_in, Path(SLC_par_in).stat().st_size)
        print(SLC_par_out, Path(SLC_par_out).stat().st_size)
        return 0, '', ''

    def multi_look_se(MLI_in: str, MLI_in_par: str, MLI_out: str, MLI_out_par: str, rlks, azlks, *args):
        shutil.copyfile(MLI_in_par, MLI_out_par)
        return 0, '', ''

    def gc_map1_se(MLI_par: str, OFF_par: str, DEM_par: str, DEM: str, DEM_seg_par: str, DEM_seg: str, lookup_table: str, *args):
        shutil.copyfile(DEM_par, DEM_seg_par)
        Path(DEM_seg).touch()
        return 0, '', ''

    def rashgt_se(hgt: str, pwr: str, width, start_hgt = None, start_pwr = None, nlines = None, pixavr = None, pixavaz = None, m_cycle = None, scale = None, exp = None, LR = None, rasf: str = None):
        Image.new('RGB', size=(50, 50), color=(155, 0, 0)).save(rasf)
        return 0, '', ''

    pgp.raspwr.side_effect = raspwr_sideffect
    pgp.raspwr.return_value = 0, '', ''

    pgp.offset_fit.side_effect = offset_fit_sideffect
    pgp.offset_fit.return_value = 0, 'final model fit std. dev. (samples) range:   0.3699  azimuth:   0.1943', ''

    pgp.SLC_copy.side_effect = SLC_copy_se
    pgp.SLC_copy.return_value = 0, '', ''

    pgp.multi_look.side_effect = multi_look_se
    pgp.multi_look.return_value = 0, '', ''

    pgp.gc_map1.side_effect = gc_map1_se
    pgp.gc_map1.return_value = 0, '', ''

    pgp.rashgt.side_effect = rashgt_se
    pgp.rashgt.return_value = 0, '', ''

    # Copy test data
    shutil.copytree(Path(__file__).parent.absolute() / 'data' / '20151127', data_dir)

    print('data_dir', data_dir)

    data = {
        'rlks': 8,
        'alks': 8,
        'dem': data_dir / '20180127_VV_8rlks_eqa.dem',
        'slc': data_dir / '20151127_VV.slc',
        'dem_par': data_dir / '20180127_VV_8rlks_eqa.dem.par',
        'slc_par': data_dir / '20151127_VV.slc.par',
        'dem_patch_window': 1024,
        'dem_rpos': None,
        'dem_azpos': None,
        'dem_offset': (0, 0),
        'dem_offset_measure': (32, 32),
        'dem_window': (256, 256),
        'dem_snr': 0.15,
        'dem_rad_max': 4,
        'dem_ovr': 1,
    }

    return pgp, data, temp_dir


def test_valid_data(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_dem, 'pg', pgp)
    monkeypatch.setattr(insar.coregister_dem, 'run_command', lambda cmds, cwd: print('run_command', cmds, 'in dir', cwd))

    with temp_dir as temp_path:
        out_dir = Path(temp_path)
        slc_outdir = out_dir / '20151127'
        dem_outdir = out_dir / '20151127' / 'dem'

        coreg = CoregisterDem(
            *data.values(),
            dem_outdir,
            slc_outdir
        )

        assert(str(coreg.dem_outdir) == str(dem_outdir))

        coreg.main()

        for path in dem_outdir.iterdir():
            print(path)

        dem_filenames = coreg.dem_filenames(dem_prefix=f"{coreg.slc.stem}_{coreg.rlks}rlks", outdir=coreg.dem_outdir)
        for name, path in dem_filenames.items():
            if path.suffix == '.dem':
                assert(Path(path).exists())


def test_copy_slc(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_dem, 'pg', pgp)

    with temp_dir as temp_path:
        out_dir = Path(temp_path)
        slc_outdir = out_dir / '20151127'
        dem_outdir = out_dir / '20151127' / 'dem'

        coreg = CoregisterDem(
            *data.values(),
            dem_outdir,
            slc_outdir
        )

        coreg.copy_slc(False)

        assert(Path(coreg.r_dem_master_mli).exists())
        assert(Path(coreg.r_dem_master_mli_par).exists())
        assert(not Path(coreg.r_dem_master_mli_bmp).exists())

        coreg.copy_slc(True)

        assert(Path(coreg.r_dem_master_mli).exists())
        assert(Path(coreg.r_dem_master_mli_par).exists())
        assert(Path(coreg.r_dem_master_mli_bmp).exists())


def test_gen_dem_rdc(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_dem, 'pg', pgp)

    with temp_dir as temp_path:
        out_dir = Path(temp_path)
        slc_outdir = out_dir / '20151127'
        dem_outdir = out_dir / '20151127' / 'dem'

        coreg = CoregisterDem(
            *data.values(),
            dem_outdir,
            slc_outdir
        )

        dem_outdir.mkdir(parents=True, exist_ok=True)
        slc_outdir.mkdir(parents=True, exist_ok=True)

        coreg.copy_slc()

        coreg.gen_dem_rdc(False)

        assert(Path(coreg.dem_pix_gam).exists())
        assert(not Path(coreg.ext_image_flt).exists())
        assert(not Path(coreg.ext_image_init_sar).exists())

        # FIXME: This branch seems to not be implemented fully/correctly (nothing ever sets ext_image)
        # coreg.gen_dem_rdc(True)

        # assert(Path(coreg.dem_pix_gam).exists())
        # assert(Path(coreg.ext_image_flt).exists())
        # assert(Path(coreg.ext_image_init_sar).exists())


def test_create_diff_par(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_dem, 'pg', pgp)

    with temp_dir as temp_path:
        outdir = Path(temp_path) / '20151127'
        outdir.mkdir(parents=True, exist_ok=True)

        coreg = CoregisterDem(
            *data.values(),
            outdir,
            outdir
        )

        coreg.geocode()

        # create_diff_par runs an external command to generate the diff, we create it explicitly...
        # which defeats the purpose of the test right now - but when it's migrated to pg.create_diff_par
        # we can remove this hack.  At the very least this is a smoke test for CoregisterDemcreate_diff_par
        monkeypatch.setattr(insar.coregister_dem, 'run_command', lambda cmds, cwd: Path(coreg.dem_diff).touch())

        coreg.create_diff_par()

        assert(Path(coreg.dem_diff).exists())


def test_offset_calc(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_dem, 'pg', pgp)

    with temp_dir as temp_path:
        out_dir = Path(temp_path)
        slc_outdir = out_dir / '20151127'
        dem_outdir = out_dir / '20151127' / 'dem'

        coreg = CoregisterDem(
            *data.values(),
            dem_outdir,
            slc_outdir
        )

        dem_outdir.mkdir(parents=True, exist_ok=True)
        slc_outdir.mkdir(parents=True, exist_ok=True)

        coreg.copy_slc()
        coreg.gen_dem_rdc()
        coreg.offset_calc()

        assert(Path(coreg.dem_lt_fine).exists())
        assert(Path(coreg.dem_master_gamma0).exists())
        assert(Path(coreg.dem_master_gamma0_bmp).exists())
        assert(Path(coreg.seamask).exists())


def test_geocode(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_dem, 'pg', pgp)

    with temp_dir as temp_path:
        outdir = Path(temp_path) / '20151127'
        outdir.mkdir(parents=True, exist_ok=True)

        coreg = CoregisterDem(
            *data.values(),
            outdir,
            outdir
        )

        coreg.geocode()

        assert(Path(coreg.rdc_dem).exists())
        assert(Path((coreg.dem_outdir / coreg.rdc_dem).with_suffix(".png")).exists())
        assert(Path(coreg.dem_rdc_sim_sar).exists())
        assert(Path(coreg.dem_rdc_inc).exists())
        assert(Path(coreg.dem_master_gamma0_eqa).exists())
        assert(Path(coreg.dem_master_gamma0_eqa_bmp).exists())
        assert(Path(coreg.dem_master_gamma0_eqa_bmp.with_suffix(".png")).exists())
        assert(Path(coreg.dem_master_gamma0_eqa_geo).exists())
        assert(Path(coreg.dem_master_sigma0_eqa).exists())
        assert(Path(coreg.dem_master_sigma0_eqa_geo).exists())

        assert(Path(coreg.dem_master_gamma0_eqa_bmp.with_suffix(".kml")).exists())


def test_look_vector(monkeypatch):
    pgp, data, temp_dir = get_test_context()
    monkeypatch.setattr(insar.coregister_dem, 'pg', pgp)

    with temp_dir as temp_path:
        outdir = Path(temp_path) / '20151127'

        outdir.mkdir(parents=True, exist_ok=True)

        coreg = CoregisterDem(
            *data.values(),
            outdir,
            outdir
        )

        coreg.look_vector()

        assert(Path(coreg.dem_lv_theta).exists())
        assert(Path(coreg.dem_lv_phi).exists())
        assert(Path(coreg.dem_lv_phi_geo).exists())


# TODO: Test more specific corner cases (what are they?)
