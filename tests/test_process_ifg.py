import pathlib
import subprocess
from unittest import mock

import insar.constant as const
from insar import process_ifg, py_gamma_ga
from insar.process_ifg import ProcessIfgException
from insar.project import ProcConfig, IfgFileNames, DEMFileNames

import structlog
import pytest


# FIXME: tweak settings to ensure working dir doesn't have to be changed for INT processing (do in workflow)
# FIXME: change all mocks to return (return_code, cout, cerr) as per decorator

# TODO: can monkeypatch be done at higher level scope to apply to multiple test funcs?
@pytest.fixture
def pg_int_mock():
    """Create basic mock of the py_gamma module for INT processing step."""
    pg_mock = mock.Mock()
    pg_mock.create_offset.return_value = 0
    pg_mock.offset_pwr.return_value = 0
    pg_mock.offset_fit.return_value = 0
    pg_mock.create_diff_par.return_value = 0
    return pg_mock


@pytest.fixture
def pc_mock():
    """Returns basic mock to simulate a ProcConfig object."""
    pc = mock.Mock(spec=ProcConfig)
    return pc


@pytest.fixture
def ic_mock():
    """Returns basic mock to simulate an IfgFileNames object."""
    ic = mock.Mock(spec=IfgFileNames)
    ic.ifg_bperp = mock.MagicMock(spec=pathlib.Path)
    return ic


def test_calc_int(monkeypatch, pg_int_mock, pc_mock, ic_mock):
    """Verify default path through the INT processing step without cleanup."""

    # craftily substitute the 'pg' py_gamma obj for a mock: avoids a missing import when testing
    # locally, or calling the real thing on Gadi...
    # TODO: monkeypatch or use unittest.mock's patch? Which is better?
    monkeypatch.setattr(process_ifg, "pg", pg_int_mock)

    ic_mock.ifg_off = mock.Mock(spec=pathlib.Path)
    ic_mock.ifg_off.exists.return_value = False  # offset not yet processed

    process_ifg.calc_int(pc_mock, ic_mock, clean_up=False)

    assert pg_int_mock.create_offset.called

    # ensure CSK sensor block / SP mode section is skipped
    assert pg_int_mock.init_offset_orbit.called is False
    assert pg_int_mock.init_offset.called is False

    # check the core processing was called
    assert pg_int_mock.offset_pwr.called
    assert pg_int_mock.offset_fit.called
    assert pg_int_mock.create_diff_par.called


def test_calc_int_with_cleanup(monkeypatch, pg_int_mock, pc_mock, ic_mock):
    monkeypatch.setattr(process_ifg, "pg", pg_int_mock)

    ic_mock.ifg_off = mock.Mock(spec=pathlib.Path)
    ic_mock.ifg_off.exists.return_value = False  # offset not yet processed

    ic_mock.ifg_offs = mock.Mock(spec=pathlib.Path)
    ic_mock.ifg_ccp = mock.Mock(spec=pathlib.Path)
    ic_mock.ifg_coffs = mock.Mock(spec=pathlib.Path)
    ic_mock.ifg_coffsets = mock.Mock(spec=pathlib.Path)

    assert ic_mock.ifg_offs.unlink.called is False
    assert ic_mock.ifg_ccp.unlink.called is False
    assert ic_mock.ifg_coffs.unlink.called is False
    assert ic_mock.ifg_coffsets.unlink.called is False

    process_ifg.calc_int(pc_mock, ic_mock, clean_up=True)

    assert ic_mock.ifg_offs.unlink.called
    assert ic_mock.ifg_ccp.unlink.called
    assert ic_mock.ifg_coffs.unlink.called
    assert ic_mock.ifg_coffsets.unlink.called


def test_error_handling_decorator(monkeypatch):
    # force all fake subprocess calls to fail
    fake_subprocess = mock.Mock(return_value=-1)

    pgi = py_gamma_ga.GammaInterface(
        install_dir="./fake-install",
        gamma_exes={"create_offset": "fake-EXE-name"},
        subprocess_func=process_ifg.decorator(fake_subprocess),
    )

    # ensure mock logger has all core error(), msg() etc logging functions
    log_mock = mock.Mock(spec=structlog.stdlib.BoundLogger)
    assert log_mock.error.called is False
    monkeypatch.setattr(process_ifg, "_LOG", log_mock)

    with pytest.raises(ProcessIfgException):
        pgi.create_offset(1, 2, 3, key="value")

    assert log_mock.error.called
    has_cout = has_cerr = False

    for c in log_mock.error.call_args:
        if const.COUT in c:
            has_cout = True

        if const.CERR in c:
            has_cerr = True

    assert has_cout
    assert has_cerr


@pytest.fixture
def pg_flat_mock():
    """Create basic mock of the py_gamma module for the INT processing step."""
    pg_mock = mock.Mock()
    ret = (0, "cout-fake-content", "cerr-fake-content")
    pg_mock.base_orbit.return_value = ret
    pg_mock.phase_sim_orb.return_value = ret
    pg_mock.SLC_diff_intf.return_value = ret
    pg_mock.base_init.return_value = ret
    pg_mock.base_add.return_value = ret
    pg_mock.phase_sim.return_value = ret

    pg_mock.gcp_phase.return_value = ret
    pg_mock.sub_phase.return_value = ret
    pg_mock.mcf.return_value = ret
    pg_mock.base_ls.return_value = ret
    pg_mock.cc_wave.return_value = ret
    pg_mock.rascc_mask.return_value = ret
    pg_mock.multi_cpx.return_value = ret
    pg_mock.multi_real.return_value = ret
    pg_mock.base_perp.return_value = ret
    pg_mock.extract_gcp.return_value = ret
    return pg_mock


@pytest.fixture
def dc_mock():
    """Default mock for DEMFileNames config."""
    dcm = mock.Mock(spec=DEMFileNames)
    return dcm


def test_generate_init_flattened_ifg(
    monkeypatch, pg_flat_mock, pc_mock, ic_mock, dc_mock
):
    monkeypatch.setattr(process_ifg, "pg", pg_flat_mock)

    assert pg_flat_mock.base_orbit.called is False
    assert pg_flat_mock.phase_sim_orb.called is False
    assert pg_flat_mock.SLC_diff_intf.called is False
    assert pg_flat_mock.base_init.called is False
    assert pg_flat_mock.base_add.called is False
    assert pg_flat_mock.phase_sim.called is False

    process_ifg.generate_init_flattened_ifg(pc_mock, ic_mock, dc_mock, clean_up=False)

    assert pg_flat_mock.base_orbit.called
    assert pg_flat_mock.phase_sim_orb.called
    assert pg_flat_mock.SLC_diff_intf.call_count == 2
    assert pg_flat_mock.base_init.called
    assert pg_flat_mock.base_add.called
    assert pg_flat_mock.phase_sim.called


def test_generate_final_flattened_ifg(
    monkeypatch, pg_flat_mock, pc_mock, ic_mock, dc_mock
):
    # test refinement of baseline model using ground control points
    monkeypatch.setattr(process_ifg, "pg", pg_flat_mock)

    assert pg_flat_mock.multi_cpx.called is False
    assert pg_flat_mock.cc_wave.called is False
    assert pg_flat_mock.rascc_mask.called is False
    assert pg_flat_mock.mcf.called is False
    assert pg_flat_mock.multi_real.called is False
    assert pg_flat_mock.sub_phase.called is False
    assert pg_flat_mock.extract_gcp.called is False
    assert pg_flat_mock.gcp_phase.called is False
    assert pg_flat_mock.base_ls.called is False
    assert pg_flat_mock.phase_sim.called is False
    assert pg_flat_mock.SLC_diff_intf.called is False
    assert pg_flat_mock.base_perp.called is False

    width10, ifg_width = 101, 99  # fake
    process_ifg.generate_final_flattened_ifg(
        pc_mock, ic_mock, dc_mock, width10, ifg_width, clean_up=False
    )

    assert pg_flat_mock.multi_cpx.called
    assert pg_flat_mock.cc_wave.call_count == 3
    assert pg_flat_mock.rascc_mask.call_count == 2
    assert pg_flat_mock.mcf.called
    assert pg_flat_mock.multi_real.called
    assert pg_flat_mock.sub_phase.call_count == 2
    assert pg_flat_mock.extract_gcp.called
    assert pg_flat_mock.gcp_phase.called
    assert pg_flat_mock.base_ls.called
    assert pg_flat_mock.phase_sim.called
    assert pg_flat_mock.SLC_diff_intf.called
    assert pg_flat_mock.base_perp.call_count == 1


@pytest.fixture
def pg_filt_mock():
    """Create basic mock of the py_gamma module for the INT processing step."""
    pgm = mock.Mock()
    pgm.adf.return_value = 0
    return pgm


def test_calc_filt(monkeypatch, pg_filt_mock, pc_mock, ic_mock):
    monkeypatch.setattr(process_ifg, "pg", pg_filt_mock)

    ic_mock.ifg_flat = mock.Mock()
    ic_mock.ifg_flat.exists.return_value = True

    assert pg_filt_mock.adf.called is False
    process_ifg.calc_filt(pc_mock, ic_mock, ifg_width=230)
    assert pg_filt_mock.adf.called


def test_calc_filt_no_flat_file(monkeypatch, pg_filt_mock, pc_mock, ic_mock):
    monkeypatch.setattr(process_ifg, "pg", pg_filt_mock)

    ic_mock.ifg_flat = mock.Mock()
    ic_mock.ifg_flat.exists.return_value = False

    with pytest.raises(ProcessIfgException):
        process_ifg.calc_filt(pc_mock, ic_mock, ifg_width=180)


@pytest.fixture
def pg_unw_mock():
    pgm = mock.Mock()
    pgm.rascc_mask.return_value = 0
    pgm.rascc_mask_thinning.return_value = 0
    pgm.mcf.return_value = 0
    pgm.interp_ad.return_value = 0
    pgm.unw_model.return_value = 0
    pgm.mask_data.return_value = 0
    return pgm


def test_calc_unw(monkeypatch, pg_unw_mock, pc_mock, ic_mock):
    # NB: (m)looks will always be 2 for Sentinel-1 ARD product generation
    monkeypatch.setattr(process_ifg, "pg", pg_unw_mock)

    # ignore the thinning step as it will be tested separately
    m_thin = mock.Mock()
    monkeypatch.setattr(process_ifg, "calc_unw_thinning", m_thin)

    pc_mock.multi_look = 2
    pc_mock.ifg_coherence_threshold = 1  # fake value
    pc_mock.ifg_unw_mask = "no"
    fake_ifg_width = 13

    assert pg_unw_mock.rascc_mask.called is False
    assert m_thin.called is False
    assert pg_unw_mock.mask_data.called is False

    process_ifg.calc_unw(pc_mock, ic_mock, fake_ifg_width, clean_up=False)

    assert pg_unw_mock.rascc_mask.called
    assert m_thin.called
    assert pg_unw_mock.mask_data.called is False


def test_calc_unw_mlooks_over_threshold_not_implemented(
    monkeypatch, pg_unw_mock, pc_mock, ic_mock
):
    monkeypatch.setattr(process_ifg, "pg", pg_unw_mock)
    pc_mock.multi_look = 5

    with pytest.raises(NotImplementedError):
        process_ifg.calc_unw(pc_mock, ic_mock, ifg_width=15, clean_up=False)


def test_calc_unw_thinning(monkeypatch, pg_unw_mock, pc_mock, ic_mock):
    monkeypatch.setattr(process_ifg, "pg", pg_unw_mock)

    pc_mock.ifg_coherence_threshold = 2.5  # fake value
    ifg_width = 37  # fake value

    assert pg_unw_mock.rascc_mask_thinning.called is False
    assert pg_unw_mock.mcf.called is False
    assert pg_unw_mock.interp_ad.called is False
    assert pg_unw_mock.unw_model.called is False

    process_ifg.calc_unw_thinning(pc_mock, ic_mock, ifg_width, clean_up=False)

    assert pg_unw_mock.rascc_mask_thinning.called
    assert pg_unw_mock.mcf.called
    assert pg_unw_mock.interp_ad.called
    assert pg_unw_mock.unw_model.called


@pytest.fixture
def pg_geocode_mock():
    """Basic mock for pygamma calls in GEOCODE"""
    pgm = mock.Mock()
    ret = (0, ["cout for pg_geocode_mock"], ["cerr for pg_geocode_mock"])

    pgm.geocode_back.return_value = ret
    pgm.mask_data.return_value = ret
    pgm.convert.return_value = ret
    pgm.kml_map.return_value = ret
    pgm.cpx_to_real.return_value = ret
    pgm.rascc.return_value = ret
    pgm.ras2ras.return_value = ret
    pgm.rasrmg.return_value = ret
    pgm.data2geotiff.return_value = ret
    return pgm


# TODO: can fixtures call other fixtures to get their setup? (e.g. mock pg inside another fixture?)
def test_geocode_unwrapped_ifg(monkeypatch, ic_mock, dc_mock, pg_geocode_mock):
    monkeypatch.setattr(process_ifg, "pg", pg_geocode_mock)

    # patch at the subprocess level for testing this part of convert() in geocode step
    m_subprocess = mock.Mock(spec=subprocess)
    m_subprocess.run.return_value = 0
    monkeypatch.setattr(process_ifg, "subprocess", m_subprocess)

    m_remove = mock.Mock()
    monkeypatch.setattr(process_ifg, "remove_files", m_remove)

    assert pg_geocode_mock.geocode_back.called is False
    assert pg_geocode_mock.mask_data.called is False
    assert pg_geocode_mock.rasrmg.called is False
    assert pg_geocode_mock.kml_map.called is False

    assert m_subprocess.run.called is False
    assert m_remove.called is False

    width_in, width_out = 5, 7  # fake values
    process_ifg.geocode_unwrapped_ifg(ic_mock, dc_mock, width_in, width_out)

    assert pg_geocode_mock.geocode_back.called
    assert pg_geocode_mock.mask_data.called
    assert pg_geocode_mock.rasrmg.called
    assert pg_geocode_mock.kml_map.called

    assert m_subprocess.run.called
    assert m_remove.called


def test_geocode_flattened_ifg(monkeypatch, ic_mock, dc_mock, pg_geocode_mock):
    monkeypatch.setattr(process_ifg, "pg", pg_geocode_mock)

    # patch convert function for testing this part of geocode step
    m_convert = mock.Mock(spec=process_ifg.convert)
    monkeypatch.setattr(process_ifg, "convert", m_convert)

    m_remove = mock.Mock()
    monkeypatch.setattr(process_ifg, "remove_files", m_remove)

    assert pg_geocode_mock.cpx_to_real.called is False
    assert pg_geocode_mock.geocode_back.called is False
    assert pg_geocode_mock.mask_data.called is False
    assert pg_geocode_mock.rasrmg.called is False
    assert pg_geocode_mock.kml_map.called is False
    assert m_convert.called is False
    assert m_remove.called is False

    width_in, width_out = 9, 13  # fake values
    process_ifg.geocode_flattened_ifg(ic_mock, dc_mock, width_in, width_out)

    assert pg_geocode_mock.cpx_to_real.called
    assert pg_geocode_mock.geocode_back.called
    assert pg_geocode_mock.mask_data.called
    assert pg_geocode_mock.rasrmg.called
    assert pg_geocode_mock.kml_map.called
    assert m_convert.called
    assert m_remove.called


def test_geocode_filtered_ifg(monkeypatch, ic_mock, dc_mock, pg_geocode_mock):
    monkeypatch.setattr(process_ifg, "pg", pg_geocode_mock)

    # patch convert function for testing this part of geocode step
    m_convert = mock.Mock(spec=process_ifg.convert)
    monkeypatch.setattr(process_ifg, "convert", m_convert)

    m_remove = mock.Mock()
    monkeypatch.setattr(process_ifg, "remove_files", m_remove)

    assert pg_geocode_mock.cpx_to_real.called is False
    assert pg_geocode_mock.geocode_back.called is False
    assert pg_geocode_mock.mask_data.called is False
    assert pg_geocode_mock.rasrmg.called is False
    assert pg_geocode_mock.kml_map.called is False
    assert m_convert.called is False
    assert m_remove.called is False

    width_in, width_out = 15, 19  # fake values
    process_ifg.geocode_filtered_ifg(ic_mock, dc_mock, width_in, width_out)

    assert pg_geocode_mock.cpx_to_real.called
    assert pg_geocode_mock.geocode_back.called
    assert pg_geocode_mock.mask_data.called
    assert pg_geocode_mock.rasrmg.called
    assert pg_geocode_mock.kml_map.called
    assert m_convert.called
    assert m_remove.called


def test_geocode_flat_coherence_file(monkeypatch, ic_mock, dc_mock, pg_geocode_mock):
    monkeypatch.setattr(process_ifg, "pg", pg_geocode_mock)

    # patch convert function for testing this part of geocode step
    m_convert = mock.Mock(spec=process_ifg.convert)
    monkeypatch.setattr(process_ifg, "convert", m_convert)

    m_remove = mock.Mock()
    monkeypatch.setattr(process_ifg, "remove_files", m_remove)

    assert pg_geocode_mock.geocode_back.called is False
    assert pg_geocode_mock.rascc.called is False
    assert pg_geocode_mock.ras2ras.called is False
    assert pg_geocode_mock.kml_map.called is False
    assert m_convert.called is False

    width_in, width_out = 33, 37  # fake values
    process_ifg.geocode_flat_coherence_file(ic_mock, dc_mock, width_in, width_out)

    assert pg_geocode_mock.geocode_back.called
    assert pg_geocode_mock.rascc.called
    assert pg_geocode_mock.ras2ras.called
    assert pg_geocode_mock.kml_map.called
    assert m_convert.called


def test_geocode_filtered_coherence_file(monkeypatch, ic_mock, dc_mock, pg_geocode_mock):
    monkeypatch.setattr(process_ifg, "pg", pg_geocode_mock)

    m_convert = mock.Mock(spec=process_ifg.convert)
    monkeypatch.setattr(process_ifg, "convert", m_convert)
    m_remove = mock.Mock()
    monkeypatch.setattr(process_ifg, "remove_files", m_remove)

    assert pg_geocode_mock.geocode_back.called is False
    assert pg_geocode_mock.rascc.called is False
    assert pg_geocode_mock.ras2ras.called is False
    assert pg_geocode_mock.kml_map.called is False
    assert m_convert.called is False

    width_in, width_out = 43, 31  # fake values
    process_ifg.geocode_filtered_coherence_file(ic_mock, dc_mock, width_in, width_out)

    assert pg_geocode_mock.geocode_back.called
    assert pg_geocode_mock.rascc.called
    assert pg_geocode_mock.ras2ras.called
    assert pg_geocode_mock.kml_map.called
    assert m_convert.called


def test_do_geocode(monkeypatch, pc_mock, ic_mock, dc_mock, pg_geocode_mock):
    """Test the full geocode step"""
    monkeypatch.setattr(process_ifg, "pg", pg_geocode_mock)

    pc_mock.ifg_geotiff.lower.return_value = "yes"

    m_geocode_unwrapped_ifg = mock.Mock()
    m_geocode_flattened_ifg = mock.Mock()
    m_geocode_filtered_ifg = mock.Mock()
    m_geocode_flat_coherence_file = mock.Mock()
    m_geocode_filtered_coherence_file = mock.Mock()
    m_remove_files = mock.Mock()

    monkeypatch.setattr(process_ifg, "geocode_unwrapped_ifg", m_geocode_unwrapped_ifg)
    monkeypatch.setattr(process_ifg, "geocode_flattened_ifg", m_geocode_flattened_ifg)
    monkeypatch.setattr(process_ifg, "geocode_filtered_ifg", m_geocode_filtered_ifg)
    monkeypatch.setattr(
        process_ifg, "geocode_flat_coherence_file", m_geocode_flat_coherence_file
    )
    monkeypatch.setattr(
        process_ifg, "geocode_filtered_coherence_file", m_geocode_filtered_coherence_file
    )
    monkeypatch.setattr(process_ifg, "remove_files", m_remove_files)

    width_in, width_out = 58, 67  # fake values
    process_ifg.do_geocode(pc_mock, ic_mock, dc_mock, width_in, width_out)

    assert m_geocode_unwrapped_ifg.called
    assert m_geocode_flattened_ifg.called
    assert m_geocode_filtered_ifg.called
    assert m_geocode_flat_coherence_file.called
    assert m_geocode_filtered_ifg.called

    assert pg_geocode_mock.data2geotiff.call_count == 5
    assert m_remove_files.call_count == len(const.TEMP_FILE_GLOBS)
