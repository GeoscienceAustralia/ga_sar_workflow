import shutil

import pytest
from insar.process_tsx_slc import process_tsx_slc, ProcessSlcException

from tests.fixtures import pgp, pgmock, temp_out_dir, test_data_dir, tsx_test_data
from tests.fixtures import TSX_TEST_DATA_DATES, TSX_TEST_DATA_SUBDIRS


DATE0 = TSX_TEST_DATA_DATES[0]


def test_tsx_slc_processing(pgp, pgmock, temp_out_dir, tsx_test_data):
    # Run SLC processing in a temp dir
    assert tsx_test_data[0].is_dir()
    process_tsx_slc(tsx_test_data[0], temp_out_dir)

    # Ensure the output is created and no errors occurred
    assert pgp.error_count == 0
    assert len(pgp.call_sequence) >= 2

    out_path = (temp_out_dir / DATE0).with_suffix(".slc")
    out_par_path = (temp_out_dir / DATE0).with_suffix(".slc.par")
    assert out_par_path.exists()
    assert out_path.exists()


def test_tsx_slc_fails_with_missing_input(pgp, pgmock, temp_out_dir, tsx_test_data):
    # try running SLC processing from a src dir that does not exist
    data_path = temp_out_dir / "21230102"

    with pytest.raises(RuntimeError):
        process_tsx_slc(data_path, temp_out_dir)

    # Ensure not a single GAMMA call occurred & no output exists
    assert len(pgp.call_sequence) == 0
    out_slc_paths = list(temp_out_dir.glob("*.slc"))
    assert out_slc_paths == []


def test_tsx_slc_fails_with_incomplete_data(pgp, pgmock, temp_out_dir, tsx_test_data):
    # copy test data & avoid modifying the source data
    data_copy = temp_out_dir / "tsx_data_copy"
    shutil.copytree(tsx_test_data[0], data_copy)

    # "Delete" important data from the set so it's incomplete
    scene_date = tsx_test_data[0].name
    annotation_xml = (data_copy / scene_date / TSX_TEST_DATA_SUBDIRS[0] / TSX_TEST_DATA_SUBDIRS[0]).with_suffix(".xml")
    shutil.move(annotation_xml, annotation_xml.with_suffix(".bak"))

    # Run the SLC processing in a temp dest dir
    with pytest.raises(ProcessSlcException):
        process_tsx_slc(data_copy, temp_out_dir)

    # Ensure not a single GAMMA call occurred & no output exists
    assert len(pgp.call_sequence) == 0
    paths = list(temp_out_dir.glob("*"))
    assert paths == [data_copy]  # nothing produced, only the data dir should exist there


def test_tsx_slc_fails_with_incomplete_missing_dir(pgp, pgmock, temp_out_dir):
    # try instance where the big ugly dir name is missing
    data = temp_out_dir / "20001234"
    data.mkdir(exist_ok=False)  # error, root dir doesn't contain the ugly subdir

    # Run the SLC processing in a temp dest dir
    with pytest.raises(ProcessSlcException):
        process_tsx_slc(data, temp_out_dir)

    assert list(temp_out_dir.glob("*")) == [data]
    assert len(pgp.call_sequence) == 0
