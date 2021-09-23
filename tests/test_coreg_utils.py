from pathlib import Path
from insar import coreg_utils

from unittest import mock

import pytest


def test_read_land_center_coords():
    shape_path = Path("tests/data/T147D_F28S_S1A.shp")
    coords = coreg_utils.read_land_center_coords(shape_path)

    ncoord, ecoord = coords
    assert ncoord == -28.129263619  # from manual extraction using ogr2ogr
    assert ecoord == 152.169759163


def test_read_land_center_coords_missing_file():
    shape_path = Path("/tmp/this-is-fake/dummy.shp")

    with pytest.raises(FileNotFoundError):
        coreg_utils.read_land_center_coords(shape_path)


def test_missing_land_centre_attrs():
    # TODO: log if there no land centre in the shapefile?
    with mock.patch("pathlib.Path.exists") as exists_mock:
        exists_mock.return_value = True

        with mock.patch("geopandas.GeoDataFrame.from_file") as from_file_mock:
            for attr in ("land_cen_l", "land_cen_1"):
                m = object()
                assert not hasattr(m, attr)  # need obj without land centre attrs
                from_file_mock.return_value = m

                path = Path("fake/file/path")
                assert coreg_utils.read_land_center_coords(path) is None


def test_zero_land_centre_attrs():
    with mock.patch("pathlib.Path.exists") as exists_mock:
        exists_mock.return_value = True

        with mock.patch("geopandas.GeoDataFrame.from_file") as from_file_mock:
            # test multiple forms of having a "0" coordinate in DBF data
            for values in (("0", "149.1"), ("35.4", "0"), ("0", "0")):
                m_dbf = mock.NonCallableMock()
                m_dbf.land_cen_l = [values[0]]
                m_dbf.land_cen_1 = [values[1]]

                from_file_mock.return_value = m_dbf

                path = Path("fake/file/path")
                assert coreg_utils.read_land_center_coords(path) is None
