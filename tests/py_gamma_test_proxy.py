from pathlib import Path
from typing import Sequence, NamedTuple, Dict, Union

PyGammaCall = NamedTuple["PyGammaCall", [("module", str), ("program", str), ("parameters", Dict[str, object])]]


class PyGammaTestProxy(object):
    call_sequence: Sequence[PyGammaCall]
    call_count: Dict[str, int]

    def __init__(self):
        self.reset_proxy()

    def reset_proxy(self):
        self.call_sequence = []
        self.call_count = {}

    def _validate(self, condition, result):
        stat, stdout, stderr = result

        # TODO: error stats?
        stat = stat if condition else -1
        # TODO: stderr?

        return stat, stdout, stderr

    def gc_map_fine(self, gc_in: str, width, DIFF_par, gc_out: str, ref_flg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map_fine", supplied_args))

        if "gc_map_fine" in self.call_count:
            self.call_count["gc_map_fine"] += 1
        else:
            self.call_count["gc_map_fine"] = 1

        result = self._validate(Path(gc_in).exists(), result)
        Path(gc_out).touch()
        valid_values = [0, 1]
        result = self._validate(ref_flg in valid_values, result)
        return result

    def diff_ls_fit(self, unw_1: str, unw_2: str, DIFF_par: str, nr, naz, mask, plot_data: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "diff_ls_fit", supplied_args))

        if "diff_ls_fit" in self.call_count:
            self.call_count["diff_ls_fit"] += 1
        else:
            self.call_count["diff_ls_fit"] = 1

        result = self._validate(Path(unw_1).exists(), result)
        result = self._validate(Path(unw_2).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(plot_data).touch()
        return result

    def WSS_mosaic(self, WSS_tab: str, MLI_par: str, WSS_data: str, type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "WSS_mosaic", supplied_args))

        if "WSS_mosaic" in self.call_count:
            self.call_count["WSS_mosaic"] += 1
        else:
            self.call_count["WSS_mosaic"] = 1

        result = self._validate(Path(WSS_tab).exists(), result)
        Path(MLI_par).touch()
        Path(WSS_data).touch()
        valid_values = [0, 1]
        result = self._validate(type in valid_values, result)
        return result

    def dispmap_vec_offset(self, DEM_par: str, DEM: str, dispmap_r: str, dispmap_az: str, lv_theta: str, lv_phi: str, dv_norm: str, dv_theta: str, dv_phi: str, dv_x: str, dv_y: str, dv_z: str, mask_angle, mode, ax_north, ax_east):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_vec_offset", supplied_args))

        if "dispmap_vec_offset" in self.call_count:
            self.call_count["dispmap_vec_offset"] += 1
        else:
            self.call_count["dispmap_vec_offset"] = 1

        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        result = self._validate(Path(dispmap_r).exists(), result)
        result = self._validate(Path(dispmap_az).exists(), result)
        result = self._validate(Path(lv_theta).exists(), result)
        result = self._validate(Path(lv_phi).exists(), result)
        Path(dv_norm).touch()
        Path(dv_theta).touch()
        Path(dv_phi).touch()
        Path(dv_x).touch()
        Path(dv_y).touch()
        Path(dv_z).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def create_diff_par(self, PAR_1: str, PAR_2: str, DIFF_par, PAR_type, iflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "create_diff_par", supplied_args))

        if "create_diff_par" in self.call_count:
            self.call_count["create_diff_par"] += 1
        else:
            self.call_count["create_diff_par"] = 1

        result = self._validate(Path(PAR_1).exists(), result)
        result = self._validate(Path(PAR_2).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(PAR_type in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(iflg in valid_values, result)
        return result

    def gc_map_grd(self, GRD_par: str, DEM_par: str, DEM: str, DEM_seg_par, DEM_seg: str, lookup_table: str, lat_ovr, lon_ovr, sim_sar: str, u: str, v: str, inc: str, psi: str, pix: str, ls_map: str, frame, ls_mode, r_ovr):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map_grd", supplied_args))

        if "gc_map_grd" in self.call_count:
            self.call_count["gc_map_grd"] += 1
        else:
            self.call_count["gc_map_grd"] = 1

        result = self._validate(Path(GRD_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        Path(DEM_seg).touch()
        Path(lookup_table).touch()
        Path(sim_sar).touch()
        Path(u).touch()
        Path(v).touch()
        Path(inc).touch()
        Path(psi).touch()
        Path(pix).touch()
        Path(ls_map).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(ls_mode in valid_values, result)
        return result

    def multi_mosaic(self, data_tab: str, data_out: str, DEM_par_out: str, mode, format_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "multi_mosaic", supplied_args))

        if "multi_mosaic" in self.call_count:
            self.call_count["multi_mosaic"] += 1
        else:
            self.call_count["multi_mosaic"] = 1

        result = self._validate(Path(data_tab).exists(), result)
        Path(data_out).touch()
        Path(DEM_par_out).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(format_flag in valid_values, result)
        return result

    def gc_map1(self, MLI_par: str, OFF_par: str, DEM_par: str, DEM: str, DEM_seg_par, DEM_seg: str, lookup_table: str, lat_ovr, lon_ovr, sim_sar: str, u: str, v: str, inc: str, psi: str, pix: str, ls_map: str, frame, ls_mode, r_ovr):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map1", supplied_args))

        if "gc_map1" in self.call_count:
            self.call_count["gc_map1"] += 1
        else:
            self.call_count["gc_map1"] = 1

        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        Path(DEM_seg).touch()
        Path(lookup_table).touch()
        Path(sim_sar).touch()
        Path(u).touch()
        Path(v).touch()
        Path(inc).touch()
        Path(psi).touch()
        Path(pix).touch()
        Path(ls_map).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(ls_mode in valid_values, result)
        return result

    def SLC_diff_intf(self, SLC_1: str, SLC_2R: str, SLC1_par: str, SLC2R_par: str, OFF_par: str, sim_unw: str, diff_int: str, rlks, azlks, sps_flg, azf_flg, rbw_min, rp1_flg, rp2_flg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "SLC_diff_intf", supplied_args))

        if "SLC_diff_intf" in self.call_count:
            self.call_count["SLC_diff_intf"] += 1
        else:
            self.call_count["SLC_diff_intf"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2R).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2R_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(sim_unw).exists(), result)
        Path(diff_int).touch()
        valid_values = [1, 0]
        result = self._validate(sps_flg in valid_values, result)
        valid_values = [1, 0]
        result = self._validate(azf_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(rp1_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(rp2_flg in valid_values, result)
        return result

    def pixel_area(self, MLI_par: str, DEM_par: str, DEM: str, lookup_table: str, ls_map: str, inc_map: str, pix_sigma0: str, pix_gamma0: str, nstep, area_fact, sigma0_ratio: str, gamma0_ratio: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "pixel_area", supplied_args))

        if "pixel_area" in self.call_count:
            self.call_count["pixel_area"] += 1
        else:
            self.call_count["pixel_area"] = 1

        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        result = self._validate(Path(lookup_table).exists(), result)
        result = self._validate(Path(ls_map).exists(), result)
        result = self._validate(Path(inc_map).exists(), result)
        Path(pix_sigma0).touch()
        Path(pix_gamma0).touch()
        Path(sigma0_ratio).touch()
        Path(gamma0_ratio).touch()
        return result

    def extract_gcp(self, DEM_rdc: str, OFF_par: str, GCP: str, nr, naz, mask: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "extract_gcp", supplied_args))

        if "extract_gcp" in self.call_count:
            self.call_count["extract_gcp"] += 1
        else:
            self.call_count["extract_gcp"] = 1

        result = self._validate(Path(DEM_rdc).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(GCP).touch()
        result = self._validate(Path(mask).exists(), result)
        return result

    def gc_map2(self, MLI_par: str, DEM_par: str, DEM: str, DEM_seg_par: str, DEM_seg: str, lookup_table: str, lat_ovr, lon_ovr, ls_map: str, ls_map_rdc: str, inc: str, res: str, offnadir: str, sim_sar: str, u: str, v: str, psi: str, pix: str, r_ovr, az_dec, mask, frame, ls_scaling, DIFF_par: str, ref_flg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map2", supplied_args))

        if "gc_map2" in self.call_count:
            self.call_count["gc_map2"] += 1
        else:
            self.call_count["gc_map2"] = 1

        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        Path(DEM_seg_par).touch()
        Path(DEM_seg).touch()
        Path(lookup_table).touch()
        Path(ls_map).touch()
        Path(ls_map_rdc).touch()
        Path(inc).touch()
        Path(res).touch()
        Path(offnadir).touch()
        Path(sim_sar).touch()
        Path(u).touch()
        Path(v).touch()
        Path(psi).touch()
        Path(pix).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(mask in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(ls_scaling in valid_values, result)
        result = self._validate(Path(DIFF_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(ref_flg in valid_values, result)
        return result

    def dem_trans(self, DEM1_par: str, DEM1: str, DEM2_par, DEM2: str, lat_ovr, lon_ovr, datum_shift, bflg, lookup_table: str, interp_mode, order):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_trans", supplied_args))

        if "dem_trans" in self.call_count:
            self.call_count["dem_trans"] += 1
        else:
            self.call_count["dem_trans"] = 1

        result = self._validate(Path(DEM1_par).exists(), result)
        result = self._validate(Path(DEM1).exists(), result)
        Path(DEM2).touch()
        valid_values = [0, 1]
        result = self._validate(datum_shift in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        Path(lookup_table).touch()
        valid_values = [0, 1, 2]
        result = self._validate(interp_mode in valid_values, result)
        return result

    def quad_sub(self, int_1: str, DIFF_par: str, int_2: str, int_type, mode):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "quad_sub", supplied_args))

        if "quad_sub" in self.call_count:
            self.call_count["quad_sub"] += 1
        else:
            self.call_count["quad_sub"] = 1

        result = self._validate(Path(int_1).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(int_2).touch()
        valid_values = [0, 1]
        result = self._validate(int_type in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def phase_sum(self, im_list: str, width, sum: str, start, nlines, pixav_x, pixav_y, zflag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "phase_sum", supplied_args))

        if "phase_sum" in self.call_count:
            self.call_count["phase_sum"] += 1
        else:
            self.call_count["phase_sum"] = 1

        result = self._validate(Path(im_list).exists(), result)
        Path(sum).touch()
        valid_values = [0, 1]
        result = self._validate(zflag in valid_values, result)
        return result

    def base_add(self, base_1: str, base_2: str, base_out: str, mode):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "base_add", supplied_args))

        if "base_add" in self.call_count:
            self.call_count["base_add"] += 1
        else:
            self.call_count["base_add"] = 1

        result = self._validate(Path(base_1).exists(), result)
        result = self._validate(Path(base_2).exists(), result)
        Path(base_out).touch()
        return result

    def gec_map_grd(self, GRD_par: str, DEM_par: str, href: str, DEM_seg_par, lookup_table: str, lat_ovr, lon_ovr):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gec_map_grd", supplied_args))

        if "gec_map_grd" in self.call_count:
            self.call_count["gec_map_grd"] += 1
        else:
            self.call_count["gec_map_grd"] = 1

        result = self._validate(Path(GRD_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(href).exists(), result)
        Path(lookup_table).touch()
        return result

    def sarpix_coord_list(self, SLC_par: str, OFF_par: str, DEM_par: str, SAR_coord: str, MAP_coord: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "sarpix_coord_list", supplied_args))

        if "sarpix_coord_list" in self.call_count:
            self.call_count["sarpix_coord_list"] += 1
        else:
            self.call_count["sarpix_coord_list"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(SAR_coord).exists(), result)
        Path(MAP_coord).touch()
        return result

    def dispmap_vec(self, DEM_par: str, dispmap: str, lv_theta: str, lv_phi: str, fv_theta: str, fv_phi: str, dv_norm: str, dv_theta: str, dv_phi: str, dv_x: str, dv_y: str, dv_z: str, mask_angle):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_vec", supplied_args))

        if "dispmap_vec" in self.call_count:
            self.call_count["dispmap_vec"] += 1
        else:
            self.call_count["dispmap_vec"] = 1

        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(dispmap).exists(), result)
        result = self._validate(Path(lv_theta).exists(), result)
        result = self._validate(Path(lv_phi).exists(), result)
        result = self._validate(Path(fv_theta).exists(), result)
        result = self._validate(Path(fv_phi).exists(), result)
        Path(dv_norm).touch()
        Path(dv_theta).touch()
        Path(dv_phi).touch()
        Path(dv_x).touch()
        Path(dv_y).touch()
        Path(dv_z).touch()
        return result

    def scale_base(self, unw_2: str, scaled_unw_2: str, baseline_1: str, SLC1_par_1: str, OFF_par_1: str, baseline_2: str, SLC1_par_2: str, OFF_par_2: str, int_type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "scale_base", supplied_args))

        if "scale_base" in self.call_count:
            self.call_count["scale_base"] += 1
        else:
            self.call_count["scale_base"] = 1

        result = self._validate(Path(unw_2).exists(), result)
        Path(scaled_unw_2).touch()
        result = self._validate(Path(baseline_1).exists(), result)
        result = self._validate(Path(SLC1_par_1).exists(), result)
        result = self._validate(Path(OFF_par_1).exists(), result)
        result = self._validate(Path(baseline_2).exists(), result)
        result = self._validate(Path(SLC1_par_2).exists(), result)
        result = self._validate(Path(OFF_par_2).exists(), result)
        return result

    def par_EORC_PALSAR_geo(self, CEOS_leader: str, MLI_par: str, DEM_par: str, CEOS_data: str, MLI: str, cal):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_EORC_PALSAR_geo", supplied_args))

        if "par_EORC_PALSAR_geo" in self.call_count:
            self.call_count["par_EORC_PALSAR_geo"] += 1
        else:
            self.call_count["par_EORC_PALSAR_geo"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        Path(MLI_par).touch()
        Path(DEM_par).touch()
        result = self._validate(Path(CEOS_data).exists(), result)
        Path(MLI).touch()
        return result

    def SLC_intf_geo(self, SLC_1: str, SLC_2: str, DEM_par: str, interf: str, DEM_par2: str, e_lks, n_lks, MLI_1: str, MLI_2: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "SLC_intf_geo", supplied_args))

        if "SLC_intf_geo" in self.call_count:
            self.call_count["SLC_intf_geo"] += 1
        else:
            self.call_count["SLC_intf_geo"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        Path(interf).touch()
        Path(DEM_par2).touch()
        Path(MLI_1).touch()
        Path(MLI_2).touch()
        return result

    def map_trans(self, DEM1_par: str, data1: str, DEM2_par, data2: str, lat_ovr, lon_ovr, interp_mode, dtype, bflg, lookup_table: str, order):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "map_trans", supplied_args))

        if "map_trans" in self.call_count:
            self.call_count["map_trans"] += 1
        else:
            self.call_count["map_trans"] = 1

        result = self._validate(Path(DEM1_par).exists(), result)
        result = self._validate(Path(data1).exists(), result)
        Path(data2).touch()
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        Path(lookup_table).touch()
        return result

    def gc_map_inversion(self, gc_map, width_in, gc_map_out, width_out, nlines_out, interp_mode, n_ovr, rad_max, nintr):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map_inversion", supplied_args))

        if "gc_map_inversion" in self.call_count:
            self.call_count["gc_map_inversion"] += 1
        else:
            self.call_count["gc_map_inversion"] = 1

        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(interp_mode in valid_values, result)
        return result

    def par_JERS_geo(self, CEOS_leader: str, CEOS_data: str, MLI_par: str, DEM_par: str, GEO: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_JERS_geo", supplied_args))

        if "par_JERS_geo" in self.call_count:
            self.call_count["par_JERS_geo"] += 1
        else:
            self.call_count["par_JERS_geo"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        result = self._validate(Path(CEOS_data).exists(), result)
        Path(MLI_par).touch()
        Path(DEM_par).touch()
        Path(GEO).touch()
        return result

    def interp_data(self, data2: str, DIFF_par: str, data2_out: str, interp_mode, dtype, order):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "interp_data", supplied_args))

        if "interp_data" in self.call_count:
            self.call_count["interp_data"] += 1
        else:
            self.call_count["interp_data"] = 1

        result = self._validate(Path(data2).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(data2_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6]
        result = self._validate(dtype in valid_values, result)
        return result

    def rdc_trans(self, MLI1_par: str, DEM_RDC: str, MLI2_par: str, lt: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "rdc_trans", supplied_args))

        if "rdc_trans" in self.call_count:
            self.call_count["rdc_trans"] += 1
        else:
            self.call_count["rdc_trans"] = 1

        result = self._validate(Path(MLI1_par).exists(), result)
        result = self._validate(Path(DEM_RDC).exists(), result)
        result = self._validate(Path(MLI2_par).exists(), result)
        Path(lt).touch()
        return result

    def SLC_interp_lt_ScanSAR(self, SLC2_tab: str, SLC2_par: str, SLC1_tab: str, SLC1_par: str, lookup_table: str, MLI1_par: str, MLI2_par: str, OFF_par: str, SLC2R_tab, SLC_2R: str, SLC2R_par: str, mode, order, SLC2R_dir):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "SLC_interp_lt_ScanSAR", supplied_args))

        if "SLC_interp_lt_ScanSAR" in self.call_count:
            self.call_count["SLC_interp_lt_ScanSAR"] += 1
        else:
            self.call_count["SLC_interp_lt_ScanSAR"] = 1

        result = self._validate(Path(SLC2_tab).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(SLC1_tab).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(lookup_table).exists(), result)
        result = self._validate(Path(MLI1_par).exists(), result)
        result = self._validate(Path(MLI2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(SLC_2R).touch()
        Path(SLC2R_par).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def par_CS_geo(self, HDF5: str, trunk: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_CS_geo", supplied_args))

        if "par_CS_geo" in self.call_count:
            self.call_count["par_CS_geo"] += 1
        else:
            self.call_count["par_CS_geo"] = 1

        result = self._validate(Path(HDF5).exists(), result)
        result = self._validate(Path(trunk).exists(), result)
        return result

    def par_RCM_geo(self, RCM_dir: str, polarization, MLI_par: str, DEM_par: str, GEO: str, dtype, ps):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_RCM_geo", supplied_args))

        if "par_RCM_geo" in self.call_count:
            self.call_count["par_RCM_geo"] += 1
        else:
            self.call_count["par_RCM_geo"] = 1

        result = self._validate(Path(RCM_dir).exists(), result)
        Path(MLI_par).touch()
        Path(DEM_par).touch()
        Path(GEO).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def offset_pwr_trackingm2(self, MLI_1: str, MLI_2: str, DIFF_par: str, offs: str, ccp: str, DIFF_par2: str, offs2: str, rwin, azwin, offsets: str, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, bw_frac, pflag, pltflg, ccs: str, std_mean):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_pwr_trackingm2", supplied_args))

        if "offset_pwr_trackingm2" in self.call_count:
            self.call_count["offset_pwr_trackingm2"] += 1
        else:
            self.call_count["offset_pwr_trackingm2"] = 1

        result = self._validate(Path(MLI_1).exists(), result)
        result = self._validate(Path(MLI_2).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(offs).touch()
        Path(ccp).touch()
        result = self._validate(Path(DIFF_par2).exists(), result)
        result = self._validate(Path(offs2).exists(), result)
        Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        Path(ccs).touch()
        return result

    def init_offset_orbitm(self, MLI1_par: str, MLI2_par: str, DIFF_par, rpos, azpos, cflag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "init_offset_orbitm", supplied_args))

        if "init_offset_orbitm" in self.call_count:
            self.call_count["init_offset_orbitm"] += 1
        else:
            self.call_count["init_offset_orbitm"] = 1

        result = self._validate(Path(MLI1_par).exists(), result)
        result = self._validate(Path(MLI2_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(cflag in valid_values, result)
        return result

    def quad_fit(self, unw: str, DIFF_par: str, dr, daz, mask, plot_data: str, model, pmodel: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "quad_fit", supplied_args))

        if "quad_fit" in self.call_count:
            self.call_count["quad_fit"] += 1
        else:
            self.call_count["quad_fit"] = 1

        result = self._validate(Path(unw).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(plot_data).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6, 7]
        result = self._validate(model in valid_values, result)
        Path(pmodel).touch()
        return result

    def stacking(self, DIFF_tab: str, width, ph_rate: str, sig_ph_rate: str, sig_ph: str, roff, loff, nr, nl, np_min, tscale):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "stacking", supplied_args))

        if "stacking" in self.call_count:
            self.call_count["stacking"] += 1
        else:
            self.call_count["stacking"] = 1

        result = self._validate(Path(DIFF_tab).exists(), result)
        Path(ph_rate).touch()
        Path(sig_ph_rate).touch()
        Path(sig_ph).touch()
        valid_values = [0, 1]
        result = self._validate(tscale in valid_values, result)
        return result

    def atm_mod2(self, diff_unw: str, hgt: str, MLI_par: str, model, dr, daz, mask: str, roff, loff, rpt: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "atm_mod2", supplied_args))

        if "atm_mod2" in self.call_count:
            self.call_count["atm_mod2"] += 1
        else:
            self.call_count["atm_mod2"] = 1

        result = self._validate(Path(diff_unw).exists(), result)
        result = self._validate(Path(hgt).exists(), result)
        result = self._validate(Path(MLI_par).exists(), result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(model in valid_values, result)
        result = self._validate(Path(mask).exists(), result)
        Path(rpt).touch()
        return result

    def par_UAVSAR_geo(self, ann: str, SLC_MLI_par: str, DEM_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_UAVSAR_geo", supplied_args))

        if "par_UAVSAR_geo" in self.call_count:
            self.call_count["par_UAVSAR_geo"] += 1
        else:
            self.call_count["par_UAVSAR_geo"] = 1

        result = self._validate(Path(ann).exists(), result)
        Path(SLC_MLI_par).touch()
        Path(DEM_par).touch()
        return result

    def dem_import(self, input_DEM: str, DEM: str, DEM_par, input_type, priority, geoid: str, geoid_par: str, geoid_type, latN_shift, lonE_shift, zflg, no_data):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_import", supplied_args))

        if "dem_import" in self.call_count:
            self.call_count["dem_import"] += 1
        else:
            self.call_count["dem_import"] = 1

        result = self._validate(Path(input_DEM).exists(), result)
        Path(DEM).touch()
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(input_type in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(priority in valid_values, result)
        result = self._validate(Path(geoid).exists(), result)
        result = self._validate(Path(geoid_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(geoid_type in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def SLC_interp_lt(self, SLC2: str, SLC1_par: str, SLC2_par: str, lookup_table: str, MLI1_par: str, MLI2_par: str, OFF_par: str, SLC_2R: str, SLC2R_par: str, blksz, mode, order):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "SLC_interp_lt", supplied_args))

        if "SLC_interp_lt" in self.call_count:
            self.call_count["SLC_interp_lt"] += 1
        else:
            self.call_count["SLC_interp_lt"] = 1

        result = self._validate(Path(SLC2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(lookup_table).exists(), result)
        result = self._validate(Path(MLI1_par).exists(), result)
        result = self._validate(Path(MLI2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(SLC_2R).touch()
        Path(SLC2R_par).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def dispmap(self, unw: str, hgt: str, MLI_par: str, OFF_par: str, disp_map: str, mode, sflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap", supplied_args))

        if "dispmap" in self.call_count:
            self.call_count["dispmap"] += 1
        else:
            self.call_count["dispmap"] = 1

        result = self._validate(Path(unw).exists(), result)
        result = self._validate(Path(hgt).exists(), result)
        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(disp_map).touch()
        return result

    def atm_mod(self, diff_unw: str, hgt: str, DIFF_par: str, model: str, dr, daz, mask: str, mode, roff, loff):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "atm_mod", supplied_args))

        if "atm_mod" in self.call_count:
            self.call_count["atm_mod"] += 1
        else:
            self.call_count["atm_mod"] = 1

        result = self._validate(Path(diff_unw).exists(), result)
        result = self._validate(Path(hgt).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(model).touch()
        result = self._validate(Path(mask).exists(), result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def dispmap_LOS(self, unw: str, width, freq, disp_map: str, sflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_LOS", supplied_args))

        if "dispmap_LOS" in self.call_count:
            self.call_count["dispmap_LOS"] += 1
        else:
            self.call_count["dispmap_LOS"] = 1

        result = self._validate(Path(unw).exists(), result)
        Path(disp_map).touch()
        return result

    def sub_phase(self, int_1: str, unw_2: str, DIFF_par: str, diff_int: str, dtype, mode):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "sub_phase", supplied_args))

        if "sub_phase" in self.call_count:
            self.call_count["sub_phase"] += 1
        else:
            self.call_count["sub_phase"] = 1

        result = self._validate(Path(int_1).exists(), result)
        result = self._validate(Path(unw_2).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(diff_int).touch()
        valid_values = [0, 1, 2]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def phase_sim(self, SLC1_par: str, OFF_par: str, baseline: str, hgt: str, sim_unw: str, ph_flag, bflag, definition: str, delta_t: str, int_mode: str, SLC2R_par: str, ph_mode):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "phase_sim", supplied_args))

        if "phase_sim" in self.call_count:
            self.call_count["phase_sim"] += 1
        else:
            self.call_count["phase_sim"] = 1

        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(baseline).exists(), result)
        result = self._validate(Path(hgt).exists(), result)
        Path(sim_unw).touch()
        valid_values = [0, 1]
        result = self._validate(ph_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(bflag in valid_values, result)
        result = self._validate(Path(definition).exists(), result)
        result = self._validate(Path(delta_t).exists(), result)
        result = self._validate(Path(int_mode).exists(), result)
        result = self._validate(Path(SLC2R_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(ph_mode in valid_values, result)
        return result

    def par_RISAT_geo(self, annotation_XML: str, GeoTIFF: str, polarization: str, DEM_par: str, MLI_par: str, MLI: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_RISAT_geo", supplied_args))

        if "par_RISAT_geo" in self.call_count:
            self.call_count["par_RISAT_geo"] += 1
        else:
            self.call_count["par_RISAT_geo"] = 1

        result = self._validate(Path(annotation_XML).exists(), result)
        result = self._validate(Path(GeoTIFF).exists(), result)
        result = self._validate(Path(polarization).exists(), result)
        Path(DEM_par).touch()
        Path(MLI_par).touch()
        Path(MLI).touch()
        return result

    def MLI_interp_lt(self, MLI_2: str, MLI1_par: str, MLI2_par: str, lookup_table: str, MLI3_par: str, MLI4_par: str, DIFF_par: str, MLI_2R: str, MLI2R_par: str, order):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "MLI_interp_lt", supplied_args))

        if "MLI_interp_lt" in self.call_count:
            self.call_count["MLI_interp_lt"] += 1
        else:
            self.call_count["MLI_interp_lt"] = 1

        result = self._validate(Path(MLI_2).exists(), result)
        result = self._validate(Path(MLI1_par).exists(), result)
        result = self._validate(Path(MLI2_par).exists(), result)
        result = self._validate(Path(lookup_table).exists(), result)
        result = self._validate(Path(MLI3_par).exists(), result)
        result = self._validate(Path(MLI4_par).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(MLI_2R).touch()
        Path(MLI2R_par).touch()
        return result

    def lk_vec_lt(self, MLI_par: str, DEM_par: str, DEM: str, lt: str, lv_theta: str, lv_phi: str, lv_ENU: str, azv_ENU: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "lk_vec_lt", supplied_args))

        if "lk_vec_lt" in self.call_count:
            self.call_count["lk_vec_lt"] += 1
        else:
            self.call_count["lk_vec_lt"] = 1

        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        result = self._validate(Path(lt).exists(), result)
        Path(lv_theta).touch()
        Path(lv_phi).touch()
        Path(lv_ENU).touch()
        Path(azv_ENU).touch()
        return result

    def coord_to_sarpix(self, SLC_MLI_par: str, OFF_par: str, DEM_par: str, north_lat: str, east_lon: str, hgt: str, DIFF_par):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "coord_to_sarpix", supplied_args))

        if "coord_to_sarpix" in self.call_count:
            self.call_count["coord_to_sarpix"] += 1
        else:
            self.call_count["coord_to_sarpix"] = 1

        result = self._validate(Path(SLC_MLI_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(north_lat).exists(), result)
        result = self._validate(Path(east_lon).exists(), result)
        result = self._validate(Path(hgt).exists(), result)
        return result

    def WSS_intf(self, SLC_1: str, SLC_2R: str, SLC1_par: str, SLC2R_par: str, OFF_par: str, interf: str, rlks, sps_flg, azf_flg, m_flg, boff, bstep, bmax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "WSS_intf", supplied_args))

        if "WSS_intf" in self.call_count:
            self.call_count["WSS_intf"] += 1
        else:
            self.call_count["WSS_intf"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2R).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2R_par).exists(), result)
        Path(OFF_par).touch()
        Path(interf).touch()
        valid_values = [1, 0]
        result = self._validate(sps_flg in valid_values, result)
        valid_values = [1, 0]
        result = self._validate(azf_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(m_flg in valid_values, result)
        return result

    def map_section(self, DEM_par: str, north1, east1, north2, east2, post_north, post_east, DEM_par2: str, lt: str, MLI_par1: str, MLI_par2: str, cflg, lt2: str, MLI_coord: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "map_section", supplied_args))

        if "map_section" in self.call_count:
            self.call_count["map_section"] += 1
        else:
            self.call_count["map_section"] = 1

        result = self._validate(Path(DEM_par).exists(), result)
        Path(DEM_par2).touch()
        result = self._validate(Path(lt).exists(), result)
        result = self._validate(Path(MLI_par1).exists(), result)
        result = self._validate(Path(MLI_par2).exists(), result)
        Path(lt2).touch()
        Path(MLI_coord).touch()
        return result

    def offset_list_fitm(self, cp_list: str, DIFF_par, DEM_par: str, lookup_table: str, lt_type, type1, type2, coffsets: str, poly_order, interact_flag, trans_list):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_list_fitm", supplied_args))

        if "offset_list_fitm" in self.call_count:
            self.call_count["offset_list_fitm"] += 1
        else:
            self.call_count["offset_list_fitm"] = 1

        result = self._validate(Path(cp_list).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(lookup_table).exists(), result)
        valid_values = [1, 2]
        result = self._validate(lt_type in valid_values, result)
        valid_values = [1, 2, 3]
        result = self._validate(type1 in valid_values, result)
        valid_values = [1, 2, 3]
        result = self._validate(type2 in valid_values, result)
        Path(coffsets).touch()
        return result

    def dem_RDC_list(self, DEM_par1: str, gc_map: str, MLI_par: str, mask: str, clist_RDC: str, clist_MAP: str, DEM_par2: str, s_north, s_east):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_RDC_list", supplied_args))

        if "dem_RDC_list" in self.call_count:
            self.call_count["dem_RDC_list"] += 1
        else:
            self.call_count["dem_RDC_list"] = 1

        result = self._validate(Path(DEM_par1).exists(), result)
        result = self._validate(Path(gc_map).exists(), result)
        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(mask).exists(), result)
        Path(clist_RDC).touch()
        Path(clist_MAP).touch()
        Path(DEM_par2).touch()
        return result

    def multi_look_geo(self, SLC: str, SLC_DEM_par: str, MLI: str, MLI_DEM_par: str, e_lks, n_lks, dtype, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "multi_look_geo", supplied_args))

        if "multi_look_geo" in self.call_count:
            self.call_count["multi_look_geo"] += 1
        else:
            self.call_count["multi_look_geo"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(SLC_DEM_par).exists(), result)
        Path(MLI).touch()
        Path(MLI_DEM_par).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def WSS_interp_lt(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, lookup_table: str, MLI1_par: str, MLI2_par: str, DIFF_par1: str, SLC_2R: str, SLC2R_par: str, DIFF_par2: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "WSS_interp_lt", supplied_args))

        if "WSS_interp_lt" in self.call_count:
            self.call_count["WSS_interp_lt"] += 1
        else:
            self.call_count["WSS_interp_lt"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(lookup_table).exists(), result)
        result = self._validate(Path(MLI1_par).exists(), result)
        result = self._validate(Path(MLI2_par).exists(), result)
        result = self._validate(Path(DIFF_par1).exists(), result)
        Path(SLC_2R).touch()
        Path(SLC2R_par).touch()
        Path(DIFF_par2).touch()
        return result

    def dispmap_vec2(self, DEM_par: str, DEM: str, dispmap1: str, lv1_theta: str, lv1_phi: str, dispmap2: str, lv2_theta: str, lv2_phi: str, dv_norm: str, dv_theta: str, dv_phi: str, dv_x: str, dv_y: str, dv_z: str, mask_angle, mode, ax_north, ax_east):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_vec2", supplied_args))

        if "dispmap_vec2" in self.call_count:
            self.call_count["dispmap_vec2"] += 1
        else:
            self.call_count["dispmap_vec2"] = 1

        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        result = self._validate(Path(dispmap1).exists(), result)
        result = self._validate(Path(lv1_theta).exists(), result)
        result = self._validate(Path(lv1_phi).exists(), result)
        result = self._validate(Path(dispmap2).exists(), result)
        result = self._validate(Path(lv2_theta).exists(), result)
        result = self._validate(Path(lv2_phi).exists(), result)
        Path(dv_norm).touch()
        Path(dv_theta).touch()
        Path(dv_phi).touch()
        Path(dv_x).touch()
        Path(dv_y).touch()
        Path(dv_z).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        return result

    def diff_ls_unw(self, int_1: str, unw_2: str, DIFF_par: str, diff_int: str, int_type, ph_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "diff_ls_unw", supplied_args))

        if "diff_ls_unw" in self.call_count:
            self.call_count["diff_ls_unw"] += 1
        else:
            self.call_count["diff_ls_unw"] = 1

        result = self._validate(Path(int_1).exists(), result)
        result = self._validate(Path(unw_2).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(diff_int).touch()
        return result

    def offset_pwr_list(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, clist_RDC: str, clist_MAP: str, offs: str, ccp: str, nx, ny, rwin, azwin, offsets: str, n_ovr, thres, bw_frac, deramp, int_filt, pflag, pltflg, ccs: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_pwr_list", supplied_args))

        if "offset_pwr_list" in self.call_count:
            self.call_count["offset_pwr_list"] += 1
        else:
            self.call_count["offset_pwr_list"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(clist_RDC).exists(), result)
        result = self._validate(Path(clist_MAP).exists(), result)
        Path(offs).touch()
        Path(ccp).touch()
        Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(int_filt in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        Path(ccs).touch()
        return result

    def geocode(self, lookup_table: str, data_in: str, width_in, data_out: str, width_out, nlines_out, interp_mode, dtype, lr_in, lr_out, n_ovr, rad_max, nintr):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "geocode", supplied_args))

        if "geocode" in self.call_count:
            self.call_count["geocode"] += 1
        else:
            self.call_count["geocode"] = 1

        result = self._validate(Path(lookup_table).exists(), result)
        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(interp_mode in valid_values, result)
        valid_values = [0, 1, 2, 3, 4, 5, 6]
        result = self._validate(dtype in valid_values, result)
        return result

    def par_TX_geo(self, annotation_XML: str, GeoTIFF: str, MLI_par: str, DEM_par: str, GEO: str, pol):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_TX_geo", supplied_args))

        if "par_TX_geo" in self.call_count:
            self.call_count["par_TX_geo"] += 1
        else:
            self.call_count["par_TX_geo"] = 1

        result = self._validate(Path(annotation_XML).exists(), result)
        result = self._validate(Path(GeoTIFF).exists(), result)
        Path(MLI_par).touch()
        Path(DEM_par).touch()
        Path(GEO).touch()
        return result

    def data2xyz(self, DEM_par: str, data: str, data_xyz: str, dflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "data2xyz", supplied_args))

        if "data2xyz" in self.call_count:
            self.call_count["data2xyz"] += 1
        else:
            self.call_count["data2xyz"] = 1

        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(data).exists(), result)
        Path(data_xyz).touch()
        valid_values = [0, 1]
        result = self._validate(dflg in valid_values, result)
        return result

    def ScanSAR_burst_diff_intf(self, SLC1_tab: str, SLC2R_tab: str, SIM_tab: str, DIFF_tab, SLCR_tab: str, DIFF_dir):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "ScanSAR_burst_diff_intf", supplied_args))

        if "ScanSAR_burst_diff_intf" in self.call_count:
            self.call_count["ScanSAR_burst_diff_intf"] += 1
        else:
            self.call_count["ScanSAR_burst_diff_intf"] = 1

        result = self._validate(Path(SLC1_tab).exists(), result)
        result = self._validate(Path(SLC2R_tab).exists(), result)
        result = self._validate(Path(SIM_tab).exists(), result)
        result = self._validate(Path(SLCR_tab).exists(), result)
        return result

    def dem_xyz(self, DEM_par: str, DEM: str, DEM_XYZ: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_xyz", supplied_args))

        if "dem_xyz" in self.call_count:
            self.call_count["dem_xyz"] += 1
        else:
            self.call_count["dem_xyz"] = 1

        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        Path(DEM_XYZ).touch()
        return result

    def look_vector(self, SLC_par: str, OFF_par: str, DEM_par: str, DEM: str, lv_theta: str, lv_phi: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "look_vector", supplied_args))

        if "look_vector" in self.call_count:
            self.call_count["look_vector"] += 1
        else:
            self.call_count["look_vector"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        Path(lv_theta).touch()
        Path(lv_phi).touch()
        return result

    def create_dem_par(self, DEM_par, SLC_par: str, terra_alt, delta_y, delta_x, EPSG, iflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "create_dem_par", supplied_args))

        if "create_dem_par" in self.call_count:
            self.call_count["create_dem_par"] += 1
        else:
            self.call_count["create_dem_par"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(iflg in valid_values, result)
        return result

    def dem_x_y_z(self, DEM_par: str, DEM: str, DEM_X: str, DEM_Y: str, DEM_Z: str, format_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_x_y_z", supplied_args))

        if "dem_x_y_z" in self.call_count:
            self.call_count["dem_x_y_z"] += 1
        else:
            self.call_count["dem_x_y_z"] = 1

        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        Path(DEM_X).touch()
        Path(DEM_Y).touch()
        Path(DEM_Z).touch()
        valid_values = [0, 1]
        result = self._validate(format_flag in valid_values, result)
        return result

    def ras_clist(self, clist: str, ras_in: str, ras_out: str, xsf, ysf, r, g, b, xs, zflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "ras_clist", supplied_args))

        if "ras_clist" in self.call_count:
            self.call_count["ras_clist"] += 1
        else:
            self.call_count["ras_clist"] = 1

        result = self._validate(Path(clist).exists(), result)
        result = self._validate(Path(ras_in).exists(), result)
        Path(ras_out).touch()
        return result

    def dispmap_sim(self, LV: str, DEM_par: str, disp_east: str, disp_north: str, disp_up: str, disp_LOS: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_sim", supplied_args))

        if "dispmap_sim" in self.call_count:
            self.call_count["dispmap_sim"] += 1
        else:
            self.call_count["dispmap_sim"] = 1

        result = self._validate(Path(LV).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(disp_east).exists(), result)
        result = self._validate(Path(disp_north).exists(), result)
        result = self._validate(Path(disp_up).exists(), result)
        Path(disp_LOS).touch()
        return result

    def dem_gradient(self, DEM_par: str, DEM: str, theta: str, phi: str, mag: str, type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_gradient", supplied_args))

        if "dem_gradient" in self.call_count:
            self.call_count["dem_gradient"] += 1
        else:
            self.call_count["dem_gradient"] = 1

        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        Path(theta).touch()
        Path(phi).touch()
        Path(mag).touch()
        valid_values = [0, 1]
        result = self._validate(type in valid_values, result)
        return result

    def offset_pwr_trackingm(self, MLI_1: str, MLI_2: str, DIFF_par: str, offs: str, ccp: str, rwin, azwin, offsets: str, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, lanczos, bw_frac, pflag, pltflg, ccs: str, std_mean):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_pwr_trackingm", supplied_args))

        if "offset_pwr_trackingm" in self.call_count:
            self.call_count["offset_pwr_trackingm"] += 1
        else:
            self.call_count["offset_pwr_trackingm"] = 1

        result = self._validate(Path(MLI_1).exists(), result)
        result = self._validate(Path(MLI_2).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(offs).touch()
        Path(ccp).touch()
        Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        Path(ccs).touch()
        return result

    def init_offsetm(self, MLI_1: str, MLI_2: str, DIFF_par, rlks, azlks, rpos, azpos, offr, offaz, thres, patch, cflag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "init_offsetm", supplied_args))

        if "init_offsetm" in self.call_count:
            self.call_count["init_offsetm"] += 1
        else:
            self.call_count["init_offsetm"] = 1

        result = self._validate(Path(MLI_1).exists(), result)
        result = self._validate(Path(MLI_2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(cflag in valid_values, result)
        return result

    def gc_map_fd(self, MLI_par: str, fdtab: str, DEM_par: str, DEM: str, DEM_seg_par, DEM_seg: str, lookup_table: str, lat_ovr, lon_ovr, sim_sar: str, u: str, v: str, inc: str, psi: str, pix: str, ls_map: str, frame, ls_mode, r_ovr):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map_fd", supplied_args))

        if "gc_map_fd" in self.call_count:
            self.call_count["gc_map_fd"] += 1
        else:
            self.call_count["gc_map_fd"] = 1

        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(fdtab).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        Path(DEM_seg).touch()
        Path(lookup_table).touch()
        Path(sim_sar).touch()
        Path(u).touch()
        Path(v).touch()
        Path(inc).touch()
        Path(psi).touch()
        Path(pix).touch()
        Path(ls_map).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(ls_mode in valid_values, result)
        return result

    def phase_sim_orb(self, SLC1_par: str, SLC2R_par: str, OFF_par: str, hgt: str, sim_unw: str, SLC_ref_par: str, definition: str, delta_t: str, int_mode: str, ph_mode):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "phase_sim_orb", supplied_args))

        if "phase_sim_orb" in self.call_count:
            self.call_count["phase_sim_orb"] += 1
        else:
            self.call_count["phase_sim_orb"] = 1

        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2R_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(hgt).exists(), result)
        Path(sim_unw).touch()
        result = self._validate(Path(SLC_ref_par).exists(), result)
        result = self._validate(Path(definition).exists(), result)
        result = self._validate(Path(delta_t).exists(), result)
        result = self._validate(Path(int_mode).exists(), result)
        valid_values = [0, 1]
        result = self._validate(ph_mode in valid_values, result)
        return result

    def gec_map(self, SLC_par: str, OFF_par: str, DEM_par: str, href: str, DEM_seg_par, lookup_table: str, lat_ovr, lon_ovr):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gec_map", supplied_args))

        if "gec_map" in self.call_count:
            self.call_count["gec_map"] += 1
        else:
            self.call_count["gec_map"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(href).exists(), result)
        Path(lookup_table).touch()
        return result

    def dh_map_orb(self, SLC1_par: str, SLC2R_par: str, OFF_par: str, hgt: str, dp: str, dpdh: str, dh: str, SLC_ref_par: str, int_mode: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dh_map_orb", supplied_args))

        if "dh_map_orb" in self.call_count:
            self.call_count["dh_map_orb"] += 1
        else:
            self.call_count["dh_map_orb"] = 1

        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2R_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(hgt).exists(), result)
        result = self._validate(Path(dp).exists(), result)
        Path(dpdh).touch()
        Path(dh).touch()
        result = self._validate(Path(SLC_ref_par).exists(), result)
        result = self._validate(Path(int_mode).exists(), result)
        return result

    def offset_subm(self, offs: str, DIFF_par: str, offs_sub: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_subm", supplied_args))

        if "offset_subm" in self.call_count:
            self.call_count["offset_subm"] += 1
        else:
            self.call_count["offset_subm"] = 1

        result = self._validate(Path(offs).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(offs_sub).touch()
        return result

    def offset_trackingm(self, offs: str, snr: str, MLI_par: str, DIFF_par: str, coffs_map: str, coffsets: str, mode, thres, poly_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_trackingm", supplied_args))

        if "offset_trackingm" in self.call_count:
            self.call_count["offset_trackingm"] += 1
        else:
            self.call_count["offset_trackingm"] = 1

        result = self._validate(Path(offs).exists(), result)
        result = self._validate(Path(snr).exists(), result)
        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(coffs_map).touch()
        Path(coffsets).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(poly_flag in valid_values, result)
        return result

    def comb_interfs(self, int_1, int_2, base_1, base_2, factor_1, factor_2, width, combi_int, combi_base, sm, Only, The):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "comb_interfs", supplied_args))

        if "comb_interfs" in self.call_count:
            self.call_count["comb_interfs"] += 1
        else:
            self.call_count["comb_interfs"] = 1

        return result

    def gc_insar(self, SLC_par: str, OFF_par: str, hgt: str, DEM_par: str, lookup_table: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_insar", supplied_args))

        if "gc_insar" in self.call_count:
            self.call_count["gc_insar"] += 1
        else:
            self.call_count["gc_insar"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(hgt).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        Path(lookup_table).touch()
        return result

    def par_KS_geo(self, HDF5: str, trunk: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_KS_geo", supplied_args))

        if "par_KS_geo" in self.call_count:
            self.call_count["par_KS_geo"] += 1
        else:
            self.call_count["par_KS_geo"] = 1

        result = self._validate(Path(HDF5).exists(), result)
        result = self._validate(Path(trunk).exists(), result)
        return result

    def offset_pwrm(self, MLI_1: str, MLI_2: str, DIFF_par: str, offs: str, ccp: str, rwin, azwin, offsets: str, n_ovr, nr, naz, thres, lanczos, bw_frac, pflag, pltflg, ccs: str, std_mean):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_pwrm", supplied_args))

        if "offset_pwrm" in self.call_count:
            self.call_count["offset_pwrm"] += 1
        else:
            self.call_count["offset_pwrm"] = 1

        result = self._validate(Path(MLI_1).exists(), result)
        result = self._validate(Path(MLI_2).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(offs).touch()
        Path(ccp).touch()
        Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        Path(ccs).touch()
        return result

    def coord_to_sarpix_list(self, SLC_par: str, OFF_par: str, DEM_par: str, MAP_coord: str, SAR_coord: str, DIFF_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "coord_to_sarpix_list", supplied_args))

        if "coord_to_sarpix_list" in self.call_count:
            self.call_count["coord_to_sarpix_list"] += 1
        else:
            self.call_count["coord_to_sarpix_list"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(MAP_coord).exists(), result)
        Path(SAR_coord).touch()
        result = self._validate(Path(DIFF_par).exists(), result)
        return result

    def offset_fitm(self, offs: str, ccp: str, DIFF_par: str, coffs: str, coffsets: str, thres, npoly, interact_mode):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_fitm", supplied_args))

        if "offset_fitm" in self.call_count:
            self.call_count["offset_fitm"] += 1
        else:
            self.call_count["offset_fitm"] = 1

        result = self._validate(Path(offs).exists(), result)
        result = self._validate(Path(ccp).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(coffs).touch()
        Path(coffsets).touch()
        valid_values = [0, 1]
        result = self._validate(interact_mode in valid_values, result)
        return result

    def dispmap_ENU(self, LV_tab: str, DISP_tab: str, SIGMA_tab: str, DEM_par: str, disp_east: str, disp_north: str, disp_up: str, sigma_east: str, sigma_north: str, sigma_up: str, chi2: str, min_obs, tol):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_ENU", supplied_args))

        if "dispmap_ENU" in self.call_count:
            self.call_count["dispmap_ENU"] += 1
        else:
            self.call_count["dispmap_ENU"] = 1

        result = self._validate(Path(LV_tab).exists(), result)
        result = self._validate(Path(DISP_tab).exists(), result)
        result = self._validate(Path(SIGMA_tab).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        Path(disp_east).touch()
        Path(disp_north).touch()
        Path(disp_up).touch()
        Path(sigma_east).touch()
        Path(sigma_north).touch()
        Path(sigma_up).touch()
        Path(chi2).touch()
        return result

    def dem_coord(self, DEM_par: str, east: str, north: str, dtype):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_coord", supplied_args))

        if "dem_coord" in self.call_count:
            self.call_count["dem_coord"] += 1
        else:
            self.call_count["dem_coord"] = 1

        result = self._validate(Path(DEM_par).exists(), result)
        Path(east).touch()
        Path(north).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def rotate_image(self, data_in: str, width_in, angle, data_out: str, width_out, nlines_out, interp_mode, dtype, order):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "rotate_image", supplied_args))

        if "rotate_image" in self.call_count:
            self.call_count["rotate_image"] += 1
        else:
            self.call_count["rotate_image"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6, 7]
        result = self._validate(interp_mode in valid_values, result)
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        return result

    def geocode_back(self, data_in: str, width_in, lookup_table: str, data_out: str, width_out, nlines_out, interp_mode, dtype, lr_in, lr_out, order, e_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "geocode_back", supplied_args))

        if "geocode_back" in self.call_count:
            self.call_count["geocode_back"] += 1
        else:
            self.call_count["geocode_back"] = 1

        result = self._validate(Path(data_in).exists(), result)
        result = self._validate(Path(lookup_table).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(e_flag in valid_values, result)
        return result

    def sarpix_coord(self, SLC_par: str, OFF_par: str, DEM_par: str, azlin, rpix, ref_hgt):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "sarpix_coord", supplied_args))

        if "sarpix_coord" in self.call_count:
            self.call_count["sarpix_coord"] += 1
        else:
            self.call_count["sarpix_coord"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        return result

    def resamp_image(self, data_in: str, width_in, xscale, yscale, data_out: str, width_out, nlines_out, interp_mode, dtype, order):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "resamp_image", supplied_args))

        if "resamp_image" in self.call_count:
            self.call_count["resamp_image"] += 1
        else:
            self.call_count["resamp_image"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        return result

    def pol2rec(self, data1: str, SLC_par1: str, data2: str, SLC_par2: str, pix_size: str, dtype, mode, xmin, nx, ymin, ny, rmax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "pol2rec", supplied_args))

        if "pol2rec" in self.call_count:
            self.call_count["pol2rec"] += 1
        else:
            self.call_count["pol2rec"] = 1

        result = self._validate(Path(data1).exists(), result)
        result = self._validate(Path(SLC_par1).exists(), result)
        Path(data2).touch()
        Path(SLC_par2).touch()
        Path(pix_size).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def WSS_interp(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, DIFF_par: str, SLC_2R: str, SLC2R_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "WSS_interp", supplied_args))

        if "WSS_interp" in self.call_count:
            self.call_count["WSS_interp"] += 1
        else:
            self.call_count["WSS_interp"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(DIFF_par).exists(), result)
        Path(SLC_2R).touch()
        Path(SLC2R_par).touch()
        return result

    def gc_GPRI_map(self, MLI_par: str, DEM_par: str, DEM: str, DEM_seg_par, DEM_seg: str, lookup_table: str, lat_ovr, lon_ovr, sim_sar: str, lv_theta: str, lv_phi: str, u: str, v: str, inc: str, psi: str, pix: str, ls_map: str, frame):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_GPRI_map", supplied_args))

        if "gc_GPRI_map" in self.call_count:
            self.call_count["gc_GPRI_map"] += 1
        else:
            self.call_count["gc_GPRI_map"] = 1

        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(DEM).exists(), result)
        Path(DEM_seg).touch()
        Path(lookup_table).touch()
        Path(sim_sar).touch()
        Path(lv_theta).touch()
        Path(lv_phi).touch()
        Path(u).touch()
        Path(v).touch()
        Path(inc).touch()
        Path(psi).touch()
        Path(pix).touch()
        Path(ls_map).touch()
        return result

    def dop_mlcc(self, SAR_par: str, PROC_par: str, signal_data: str, plot_data: str, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "dop_mlcc", supplied_args))

        if "dop_mlcc" in self.call_count:
            self.call_count["dop_mlcc"] += 1
        else:
            self.call_count["dop_mlcc"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(plot_data).touch()
        return result

    def ERS_proc_ASF_2000(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_ASF_2000", supplied_args))

        if "ERS_proc_ASF_2000" in self.call_count:
            self.call_count["ERS_proc_ASF_2000"] += 1
        else:
            self.call_count["ERS_proc_ASF_2000"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def ERS_proc_ASF(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_ASF", supplied_args))

        if "ERS_proc_ASF" in self.call_count:
            self.call_count["ERS_proc_ASF"] += 1
        else:
            self.call_count["ERS_proc_ASF"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def doppler_2d(self, SAR_par: str, PROC_par: str, signal_data: str, dop2d: str, loff, blsz, nbl, a2_flg, b0_flg, b1_flg, c0_flg, ambig_flag, namb):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "doppler_2d", supplied_args))

        if "doppler_2d" in self.call_count:
            self.call_count["doppler_2d"] += 1
        else:
            self.call_count["doppler_2d"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(dop2d).touch()
        return result

    def ERS_proc_ESA(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_ESA", supplied_args))

        if "ERS_proc_ESA" in self.call_count:
            self.call_count["ERS_proc_ESA"] += 1
        else:
            self.call_count["ERS_proc_ESA"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def az_proc(self, SAR_par: str, PROC_par: str, rc_data: str, SLC: str, az_patch, SLC_format, cal_fact, SLC_type, kaiser, npatch):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "az_proc", supplied_args))

        if "az_proc" in self.call_count:
            self.call_count["az_proc"] += 1
        else:
            self.call_count["az_proc"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(rc_data).exists(), result)
        Path(SLC).touch()
        valid_values = [0, 1]
        result = self._validate(SLC_format in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(SLC_type in valid_values, result)
        return result

    def dop_interf(self, SAR_par1_in: str, PROC_par1_in: str, PROC_par2_in: str, PROC_par1_out: str, PROC_par2_out: str, dop: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "dop_interf", supplied_args))

        if "dop_interf" in self.call_count:
            self.call_count["dop_interf"] += 1
        else:
            self.call_count["dop_interf"] = 1

        result = self._validate(Path(SAR_par1_in).exists(), result)
        result = self._validate(Path(PROC_par1_in).exists(), result)
        result = self._validate(Path(PROC_par2_in).exists(), result)
        Path(PROC_par1_out).touch()
        Path(PROC_par2_out).touch()
        Path(dop).touch()
        return result

    def CS_proc(self, HDF5: str, SAR_par: str, PROC_par: str, raw_out: str, loff, nl):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "CS_proc", supplied_args))

        if "CS_proc" in self.call_count:
            self.call_count["CS_proc"] += 1
        else:
            self.call_count["CS_proc"] = 1

        result = self._validate(Path(HDF5).exists(), result)
        Path(SAR_par).touch()
        Path(PROC_par).touch()
        Path(raw_out).touch()
        return result

    def azsp_SLC(self, SAR_par: str, PROC_par: str, SAR_data: str, spectrum, loff, roff, nsub, data_format, pltflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "azsp_SLC", supplied_args))

        if "azsp_SLC" in self.call_count:
            self.call_count["azsp_SLC"] += 1
        else:
            self.call_count["azsp_SLC"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(SAR_data).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_format in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def pre_rc_JERS(self, SAR_par: str, PROC_par: str, rspec: str, signal_data: str, rc_data: str, prefilt_dec, kaiser, filt_lm):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "pre_rc_JERS", supplied_args))

        if "pre_rc_JERS" in self.call_count:
            self.call_count["pre_rc_JERS"] += 1
        else:
            self.call_count["pre_rc_JERS"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(rspec).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(rc_data).touch()
        return result

    def ERS_proc_NASDA(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_NASDA", supplied_args))

        if "ERS_proc_NASDA" in self.call_count:
            self.call_count["ERS_proc_NASDA"] += 1
        else:
            self.call_count["ERS_proc_NASDA"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def ERS_proc_ARG(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_ARG", supplied_args))

        if "ERS_proc_ARG" in self.call_count:
            self.call_count["ERS_proc_ARG"] += 1
        else:
            self.call_count["ERS_proc_ARG"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def PALSAR_proc(self, CEOS_SAR_leader: str, SAR_par: str, PROC_par: str, CEOS_raw_data: str, raw_out: str, TX_POL, RX_POL):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "PALSAR_proc", supplied_args))

        if "PALSAR_proc" in self.call_count:
            self.call_count["PALSAR_proc"] += 1
        else:
            self.call_count["PALSAR_proc"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(SAR_par).touch()
        Path(PROC_par).touch()
        result = self._validate(Path(CEOS_raw_data).exists(), result)
        Path(raw_out).touch()
        valid_values = [0, 1]
        result = self._validate(TX_POL in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(RX_POL in valid_values, result)
        return result

    def rspec_real(self, SAR_par: str, PROC_par: str, signal_data: str, range_spec: str, loff, nlspec, nrfft, pltflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "rspec_real", supplied_args))

        if "rspec_real" in self.call_count:
            self.call_count["rspec_real"] += 1
        else:
            self.call_count["rspec_real"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(range_spec).touch()
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def rspec_IQ(self, SAR_par: str, PROC_par: str, signal_data: str, range_spec: str, loff, nlspec, nrfft, roff, pltflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "rspec_IQ", supplied_args))

        if "rspec_IQ" in self.call_count:
            self.call_count["rspec_IQ"] += 1
        else:
            self.call_count["rspec_IQ"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(range_spec).touch()
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def af(self, SAR_par: str, PROC_par: str, SLC: str, rwin, azwin, dr, daz, thres, update_flg, a1_flg, b0_flg, offsets: str, dac_flg, n_ovr, roff, azoff):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "af", supplied_args))

        if "af" in self.call_count:
            self.call_count["af"] += 1
        else:
            self.call_count["af"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(SLC).exists(), result)
        valid_values = [0, 1]
        result = self._validate(update_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(a1_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(b0_flg in valid_values, result)
        Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(dac_flg in valid_values, result)
        return result

    def hist_IQ(self, SAR_par: str, PROC_par: str, signal_data: str, histogram: str, loff, nl):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "hist_IQ", supplied_args))

        if "hist_IQ" in self.call_count:
            self.call_count["hist_IQ"] += 1
        else:
            self.call_count["hist_IQ"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(histogram).touch()
        return result

    def PALSAR_burst_sync(self, SAR_par1: str, PROC_par1: str, raw1: str, SAR_par2: str, PROC_par2: str, raw2: str, PROC_par1_out: str, raw1_out: str, PROC_par2_out: str, raw2_out: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "PALSAR_burst_sync", supplied_args))

        if "PALSAR_burst_sync" in self.call_count:
            self.call_count["PALSAR_burst_sync"] += 1
        else:
            self.call_count["PALSAR_burst_sync"] = 1

        result = self._validate(Path(SAR_par1).exists(), result)
        result = self._validate(Path(PROC_par1).exists(), result)
        result = self._validate(Path(raw1).exists(), result)
        result = self._validate(Path(SAR_par2).exists(), result)
        result = self._validate(Path(PROC_par2).exists(), result)
        result = self._validate(Path(raw2).exists(), result)
        Path(PROC_par1_out).touch()
        Path(raw1_out).touch()
        Path(PROC_par2_out).touch()
        Path(raw2_out).touch()
        return result

    def swap_IQ(self, SAR_par: str, raw_IQ: str, raw_IQ_swap: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "swap_IQ", supplied_args))

        if "swap_IQ" in self.call_count:
            self.call_count["swap_IQ"] += 1
        else:
            self.call_count["swap_IQ"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(raw_IQ).exists(), result)
        Path(raw_IQ_swap).touch()
        return result

    def dop_ambig(self, SAR_par: str, PROC_par: str, signal_data: str, algorithm, loff, output_plot: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "dop_ambig", supplied_args))

        if "dop_ambig" in self.call_count:
            self.call_count["dop_ambig"] += 1
        else:
            self.call_count["dop_ambig"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        valid_values = [1, 2]
        result = self._validate(algorithm in valid_values, result)
        Path(output_plot).touch()
        return result

    def PRC_proc(self, PROC_par: str, PRC, nstate):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "PRC_proc", supplied_args))

        if "PRC_proc" in self.call_count:
            self.call_count["PRC_proc"] += 1
        else:
            self.call_count["PRC_proc"] = 1

        result = self._validate(Path(PROC_par).exists(), result)
        return result

    def rc_fmcw(self, SAR_par: str, PROC_par: str, signal_data: str, rc_data: str, nrc_off, nrc_samp, loff, nl, kaiser):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "rc_fmcw", supplied_args))

        if "rc_fmcw" in self.call_count:
            self.call_count["rc_fmcw"] += 1
        else:
            self.call_count["rc_fmcw"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(rc_data).touch()
        return result

    def RSAT_lks(self, SLC_PROC_par: str, MLI_PROC_par: str, SLC_image: str, MLI_image: str, kaiser):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "RSAT_lks", supplied_args))

        if "RSAT_lks" in self.call_count:
            self.call_count["RSAT_lks"] += 1
        else:
            self.call_count["RSAT_lks"] = 1

        result = self._validate(Path(SLC_PROC_par).exists(), result)
        Path(MLI_PROC_par).touch()
        result = self._validate(Path(SLC_image).exists(), result)
        Path(MLI_image).touch()
        return result

    def extract_psd(self, range_spec, num, output):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "extract_psd", supplied_args))

        if "extract_psd" in self.call_count:
            self.call_count["extract_psd"] += 1
        else:
            self.call_count["extract_psd"] = 1

        return result

    def doppler_real(self, SAR_par: str, PROC_par: str, signal_data: str, doppler: str, loff, nsub, ambig_flag, namb):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "doppler_real", supplied_args))

        if "doppler_real" in self.call_count:
            self.call_count["doppler_real"] += 1
        else:
            self.call_count["doppler_real"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(doppler).touch()
        return result

    def pre_rc(self, SAR_par: str, PROC_par: str, signal_data: str, rc_data: str, prefilt_dec, loff, nl, nr_samp, kaiser, filt_lm, nr_ext, fr_ext, pre_ext, post_ext, RFI_filt, RFI_thres, fc_offset, win_bw):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "pre_rc", supplied_args))

        if "pre_rc" in self.call_count:
            self.call_count["pre_rc"] += 1
        else:
            self.call_count["pre_rc"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(rc_data).touch()
        valid_values = [0, 1]
        result = self._validate(RFI_filt in valid_values, result)
        return result

    def PALSAR_proc_WB(self, CEOS_SAR_leader: str, SAR_par: str, PROC_par: str, CEOS_raw_data: str, beam: str, raw_out: str, prf: str, wflg: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "PALSAR_proc_WB", supplied_args))

        if "PALSAR_proc_WB" in self.call_count:
            self.call_count["PALSAR_proc_WB"] += 1
        else:
            self.call_count["PALSAR_proc_WB"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(SAR_par).touch()
        Path(PROC_par).touch()
        result = self._validate(Path(CEOS_raw_data).exists(), result)
        result = self._validate(Path(beam).exists(), result)
        Path(raw_out).touch()
        result = self._validate(Path(prf).exists(), result)
        result = self._validate(Path(wflg).exists(), result)
        return result

    def ORRM_proc(self, PROC_par, ORRM: str, nstate):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ORRM_proc", supplied_args))

        if "ORRM_proc" in self.call_count:
            self.call_count["ORRM_proc"] += 1
        else:
            self.call_count["ORRM_proc"] = 1

        result = self._validate(Path(ORRM).exists(), result)
        return result

    def ERS_proc_UK(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_UK", supplied_args))

        if "ERS_proc_UK" in self.call_count:
            self.call_count["ERS_proc_UK"] += 1
        else:
            self.call_count["ERS_proc_UK"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def RSAT_raw(self, CEOS_ldr: str, SAR_par: str, PROC_par: str, raw_data_files: str, raw_out: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "RSAT_raw", supplied_args))

        if "RSAT_raw" in self.call_count:
            self.call_count["RSAT_raw"] += 1
        else:
            self.call_count["RSAT_raw"] = 1

        result = self._validate(Path(CEOS_ldr).exists(), result)
        Path(SAR_par).touch()
        Path(PROC_par).touch()
        result = self._validate(Path(raw_data_files).exists(), result)
        Path(raw_out).touch()
        return result

    def ERS_ENVISAT_proc(self, L0: str, SAR_par, PROC_par: str, raw: str, loff, nl, swst_flg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_ENVISAT_proc", supplied_args))

        if "ERS_ENVISAT_proc" in self.call_count:
            self.call_count["ERS_ENVISAT_proc"] += 1
        else:
            self.call_count["ERS_ENVISAT_proc"] = 1

        result = self._validate(Path(L0).exists(), result)
        Path(PROC_par).touch()
        Path(raw).touch()
        valid_values = [0, 1]
        result = self._validate(swst_flg in valid_values, result)
        return result

    def multi_SLC(self, SLC_PROC_par: str, MLI_PROC_par: str, SLC: str, MLI: str, rlks, azlks, slc_format):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "multi_SLC", supplied_args))

        if "multi_SLC" in self.call_count:
            self.call_count["multi_SLC"] += 1
        else:
            self.call_count["multi_SLC"] = 1

        result = self._validate(Path(SLC_PROC_par).exists(), result)
        Path(MLI_PROC_par).touch()
        result = self._validate(Path(SLC).exists(), result)
        Path(MLI).touch()
        valid_values = [0, 1]
        result = self._validate(slc_format in valid_values, result)
        return result

    def JERS_acs(self, USER_HEADER: str, SEG_DESCR: str, ORBIT_DATA: str, SENSOR_DATA: str, track: str, SAR_par: str, PROC_par: str, raw_out: str, loff, nl, nsx, fsx, terra_alt, deskew):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "JERS_acs", supplied_args))

        if "JERS_acs" in self.call_count:
            self.call_count["JERS_acs"] += 1
        else:
            self.call_count["JERS_acs"] = 1

        result = self._validate(Path(USER_HEADER).exists(), result)
        result = self._validate(Path(SEG_DESCR).exists(), result)
        result = self._validate(Path(ORBIT_DATA).exists(), result)
        result = self._validate(Path(SENSOR_DATA).exists(), result)
        Path(track).touch()
        Path(SAR_par).touch()
        Path(PROC_par).touch()
        Path(raw_out).touch()
        return result

    def DORIS_proc(self, PROC_par, DOR: str, nstate):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "DORIS_proc", supplied_args))

        if "DORIS_proc" in self.call_count:
            self.call_count["DORIS_proc"] += 1
        else:
            self.call_count["DORIS_proc"] = 1

        result = self._validate(Path(DOR).exists(), result)
        return result

    def ASAR_AP_proc(self, L0: str, INS: str, SAR_par1: str, SAR_par2: str, PROC_par1: str, PROC_par2: str, raw1: str, raw2: str, ant_gain1: str, ant_gain2: str, loff, nl, roff, nr, refer):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ASAR_AP_proc", supplied_args))

        if "ASAR_AP_proc" in self.call_count:
            self.call_count["ASAR_AP_proc"] += 1
        else:
            self.call_count["ASAR_AP_proc"] = 1

        result = self._validate(Path(L0).exists(), result)
        result = self._validate(Path(INS).exists(), result)
        Path(SAR_par1).touch()
        Path(SAR_par2).touch()
        Path(PROC_par1).touch()
        Path(PROC_par2).touch()
        Path(raw1).touch()
        Path(raw2).touch()
        result = self._validate(Path(ant_gain1).exists(), result)
        result = self._validate(Path(ant_gain2).exists(), result)
        return result

    def JERS_proc_ASF(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "JERS_proc_ASF", supplied_args))

        if "JERS_proc_ASF" in self.call_count:
            self.call_count["JERS_proc_ASF"] += 1
        else:
            self.call_count["JERS_proc_ASF"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def rc_real(self, SAR_par: str, PROC_par: str, signal_data: str, rc_data: str, loff, nl, kaiser, nr_ext, fr_ext, r_chirp: str, rfi_filt, rfi_thres):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "rc_real", supplied_args))

        if "rc_real" in self.call_count:
            self.call_count["rc_real"] += 1
        else:
            self.call_count["rc_real"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(rc_data).touch()
        result = self._validate(Path(r_chirp).exists(), result)
        return result

    def ERS_proc_CRISP(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_CRISP", supplied_args))

        if "ERS_proc_CRISP" in self.call_count:
            self.call_count["ERS_proc_CRISP"] += 1
        else:
            self.call_count["ERS_proc_CRISP"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def prefilt(self, SAR_par: str, PROC_par: str, rc_data: str, prefilt_out: str, prefilt_dec, filt_lm):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "prefilt", supplied_args))

        if "prefilt" in self.call_count:
            self.call_count["prefilt"] += 1
        else:
            self.call_count["prefilt"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(rc_data).exists(), result)
        Path(prefilt_out).touch()
        return result

    def doppler(self, SAR_par: str, PROC_par: str, signal_data: str, doppler: str, loff, nsub, ambig_flag, namb, order, pltflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "doppler", supplied_args))

        if "doppler" in self.call_count:
            self.call_count["doppler"] += 1
        else:
            self.call_count["doppler"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(doppler).touch()
        valid_values = [0, 1, 2]
        result = self._validate(ambig_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def DELFT_proc2(self, PROC_par: str, DELFT_dir, nstate, interval, ODR):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "DELFT_proc2", supplied_args))

        if "DELFT_proc2" in self.call_count:
            self.call_count["DELFT_proc2"] += 1
        else:
            self.call_count["DELFT_proc2"] = 1

        result = self._validate(Path(PROC_par).exists(), result)
        return result

    def copy(self, infile: str, outfile: str, lbytes, start, nlines, offset, file_ldr, offb, nbyte):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "copy", supplied_args))

        if "copy" in self.call_count:
            self.call_count["copy"] += 1
        else:
            self.call_count["copy"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        return result

    def ERS_proc_ESRIN_ACS(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_ESRIN_ACS", supplied_args))

        if "ERS_proc_ESRIN_ACS" in self.call_count:
            self.call_count["ERS_proc_ESRIN_ACS"] += 1
        else:
            self.call_count["ERS_proc_ESRIN_ACS"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def ERS_proc_ASF_91(self, CEOS_SAR_leader: str, CEOS_trailer: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_ASF_91", supplied_args))

        if "ERS_proc_ASF_91" in self.call_count:
            self.call_count["ERS_proc_ASF_91"] += 1
        else:
            self.call_count["ERS_proc_ASF_91"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        result = self._validate(Path(CEOS_trailer).exists(), result)
        Path(PROC_par).touch()
        return result

    def ORB_prop(self, PROC_par: str, nstate, interval, extra):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ORB_prop", supplied_args))

        if "ORB_prop" in self.call_count:
            self.call_count["ORB_prop"] += 1
        else:
            self.call_count["ORB_prop"] = 1

        result = self._validate(Path(PROC_par).exists(), result)
        return result

    def ERS_fix(self, ERS_PAF, SAR_par: str, PROC_par: str, cc_flag, raw, output_file: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_fix", supplied_args))

        if "ERS_fix" in self.call_count:
            self.call_count["ERS_fix"] += 1
        else:
            self.call_count["ERS_fix"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        Path(output_file).touch()
        return result

    def ERS_proc_CCRS(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_CCRS", supplied_args))

        if "ERS_proc_CCRS" in self.call_count:
            self.call_count["ERS_proc_CCRS"] += 1
        else:
            self.call_count["ERS_proc_CCRS"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def ptarg(self, SLC: str, width: str, r_samp: str, az_samp: str, ptr_image: str, r_plot: str, az_plot: str, data_format, win, pltflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ptarg", supplied_args))

        if "ptarg" in self.call_count:
            self.call_count["ptarg"] += 1
        else:
            self.call_count["ptarg"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(width).exists(), result)
        result = self._validate(Path(r_samp).exists(), result)
        result = self._validate(Path(az_samp).exists(), result)
        Path(ptr_image).touch()
        Path(r_plot).touch()
        Path(az_plot).touch()
        valid_values = [0, 1]
        result = self._validate(data_format in valid_values, result)
        valid_values = [0, 1, 2]
        result = self._validate(pltflg in valid_values, result)
        return result

    def multi_GRD_SLC(self, SLC_PROC_par: str, GRD_PROC_par: str, SLC_image: str, GRD_image: str, rlks, azlks, interp_mode, sample_spacing, gr_start, t_start, t_end):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "multi_GRD_SLC", supplied_args))

        if "multi_GRD_SLC" in self.call_count:
            self.call_count["multi_GRD_SLC"] += 1
        else:
            self.call_count["multi_GRD_SLC"] = 1

        result = self._validate(Path(SLC_PROC_par).exists(), result)
        Path(GRD_PROC_par).touch()
        result = self._validate(Path(SLC_image).exists(), result)
        Path(GRD_image).touch()
        valid_values = [0, 1]
        result = self._validate(interp_mode in valid_values, result)
        return result

    def create_sar_par(self, SAR_par):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "create_sar_par", supplied_args))

        if "create_sar_par" in self.call_count:
            self.call_count["create_sar_par"] += 1
        else:
            self.call_count["create_sar_par"] = 1

        return result

    def rspec_JERS(self, SAR_par: str, PROC_par: str, signal_data: str, range_spec: str, nr_samp, nl_spec, loff, nlines, nr_ext, fr_ext):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "rspec_JERS", supplied_args))

        if "rspec_JERS" in self.call_count:
            self.call_count["rspec_JERS"] += 1
        else:
            self.call_count["rspec_JERS"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(range_spec).touch()
        return result

    def azsp_IQ(self, SAR_par: str, PROC_par: str, signal_data: str, spectrum: str, loff, roff, nsub, ambig_flg, namb, pltflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "azsp_IQ", supplied_args))

        if "azsp_IQ" in self.call_count:
            self.call_count["azsp_IQ"] += 1
        else:
            self.call_count["azsp_IQ"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(spectrum).touch()
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def ERS_proc_ACRES(self, CEOS_SAR_leader: str, PROC_par: str, type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_ACRES", supplied_args))

        if "ERS_proc_ACRES" in self.call_count:
            self.call_count["ERS_proc_ACRES"] += 1
        else:
            self.call_count["ERS_proc_ACRES"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        valid_values = [0, 1]
        result = self._validate(type in valid_values, result)
        return result

    def pre_rc_RSAT(self, SAR_par: str, PROC_par: str, signal_data: str, rc_data: str, prefilt_dec, loff, nl, nr_samp, kaiser, filt_lm, nr_ext, fr_ext):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "pre_rc_RSAT", supplied_args))

        if "pre_rc_RSAT" in self.call_count:
            self.call_count["pre_rc_RSAT"] += 1
        else:
            self.call_count["pre_rc_RSAT"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(signal_data).exists(), result)
        Path(rc_data).touch()
        return result

    def autof(self, SAR_par: str, PROC_par: str, rc_data: str, autofocus: str, SNR_min, prefilter, auto_az, az_offset, auto_bins, dop_ambig):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "autof", supplied_args))

        if "autof" in self.call_count:
            self.call_count["autof"] += 1
        else:
            self.call_count["autof"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        result = self._validate(Path(rc_data).exists(), result)
        Path(autofocus).touch()
        valid_values = [0, 1]
        result = self._validate(dop_ambig in valid_values, result)
        return result

    def cat_raw(self, RAW_list: str, SAR_par: str, PROC_par: str, RAW_out: str, fill, loff, nl):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "cat_raw", supplied_args))

        if "cat_raw" in self.call_count:
            self.call_count["cat_raw"] += 1
        else:
            self.call_count["cat_raw"] = 1

        result = self._validate(Path(RAW_list).exists(), result)
        Path(SAR_par).touch()
        Path(PROC_par).touch()
        Path(RAW_out).touch()
        valid_values = [0, 1]
        result = self._validate(fill in valid_values, result)
        return result

    def ERS_proc_ASI(self, CEOS_SAR_leader: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_ASI", supplied_args))

        if "ERS_proc_ASI" in self.call_count:
            self.call_count["ERS_proc_ASI"] += 1
        else:
            self.call_count["ERS_proc_ASI"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PROC_par).touch()
        return result

    def JERS_proc(self, CEOS_SAR_ldr: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "JERS_proc", supplied_args))

        if "JERS_proc" in self.call_count:
            self.call_count["JERS_proc"] += 1
        else:
            self.call_count["JERS_proc"] = 1

        result = self._validate(Path(CEOS_SAR_ldr).exists(), result)
        Path(PROC_par).touch()
        return result

    def create_proc_par(self, SAR_par: str, PROC_par):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "create_proc_par", supplied_args))

        if "create_proc_par" in self.call_count:
            self.call_count["create_proc_par"] += 1
        else:
            self.call_count["create_proc_par"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        return result

    def SIRC_proc(self, CEOS_SAR_leader: str, SAR_par: str, PROC_par: str, UTC_MET):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "SIRC_proc", supplied_args))

        if "SIRC_proc" in self.call_count:
            self.call_count["SIRC_proc"] += 1
        else:
            self.call_count["SIRC_proc"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(SAR_par).touch()
        Path(PROC_par).touch()
        return result

    def ASAR_IM_proc(self, L0: str, INS: str, SAR_par: str, PROC_par: str, raw: str, ant_gain: str, loff, nl):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ASAR_IM_proc", supplied_args))

        if "ASAR_IM_proc" in self.call_count:
            self.call_count["ASAR_IM_proc"] += 1
        else:
            self.call_count["ASAR_IM_proc"] = 1

        result = self._validate(Path(L0).exists(), result)
        result = self._validate(Path(INS).exists(), result)
        Path(SAR_par).touch()
        Path(PROC_par).touch()
        Path(raw).touch()
        result = self._validate(Path(ant_gain).exists(), result)
        return result

    def dishgt(self, hgt: str, pwr: str, width, start_hgt, start_pwr, nlines, m_cycle, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dishgt", supplied_args))

        if "dishgt" in self.call_count:
            self.call_count["dishgt"] += 1
        else:
            self.call_count["dishgt"] = 1

        result = self._validate(Path(hgt).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        return result

    def rasdt_pwr(self, data: str, pwr: str, width, start_data, start_pwr, nlines, pixavr, pixavaz, cycle, scale, exp, LR, rasf: str, cc: str, start_cc, cc_min):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasdt_pwr", supplied_args))

        if "rasdt_pwr" in self.call_count:
            self.call_count["rasdt_pwr"] += 1
        else:
            self.call_count["rasdt_pwr"] = 1

        result = self._validate(Path(data).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        result = self._validate(Path(cc).exists(), result)
        return result

    def rasdt_cmap(self, data: str, pwr: str, width, start_data, start_pwr, nlines, pixavr, pixavaz, min, max, mflg, cmap, scale, exp, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasdt_cmap", supplied_args))

        if "rasdt_cmap" in self.call_count:
            self.call_count["rasdt_cmap"] += 1
        else:
            self.call_count["rasdt_cmap"] = 1

        result = self._validate(Path(data).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        valid_values = [0, 1]
        result = self._validate(mflg in valid_values, result)
        Path(rasf).touch()
        return result

    def dis2hgt(self, hgt1: str, hgt2: str, width1, width2, start, nlines, roff, azoff, m_cycle):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2hgt", supplied_args))

        if "dis2hgt" in self.call_count:
            self.call_count["dis2hgt"] += 1
        else:
            self.call_count["dis2hgt"] = 1

        result = self._validate(Path(hgt1).exists(), result)
        result = self._validate(Path(hgt2).exists(), result)
        return result

    def gcp_ras(self, ras: str, GCP: str, mag, win_sz):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "gcp_ras", supplied_args))

        if "gcp_ras" in self.call_count:
            self.call_count["gcp_ras"] += 1
        else:
            self.call_count["gcp_ras"] = 1

        result = self._validate(Path(ras).exists(), result)
        Path(GCP).touch()
        return result

    def cpx_math(self, d1: str, d2: str, d_out: str, width, mode, roff, loff, nr, nl, c_re, c_im, zflg, rflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "cpx_math", supplied_args))

        if "cpx_math" in self.call_count:
            self.call_count["cpx_math"] += 1
        else:
            self.call_count["cpx_math"] = 1

        result = self._validate(Path(d1).exists(), result)
        result = self._validate(Path(d2).exists(), result)
        Path(d_out).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def disflag(self, flag: str, width, start, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disflag", supplied_args))

        if "disflag" in self.call_count:
            self.call_count["disflag"] += 1
        else:
            self.call_count["disflag"] = 1

        result = self._validate(Path(flag).exists(), result)
        return result

    def set_value(self, PAR_in: str, PAR_out: str, keyword, value, new_key):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "set_value", supplied_args))

        if "set_value" in self.call_count:
            self.call_count["set_value"] += 1
        else:
            self.call_count["set_value"] = 1

        result = self._validate(Path(PAR_in).exists(), result)
        Path(PAR_out).touch()
        valid_values = [0, 1]
        result = self._validate(new_key in valid_values, result)
        return result

    def cpx_to_real(self, cpx: str, real: str, width, type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "cpx_to_real", supplied_args))

        if "cpx_to_real" in self.call_count:
            self.call_count["cpx_to_real"] += 1
        else:
            self.call_count["cpx_to_real"] = 1

        result = self._validate(Path(cpx).exists(), result)
        Path(real).touch()
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(type in valid_values, result)
        return result

    def rascc(self, cc: str, pwr: str, width, start_cc, start_pwr, nlines, pixavr, pixavaz, cmin, cmax, scale, exp, LR, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rascc", supplied_args))

        if "rascc" in self.call_count:
            self.call_count["rascc"] += 1
        else:
            self.call_count["rascc"] = 1

        result = self._validate(Path(cc).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        return result

    def disgbyte(self, image: str, width, start, nlines, scale):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disgbyte", supplied_args))

        if "disgbyte" in self.call_count:
            self.call_count["disgbyte"] += 1
        else:
            self.call_count["disgbyte"] = 1

        result = self._validate(Path(image).exists(), result)
        return result

    def dis_dB(self, pwr: str, width, start, nlines, min_dB, max_dB):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis_dB", supplied_args))

        if "dis_dB" in self.call_count:
            self.call_count["dis_dB"] += 1
        else:
            self.call_count["dis_dB"] = 1

        result = self._validate(Path(pwr).exists(), result)
        return result

    def discc(self, cc: str, pwr: str, width, start_cc, start_pwr, nlines, cmin, cmax, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "discc", supplied_args))

        if "discc" in self.call_count:
            self.call_count["discc"] += 1
        else:
            self.call_count["discc"] = 1

        result = self._validate(Path(cc).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        return result

    def svg_map(self, image: str, dem_par: str, svg: str, font, fsize, color, gcolor, majorx, majory, minorx, minory, thick, grid, gopac, gdash):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "svg_map", supplied_args))

        if "svg_map" in self.call_count:
            self.call_count["svg_map"] += 1
        else:
            self.call_count["svg_map"] = 1

        result = self._validate(Path(image).exists(), result)
        result = self._validate(Path(dem_par).exists(), result)
        Path(svg).touch()
        return result

    def gbyte2float(self, infile: str, outfile: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "gbyte2float", supplied_args))

        if "gbyte2float" in self.call_count:
            self.call_count["gbyte2float"] += 1
        else:
            self.call_count["gbyte2float"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        return result

    def real_to_cpx(self, data1: str, data2: str, cpx: str, width, type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "real_to_cpx", supplied_args))

        if "real_to_cpx" in self.call_count:
            self.call_count["real_to_cpx"] += 1
        else:
            self.call_count["real_to_cpx"] = 1

        result = self._validate(Path(data1).exists(), result)
        result = self._validate(Path(data2).exists(), result)
        Path(cpx).touch()
        return result

    def rasmph_pwr(self, cpx: str, pwr: str, width, start_cpx, start_pwr, nlines, pixavr, pixavaz, scale, exp, LR, rasf: str, cc: str, start_cc, cc_min):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasmph_pwr", supplied_args))

        if "rasmph_pwr" in self.call_count:
            self.call_count["rasmph_pwr"] += 1
        else:
            self.call_count["rasmph_pwr"] = 1

        result = self._validate(Path(cpx).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        result = self._validate(Path(cc).exists(), result)
        return result

    def data2geotiff(self, DEM_par: str, data: str, type, GeoTIFF: str, nodata):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "data2geotiff", supplied_args))

        if "data2geotiff" in self.call_count:
            self.call_count["data2geotiff"] += 1
        else:
            self.call_count["data2geotiff"] = 1

        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(data).exists(), result)
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(type in valid_values, result)
        Path(GeoTIFF).touch()
        return result

    def disras(self, ras: str, mag, win_sz):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disras", supplied_args))

        if "disras" in self.call_count:
            self.call_count["disras"] += 1
        else:
            self.call_count["disras"] = 1

        result = self._validate(Path(ras).exists(), result)
        return result

    def disrmg(self, unw: str, pwr: str, width, start_unw, start_pwr, nlines, ph_scale, scale, exp, ph_offset):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disrmg", supplied_args))

        if "disrmg" in self.call_count:
            self.call_count["disrmg"] += 1
        else:
            self.call_count["disrmg"] = 1

        result = self._validate(Path(unw).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        return result

    def dis2rmg(self, unw1: str, unw2: str, width1, width2, start, nlines, roff, azoff, ph_scale, ph_offset):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2rmg", supplied_args))

        if "dis2rmg" in self.call_count:
            self.call_count["dis2rmg"] += 1
        else:
            self.call_count["dis2rmg"] = 1

        result = self._validate(Path(unw1).exists(), result)
        result = self._validate(Path(unw2).exists(), result)
        return result

    def dis2ras(self, ras1: str, ras2: str, mag, win_sz):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2ras", supplied_args))

        if "dis2ras" in self.call_count:
            self.call_count["dis2ras"] += 1
        else:
            self.call_count["dis2ras"] = 1

        result = self._validate(Path(ras1).exists(), result)
        result = self._validate(Path(ras2).exists(), result)
        return result

    def dismph_pwr(self, cpx: str, pwr: str, width, start_cpx, start_pwr, nlines, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph_pwr", supplied_args))

        if "dismph_pwr" in self.call_count:
            self.call_count["dismph_pwr"] += 1
        else:
            self.call_count["dismph_pwr"] = 1

        result = self._validate(Path(cpx).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        return result

    def vec_math(self, d1: str, d2: str, d_out: str, width, mode, c1, c2, c3, nflg, nl):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "vec_math", supplied_args))

        if "vec_math" in self.call_count:
            self.call_count["vec_math"] += 1
        else:
            self.call_count["vec_math"] = 1

        result = self._validate(Path(d1).exists(), result)
        result = self._validate(Path(d2).exists(), result)
        Path(d_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6]
        result = self._validate(mode in valid_values, result)
        return result

    def dismph_fft(self, cpx: str, width, start, nlines, scale, exp, nfft, mag, data_type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph_fft", supplied_args))

        if "dismph_fft" in self.call_count:
            self.call_count["dismph_fft"] += 1
        else:
            self.call_count["dismph_fft"] = 1

        result = self._validate(Path(cpx).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def svg_arrow(self, dv_norm: str, dv_phi: str, width, svg: str, image, norm, gridx, gridy, color, thick, head):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "svg_arrow", supplied_args))

        if "svg_arrow" in self.call_count:
            self.call_count["svg_arrow"] += 1
        else:
            self.call_count["svg_arrow"] = 1

        result = self._validate(Path(dv_norm).exists(), result)
        result = self._validate(Path(dv_phi).exists(), result)
        Path(svg).touch()
        return result

    def ras3pwr(self, d1: str, d2: str, d3: str, width, start, nlines, pixavr, pixavaz, scale1, scale2, scale3, exp, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras3pwr", supplied_args))

        if "ras3pwr" in self.call_count:
            self.call_count["ras3pwr"] += 1
        else:
            self.call_count["ras3pwr"] = 1

        result = self._validate(Path(d1).exists(), result)
        result = self._validate(Path(d2).exists(), result)
        result = self._validate(Path(d3).exists(), result)
        Path(rasf).touch()
        return result

    def disdt_pwr(self, data: str, pwr: str, width, start_data, start_pwr, nlines, cycle, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disdt_pwr", supplied_args))

        if "disdt_pwr" in self.call_count:
            self.call_count["disdt_pwr"] += 1
        else:
            self.call_count["disdt_pwr"] = 1

        result = self._validate(Path(data).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        return result

    def real_to_vec(self, cmp1: str, cmp2: str, cmp3: str, width, vec: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "real_to_vec", supplied_args))

        if "real_to_vec" in self.call_count:
            self.call_count["real_to_vec"] += 1
        else:
            self.call_count["real_to_vec"] = 1

        result = self._validate(Path(cmp1).exists(), result)
        result = self._validate(Path(cmp2).exists(), result)
        result = self._validate(Path(cmp3).exists(), result)
        Path(vec).touch()
        return result

    def float2short(self, infile: str, outfile: str, scale, exp, neg, output):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "float2short", supplied_args))

        if "float2short" in self.call_count:
            self.call_count["float2short"] += 1
        else:
            self.call_count["float2short"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        return result

    def disdt_pwr24(self, data: str, pwr: str, width, start_data, start_pwr, nlines, cycle, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disdt_pwr24", supplied_args))

        if "disdt_pwr24" in self.call_count:
            self.call_count["disdt_pwr24"] += 1
        else:
            self.call_count["disdt_pwr24"] = 1

        result = self._validate(Path(data).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        return result

    def rashgt_shd(self, hgt: str, data: str, width, col_post, row_post, start, nlines, pixavr, pixavaz, theta0, phi0, color0, cycle, LR, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rashgt_shd", supplied_args))

        if "rashgt_shd" in self.call_count:
            self.call_count["rashgt_shd"] += 1
        else:
            self.call_count["rashgt_shd"] = 1

        result = self._validate(Path(hgt).exists(), result)
        result = self._validate(Path(data).exists(), result)
        Path(rasf).touch()
        return result

    def ras2ras(self, ras_in: str, ras_out: str, cmap, force24):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras2ras", supplied_args))

        if "ras2ras" in self.call_count:
            self.call_count["ras2ras"] += 1
        else:
            self.call_count["ras2ras"] = 1

        result = self._validate(Path(ras_in).exists(), result)
        Path(ras_out).touch()
        valid_values = [0, 1]
        result = self._validate(force24 in valid_values, result)
        return result

    def ras24_float(self, f1: str, f2: str, f3: str, width, rasf: str, color_model, h0, hrange, imin, imax, sat_min, sat_max, Input, sc1, A1, B1, cyclic1, sc2, A2, B2, sc3, A3, B3, General, start_f1, start_f2, start_f3, nlines, pixavr, pixavaz, LR):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras24_float", supplied_args))

        if "ras24_float" in self.call_count:
            self.call_count["ras24_float"] += 1
        else:
            self.call_count["ras24_float"] = 1

        result = self._validate(Path(f1).exists(), result)
        result = self._validate(Path(f2).exists(), result)
        result = self._validate(Path(f3).exists(), result)
        Path(rasf).touch()
        valid_values = [0, 1, 2]
        result = self._validate(color_model in valid_values, result)
        return result

    def double2float(self, infile: str, outfile: str, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "double2float", supplied_args))

        if "double2float" in self.call_count:
            self.call_count["double2float"] += 1
        else:
            self.call_count["double2float"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        return result

    def float2uchar(self, infile, outfile, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "float2uchar", supplied_args))

        if "float2uchar" in self.call_count:
            self.call_count["float2uchar"] += 1
        else:
            self.call_count["float2uchar"] = 1

        return result

    def float_math(self, d1: str, d2: str, d_out: str, width, mode, roff, loff, nr, nl, c0, zflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "float_math", supplied_args))

        if "float_math" in self.call_count:
            self.call_count["float_math"] += 1
        else:
            self.call_count["float_math"] = 1

        result = self._validate(Path(d1).exists(), result)
        result = self._validate(Path(d2).exists(), result)
        Path(d_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def cpd(self, din: str, dout: str, width, dtype, xoff, nx, yoff, ny):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "cpd", supplied_args))

        if "cpd" in self.call_count:
            self.call_count["cpd"] += 1
        else:
            self.call_count["cpd"] = 1

        result = self._validate(Path(din).exists(), result)
        Path(dout).touch()
        valid_values = [0, 1, 2]
        result = self._validate(dtype in valid_values, result)
        return result

    def disras_dem_par(self, ras: str, DEM_par: str, mag, win_sz):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disras_dem_par", supplied_args))

        if "disras_dem_par" in self.call_count:
            self.call_count["disras_dem_par"] += 1
        else:
            self.call_count["disras_dem_par"] = 1

        result = self._validate(Path(ras).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        return result

    def dismph(self, cpx: str, width, start, nlines, scale, exp, dtype):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph", supplied_args))

        if "dismph" in self.call_count:
            self.call_count["dismph"] += 1
        else:
            self.call_count["dismph"] = 1

        result = self._validate(Path(cpx).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def rasshd(self, DEM: str, width, col_post, row_post, start, nlines, pixavr, pixavaz, theta0, phi0, LR, rasf: str, dtype, zero_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasshd", supplied_args))

        if "rasshd" in self.call_count:
            self.call_count["rasshd"] += 1
        else:
            self.call_count["rasshd"] = 1

        result = self._validate(Path(DEM).exists(), result)
        Path(rasf).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def dis2mph(self, cpx1: str, cpx2: str, width1, width2, start, nlines, roff, azoff, scale, exp, sc_abs1, sc_abs2, dtype):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2mph", supplied_args))

        if "dis2mph" in self.call_count:
            self.call_count["dis2mph"] += 1
        else:
            self.call_count["dis2mph"] = 1

        result = self._validate(Path(cpx1).exists(), result)
        result = self._validate(Path(cpx2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def dis2SLC(self, SLC1: str, SLC2: str, width1, width2, start, nlines, roff, azoff, scale, exp, data_type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2SLC", supplied_args))

        if "dis2SLC" in self.call_count:
            self.call_count["dis2SLC"] += 1
        else:
            self.call_count["dis2SLC"] = 1

        result = self._validate(Path(SLC1).exists(), result)
        result = self._validate(Path(SLC2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def raspwr(self, pwr: str, width, start, nlines, pixavr, pixavaz, scale, exp, LR, rasf: str, data_type, hdrsz):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "raspwr", supplied_args))

        if "raspwr" in self.call_count:
            self.call_count["raspwr"] += 1
        else:
            self.call_count["raspwr"] = 1

        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        valid_values = [0, 1, 2]
        result = self._validate(data_type in valid_values, result)
        return result

    def disSLC(self, SLC: str, width, start, nlines, scale, exp, data_type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disSLC", supplied_args))

        if "disSLC" in self.call_count:
            self.call_count["disSLC"] += 1
        else:
            self.call_count["disSLC"] = 1

        result = self._validate(Path(SLC).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def ras8_colormap(self, model, h0, hrange, ival, sat, cm: str, cm_ras: str, width, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras8_colormap", supplied_args))

        if "ras8_colormap" in self.call_count:
            self.call_count["ras8_colormap"] += 1
        else:
            self.call_count["ras8_colormap"] = 1

        valid_values = [0, 1, 2, 3]
        result = self._validate(model in valid_values, result)
        Path(cm).touch()
        Path(cm_ras).touch()
        return result

    def cp_data(self, infile: str, outfile: str, lbytes, start, nlines, offset, file_ldr, offb, nbyte):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "cp_data", supplied_args))

        if "cp_data" in self.call_count:
            self.call_count["cp_data"] += 1
        else:
            self.call_count["cp_data"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        return result

    def ras_cpt(self, data: str, width, cpt: str, color_model, start, nlines, pixavr, pixavaz, LR, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras_cpt", supplied_args))

        if "ras_cpt" in self.call_count:
            self.call_count["ras_cpt"] += 1
        else:
            self.call_count["ras_cpt"] = 1

        result = self._validate(Path(data).exists(), result)
        result = self._validate(Path(cpt).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(color_model in valid_values, result)
        Path(rasf).touch()
        return result

    def uchar2float(self, data_in: str, data_out, scale, exp, offset):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "uchar2float", supplied_args))

        if "uchar2float" in self.call_count:
            self.call_count["uchar2float"] += 1
        else:
            self.call_count["uchar2float"] = 1

        result = self._validate(Path(data_in).exists(), result)
        return result

    def rasdt_pwr24(self, data: str, pwr: str, width, start_data, start_pwr, nlines, pixavr, pixavaz, cycle, scale, exp, LR, rasf: str, cc, start_cc, cc_min):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasdt_pwr24", supplied_args))

        if "rasdt_pwr24" in self.call_count:
            self.call_count["rasdt_pwr24"] += 1
        else:
            self.call_count["rasdt_pwr24"] = 1

        result = self._validate(Path(data).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        return result

    def float2gbyte(self, infile: str, outfile: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "float2gbyte", supplied_args))

        if "float2gbyte" in self.call_count:
            self.call_count["float2gbyte"] += 1
        else:
            self.call_count["float2gbyte"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        return result

    def disbyte(self, image: str, width, start, nlines, scale):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disbyte", supplied_args))

        if "disbyte" in self.call_count:
            self.call_count["disbyte"] += 1
        else:
            self.call_count["disbyte"] = 1

        result = self._validate(Path(image).exists(), result)
        return result

    def rasrmg(self, unw: str, pwr: str, width, start_unw, start_pwr, nlines, pixavr, pixavaz, ph_scale, scale, exp, ph_offset, LR, rasf: str, cc: str, start_cc, cc_min):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasrmg", supplied_args))

        if "rasrmg" in self.call_count:
            self.call_count["rasrmg"] += 1
        else:
            self.call_count["rasrmg"] = 1

        result = self._validate(Path(unw).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        result = self._validate(Path(cc).exists(), result)
        return result

    def kml_plan(self, MLI_par: str, DEM_par: str, lookup_table: str, kml: str, geoid: str, geoid_par: str, extension, flight_path, t_event, pt_list: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "kml_plan", supplied_args))

        if "kml_plan" in self.call_count:
            self.call_count["kml_plan"] += 1
        else:
            self.call_count["kml_plan"] = 1

        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        result = self._validate(Path(lookup_table).exists(), result)
        Path(kml).touch()
        result = self._validate(Path(geoid).exists(), result)
        result = self._validate(Path(geoid_par).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(flight_path in valid_values, result)
        result = self._validate(Path(pt_list).exists(), result)
        return result

    def kml_map(self, image: str, dem_par: str, kml: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "kml_map", supplied_args))

        if "kml_map" in self.call_count:
            self.call_count["kml_map"] += 1
        else:
            self.call_count["kml_map"] = 1

        result = self._validate(Path(image).exists(), result)
        result = self._validate(Path(dem_par).exists(), result)
        Path(kml).touch()
        return result

    def dismph_pk(self, cpx: str, width, start, nlines, scale, exp, data_type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph_pk", supplied_args))

        if "dismph_pk" in self.call_count:
            self.call_count["dismph_pk"] += 1
        else:
            self.call_count["dismph_pk"] = 1

        result = self._validate(Path(cpx).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def ras_dB(self, pwr: str, width, start, nlines, pixavr, pixavaz, min_dB, max_dB, dB_offset, LR, rasf: str, abs_flag, inverse, channel):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras_dB", supplied_args))

        if "ras_dB" in self.call_count:
            self.call_count["ras_dB"] += 1
        else:
            self.call_count["ras_dB"] = 1

        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        valid_values = [0, 1]
        result = self._validate(abs_flag in valid_values, result)
        valid_values = [1, 2, 3]
        result = self._validate(channel in valid_values, result)
        return result

    def thres_data(self, data_in: str, width, data_out: str, t_data: str, t_min, t_max, dtype):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "thres_data", supplied_args))

        if "thres_data" in self.call_count:
            self.call_count["thres_data"] += 1
        else:
            self.call_count["thres_data"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        result = self._validate(Path(t_data).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(dtype in valid_values, result)
        return result

    def distree(self, flag: str, unw: str, cpx: str, width, start, nlines, ph_scale):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "distree", supplied_args))

        if "distree" in self.call_count:
            self.call_count["distree"] += 1
        else:
            self.call_count["distree"] = 1

        result = self._validate(Path(flag).exists(), result)
        result = self._validate(Path(unw).exists(), result)
        result = self._validate(Path(cpx).exists(), result)
        return result

    def dis_linear(self, pwr: str, width, start, nlines, min, max):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis_linear", supplied_args))

        if "dis_linear" in self.call_count:
            self.call_count["dis_linear"] += 1
        else:
            self.call_count["dis_linear"] = 1

        result = self._validate(Path(pwr).exists(), result)
        return result

    def disdem_par(self, DEM: str, DEM_par: str, start, nlines, exaggerate, theta0, phi0):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disdem_par", supplied_args))

        if "disdem_par" in self.call_count:
            self.call_count["disdem_par"] += 1
        else:
            self.call_count["disdem_par"] = 1

        result = self._validate(Path(DEM).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        return result

    def gcp_2ras(self, ras1: str, ras2: str, gcp: str, mag, win_sz):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "gcp_2ras", supplied_args))

        if "gcp_2ras" in self.call_count:
            self.call_count["gcp_2ras"] += 1
        else:
            self.call_count["gcp_2ras"] = 1

        result = self._validate(Path(ras1).exists(), result)
        result = self._validate(Path(ras2).exists(), result)
        Path(gcp).touch()
        return result

    def rashgt(self, hgt: str, pwr: str, width, start_hgt, start_pwr, nlines, pixavr, pixavaz, m_cycle, scale, exp, LR, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rashgt", supplied_args))

        if "rashgt" in self.call_count:
            self.call_count["rashgt"] += 1
        else:
            self.call_count["rashgt"] = 1

        result = self._validate(Path(hgt).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        return result

    def ras_linear(self, pwr: str, width, start, nlines, pixavr, pixavaz, min, max, LR, rasf: str, inverse, channel):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras_linear", supplied_args))

        if "ras_linear" in self.call_count:
            self.call_count["ras_linear"] += 1
        else:
            self.call_count["ras_linear"] = 1

        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        valid_values = [1, 2, 3]
        result = self._validate(channel in valid_values, result)
        return result

    def tree_edit(self, flag: str, ras: str, mag, win_sz):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "tree_edit", supplied_args))

        if "tree_edit" in self.call_count:
            self.call_count["tree_edit"] += 1
        else:
            self.call_count["tree_edit"] = 1

        result = self._validate(Path(flag).exists(), result)
        result = self._validate(Path(ras).exists(), result)
        return result

    def fill(self, d1: str, d2: str, dout: str, width):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "fill", supplied_args))

        if "fill" in self.call_count:
            self.call_count["fill"] += 1
        else:
            self.call_count["fill"] = 1

        result = self._validate(Path(d1).exists(), result)
        result = self._validate(Path(d2).exists(), result)
        Path(dout).touch()
        return result

    def dis2pwr(self, pwr1: str, pwr2: str, width1, width2, start, nlines, roff, azoff, scale, exp, dtype, sc_abs1, sc_abs2):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2pwr", supplied_args))

        if "dis2pwr" in self.call_count:
            self.call_count["dis2pwr"] += 1
        else:
            self.call_count["dis2pwr"] = 1

        result = self._validate(Path(pwr1).exists(), result)
        result = self._validate(Path(pwr2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def dispwr(self, pwr: str, width, start, nlines, scale, exp, dtype, sc_abs):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dispwr", supplied_args))

        if "dispwr" in self.call_count:
            self.call_count["dispwr"] += 1
        else:
            self.call_count["dispwr"] = 1

        result = self._validate(Path(pwr).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def rasSLC(self, SLC: str, width, start, nlines, pixavr, pixavaz, scale, exp, LR, data_type, hdrsz, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasSLC", supplied_args))

        if "rasSLC" in self.call_count:
            self.call_count["rasSLC"] += 1
        else:
            self.call_count["rasSLC"] = 1

        result = self._validate(Path(SLC).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        Path(rasf).touch()
        return result

    def rasmph_pwr24(self, cpx: str, pwr: str, width, start_cpx, start_pwr, nlines, pixavr, pixavaz, scale, exp, LR, rasf: str, cc, start_cc, cc_min):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasmph_pwr24", supplied_args))

        if "rasmph_pwr24" in self.call_count:
            self.call_count["rasmph_pwr24"] += 1
        else:
            self.call_count["rasmph_pwr24"] = 1

        result = self._validate(Path(cpx).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        return result

    def dismph_ub(self, cpx: str, width, start, nlines, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph_ub", supplied_args))

        if "dismph_ub" in self.call_count:
            self.call_count["dismph_ub"] += 1
        else:
            self.call_count["dismph_ub"] = 1

        result = self._validate(Path(cpx).exists(), result)
        return result

    def short2float(self, infile: str, outfile: str, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "short2float", supplied_args))

        if "short2float" in self.call_count:
            self.call_count["short2float"] += 1
        else:
            self.call_count["short2float"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        return result

    def ras8_color_scale(self, rasf: str, color_model, h0, hrange, ival, sat, chip_width, gap, chip_height, nval):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras8_color_scale", supplied_args))

        if "ras8_color_scale" in self.call_count:
            self.call_count["ras8_color_scale"] += 1
        else:
            self.call_count["ras8_color_scale"] = 1

        Path(rasf).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(color_model in valid_values, result)
        return result

    def rasmph(self, cpx: str, width, start, nlines, pixavr, pixavaz, scale, exp, LR, rasf: str, data_type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasmph", supplied_args))

        if "rasmph" in self.call_count:
            self.call_count["rasmph"] += 1
        else:
            self.call_count["rasmph"] = 1

        result = self._validate(Path(cpx).exists(), result)
        Path(rasf).touch()
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def get_value(self, PAR_in, keyword):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "get_value", supplied_args))

        if "get_value" in self.call_count:
            self.call_count["get_value"] += 1
        else:
            self.call_count["get_value"] = 1

        return result

    def polyras(self, ras: str, mag, win_sz, poly_file: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "polyras", supplied_args))

        if "polyras" in self.call_count:
            self.call_count["polyras"] += 1
        else:
            self.call_count["polyras"] = 1

        result = self._validate(Path(ras).exists(), result)
        Path(poly_file).touch()
        return result

    def rastree(self, flag: str, unw: str, cpx, width, start, nlines, ph_scale, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rastree", supplied_args))

        if "rastree" in self.call_count:
            self.call_count["rastree"] += 1
        else:
            self.call_count["rastree"] = 1

        result = self._validate(Path(flag).exists(), result)
        result = self._validate(Path(unw).exists(), result)
        Path(rasf).touch()
        return result

    def disshd(self, DEM: str, width, col_post, row_post, start, nlines, theta0, phi0, data_type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disshd", supplied_args))

        if "disshd" in self.call_count:
            self.call_count["disshd"] += 1
        else:
            self.call_count["disshd"] = 1

        result = self._validate(Path(DEM).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def dis2gbyte(self, image1: str, image2: str, width1, width2, start, nlines, roff, azoff, scale):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2gbyte", supplied_args))

        if "dis2gbyte" in self.call_count:
            self.call_count["dis2gbyte"] += 1
        else:
            self.call_count["dis2gbyte"] = 1

        result = self._validate(Path(image1).exists(), result)
        result = self._validate(Path(image2).exists(), result)
        return result

    def ras8_float(self, f1: str, f2: str, width, rasf: str, color_model, h0, hrange, imin, imax, sat, Image, sc1, A1, B1, cyclic1, sc2, A2, B2, General, start_f1, start_f2, nlines, pixavr, pixavaz, LR):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras8_float", supplied_args))

        if "ras8_float" in self.call_count:
            self.call_count["ras8_float"] += 1
        else:
            self.call_count["ras8_float"] = 1

        result = self._validate(Path(f1).exists(), result)
        result = self._validate(Path(f2).exists(), result)
        Path(rasf).touch()
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(color_model in valid_values, result)
        return result

    def flip(self, infile: str, outfile: str, width, format, sense):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "flip", supplied_args))

        if "flip" in self.call_count:
            self.call_count["flip"] += 1
        else:
            self.call_count["flip"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        return result

    def create_array(self, output: str, width, nlines, dtype, val, val_im):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "create_array", supplied_args))

        if "create_array" in self.call_count:
            self.call_count["create_array"] += 1
        else:
            self.call_count["create_array"] = 1

        Path(output).touch()
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        return result

    def ascii2float(self, data_in: str, width, data_out: str, loff, nl, coff, nv):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ascii2float", supplied_args))

        if "ascii2float" in self.call_count:
            self.call_count["ascii2float"] += 1
        else:
            self.call_count["ascii2float"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        return result

    def data2tiff(self, data: str, width, type, TIFF: str, nodata, xspacing, yspacing):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "data2tiff", supplied_args))

        if "data2tiff" in self.call_count:
            self.call_count["data2tiff"] += 1
        else:
            self.call_count["data2tiff"] = 1

        result = self._validate(Path(data).exists(), result)
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(type in valid_values, result)
        Path(TIFF).touch()
        return result

    def mapshd(self, DEM: str, width, col_post, row_post, theta0, phi0, shade: str, dtype, zero_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "mapshd", supplied_args))

        if "mapshd" in self.call_count:
            self.call_count["mapshd"] += 1
        else:
            self.call_count["mapshd"] = 1

        result = self._validate(Path(DEM).exists(), result)
        Path(shade).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def kml_pt(self, table: str, lat_col, lon_col, val1_col, val1_label, val2_col, val2_label, val3_col, val3_label, id_col, kml: str, icon_URL, logo_URL, legend_URL, color_model, h0, hrange, imin, imax, sat_min, sat_max, sc1, A1, B1, cyclic1, sc2, A2, B2, sc3, A3, B3):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "kml_pt", supplied_args))

        if "kml_pt" in self.call_count:
            self.call_count["kml_pt"] += 1
        else:
            self.call_count["kml_pt"] = 1

        result = self._validate(Path(table).exists(), result)
        Path(kml).touch()
        valid_values = [0, 1, 2]
        result = self._validate(color_model in valid_values, result)
        return result

    def vec_to_real(self, vec: str, width, index, cmp: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "vec_to_real", supplied_args))

        if "vec_to_real" in self.call_count:
            self.call_count["vec_to_real"] += 1
        else:
            self.call_count["vec_to_real"] = 1

        result = self._validate(Path(vec).exists(), result)
        valid_values = [1, 2, 3]
        result = self._validate(index in valid_values, result)
        Path(cmp).touch()
        return result

    def dismph_pwr24(self, cpx: str, pwr: str, width, start_cpx, start_pwr, nlines, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph_pwr24", supplied_args))

        if "dismph_pwr24" in self.call_count:
            self.call_count["dismph_pwr24"] += 1
        else:
            self.call_count["dismph_pwr24"] = 1

        result = self._validate(Path(cpx).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        return result

    def swap_bytes(self, infile: str, outfile: str, swap_type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "swap_bytes", supplied_args))

        if "swap_bytes" in self.call_count:
            self.call_count["swap_bytes"] += 1
        else:
            self.call_count["swap_bytes"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        valid_values = [2, 4, 8]
        result = self._validate(swap_type in valid_values, result)
        return result

    def float2double(self, infile: str, outfile: str, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "float2double", supplied_args))

        if "float2double" in self.call_count:
            self.call_count["float2double"] += 1
        else:
            self.call_count["float2double"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        return result

    def dis2byte(self, image1: str, image2: str, width1, width2, start, nlines, roff, azoff, scale):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2byte", supplied_args))

        if "dis2byte" in self.call_count:
            self.call_count["dis2byte"] += 1
        else:
            self.call_count["dis2byte"] = 1

        result = self._validate(Path(image1).exists(), result)
        result = self._validate(Path(image2).exists(), result)
        return result

    def svg_poly(self, image: str, dem_par: str, poly: str, svg: str, width, nlines, mode, thick, lcolor, lopac, pcolor, popac, tcolor, font, fsize):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "svg_poly", supplied_args))

        if "svg_poly" in self.call_count:
            self.call_count["svg_poly"] += 1
        else:
            self.call_count["svg_poly"] = 1

        result = self._validate(Path(image).exists(), result)
        result = self._validate(Path(dem_par).exists(), result)
        result = self._validate(Path(poly).exists(), result)
        Path(svg).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        result = self._validate(mode in valid_values, result)
        return result

    def rasbyte(self, raw: str, width, start, nlines, pixavr, pixavaz, scale, LR, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasbyte", supplied_args))

        if "rasbyte" in self.call_count:
            self.call_count["rasbyte"] += 1
        else:
            self.call_count["rasbyte"] = 1

        result = self._validate(Path(raw).exists(), result)
        Path(rasf).touch()
        return result

    def ras_cpt_scale(self, rasf: str, cpt: str, color_model, width, nlines, start_value, end_value):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras_cpt_scale", supplied_args))

        if "ras_cpt_scale" in self.call_count:
            self.call_count["ras_cpt_scale"] += 1
        else:
            self.call_count["ras_cpt_scale"] = 1

        Path(rasf).touch()
        result = self._validate(Path(cpt).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(color_model in valid_values, result)
        return result

    def dis2cc(self, cc1: str, cc2: str, width1, width2, start, nlines, roff, azoff, cmin, cmax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2cc", supplied_args))

        if "dis2cc" in self.call_count:
            self.call_count["dis2cc"] += 1
        else:
            self.call_count["dis2cc"] = 1

        result = self._validate(Path(cc1).exists(), result)
        result = self._validate(Path(cc2).exists(), result)
        return result

    def float2ascii(self, din: str, width, data_out: str, loff, nl):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "float2ascii", supplied_args))

        if "float2ascii" in self.call_count:
            self.call_count["float2ascii"] += 1
        else:
            self.call_count["float2ascii"] = 1

        result = self._validate(Path(din).exists(), result)
        Path(data_out).touch()
        return result

    def replace_values(self, f_in: str, value, new_value, f_out: str, width, rpl_flg, dtype, zflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "replace_values", supplied_args))

        if "replace_values" in self.call_count:
            self.call_count["replace_values"] += 1
        else:
            self.call_count["replace_values"] = 1

        result = self._validate(Path(f_in).exists(), result)
        Path(f_out).touch()
        valid_values = [0, 1, 2]
        result = self._validate(rpl_flg in valid_values, result)
        valid_values = [2, 4]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def sbi_offset(self, sbi_unw: str, SLCf_par: str, SLCb_par: str, OFF_par: str, az_offset: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "sbi_offset", supplied_args))

        if "sbi_offset" in self.call_count:
            self.call_count["sbi_offset"] += 1
        else:
            self.call_count["sbi_offset"] = 1

        result = self._validate(Path(sbi_unw).exists(), result)
        result = self._validate(Path(SLCf_par).exists(), result)
        result = self._validate(Path(SLCb_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(az_offset).touch()
        return result

    def par_ASF_RSAT_SS(self, CEOS_leader: str, CEOS_data: str, GRD_par: str, GRD: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ASF_RSAT_SS", supplied_args))

        if "par_ASF_RSAT_SS" in self.call_count:
            self.call_count["par_ASF_RSAT_SS"] += 1
        else:
            self.call_count["par_ASF_RSAT_SS"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        result = self._validate(Path(CEOS_data).exists(), result)
        Path(GRD_par).touch()
        Path(GRD).touch()
        return result

    def par_ASF_SLC(self, CEOS_leader: str, SLC_par: str, CEOS_data: str, SLC: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ASF_SLC", supplied_args))

        if "par_ASF_SLC" in self.call_count:
            self.call_count["par_ASF_SLC"] += 1
        else:
            self.call_count["par_ASF_SLC"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        Path(SLC_par).touch()
        result = self._validate(Path(CEOS_data).exists(), result)
        Path(SLC).touch()
        return result

    def par_KS_SLC(self, HDF5: str, trunk: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_KS_SLC", supplied_args))

        if "par_KS_SLC" in self.call_count:
            self.call_count["par_KS_SLC"] += 1
        else:
            self.call_count["par_KS_SLC"] = 1

        result = self._validate(Path(HDF5).exists(), result)
        Path(trunk).touch()
        return result

    def ScanSAR_full_aperture_SLC(self, SLC1_tab: str, SLC2_tab, SLCR_tab: str, SLC2_dir, vmode, wflg, imode, order, dtype, n_ovr):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_full_aperture_SLC", supplied_args))

        if "ScanSAR_full_aperture_SLC" in self.call_count:
            self.call_count["ScanSAR_full_aperture_SLC"] += 1
        else:
            self.call_count["ScanSAR_full_aperture_SLC"] = 1

        result = self._validate(Path(SLC1_tab).exists(), result)
        result = self._validate(Path(SLCR_tab).exists(), result)
        valid_values = [0, 1]
        result = self._validate(vmode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(wflg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(imode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def init_offset(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, rlks, azlks, rpos, azpos, offr, offaz, thres, rwin, azwin, cflag, deramp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "init_offset", supplied_args))

        if "init_offset" in self.call_count:
            self.call_count["init_offset"] += 1
        else:
            self.call_count["init_offset"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(cflag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(deramp in valid_values, result)
        return result

    def bridge(self, int: str, flag: str, unw, bridge: str, width, xmin, xmax, ymin, ymax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "bridge", supplied_args))

        if "bridge" in self.call_count:
            self.call_count["bridge"] += 1
        else:
            self.call_count["bridge"] = 1

        result = self._validate(Path(int).exists(), result)
        result = self._validate(Path(flag).exists(), result)
        result = self._validate(Path(bridge).exists(), result)
        return result

    def par_ERSDAC_PALSAR(self, ERSDAC_SLC_par, SLC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ERSDAC_PALSAR", supplied_args))

        if "par_ERSDAC_PALSAR" in self.call_count:
            self.call_count["par_ERSDAC_PALSAR"] += 1
        else:
            self.call_count["par_ERSDAC_PALSAR"] = 1

        Path(SLC_par).touch()
        return result

    def SR_to_GRD(self, MLI_par: str, OFF_par: str, GRD_par, in_file: str, out_file: str, rlks, azlks, interp_mode, grd_rsp, grd_azsp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SR_to_GRD", supplied_args))

        if "SR_to_GRD" in self.call_count:
            self.call_count["SR_to_GRD"] += 1
        else:
            self.call_count["SR_to_GRD"] = 1

        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(in_file).exists(), result)
        Path(out_file).touch()
        valid_values = [0, 1, 2]
        result = self._validate(interp_mode in valid_values, result)
        return result

    def par_ACS_ERS(self, CEOS_SAR_leader: str, SLC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ACS_ERS", supplied_args))

        if "par_ACS_ERS" in self.call_count:
            self.call_count["par_ACS_ERS"] += 1
        else:
            self.call_count["par_ACS_ERS"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(SLC_par).touch()
        return result

    def offset_fit(self, offs: str, ccp: str, OFF_par: str, coffs: str, coffsets: str, thres, npoly, interact_mode):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_fit", supplied_args))

        if "offset_fit" in self.call_count:
            self.call_count["offset_fit"] += 1
        else:
            self.call_count["offset_fit"] = 1

        result = self._validate(Path(offs).exists(), result)
        result = self._validate(Path(ccp).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(coffs).touch()
        Path(coffsets).touch()
        valid_values = [0, 1]
        result = self._validate(interact_mode in valid_values, result)
        return result

    def radcal_PRI(self, PRI: str, PRI_PAR: str, GRD: str, GRD_PAR: str, K_dB, inc_ref, roff, nr, loff, nl):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "radcal_PRI", supplied_args))

        if "radcal_PRI" in self.call_count:
            self.call_count["radcal_PRI"] += 1
        else:
            self.call_count["radcal_PRI"] = 1

        result = self._validate(Path(PRI).exists(), result)
        result = self._validate(Path(PRI_PAR).exists(), result)
        Path(GRD).touch()
        Path(GRD_PAR).touch()
        return result

    def gcp_phase(self, unw: str, OFF_par: str, gcp: str, gcp_ph: str, win_sz):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "gcp_phase", supplied_args))

        if "gcp_phase" in self.call_count:
            self.call_count["gcp_phase"] += 1
        else:
            self.call_count["gcp_phase"] = 1

        result = self._validate(Path(unw).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(gcp).exists(), result)
        Path(gcp_ph).touch()
        return result

    def par_MSP(self, SAR_par: str, PROC_par: str, SLC_MLI_par: str, image_format):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_MSP", supplied_args))

        if "par_MSP" in self.call_count:
            self.call_count["par_MSP"] += 1
        else:
            self.call_count["par_MSP"] = 1

        result = self._validate(Path(SAR_par).exists(), result)
        result = self._validate(Path(PROC_par).exists(), result)
        Path(SLC_MLI_par).touch()
        valid_values = [0, 1, 2]
        result = self._validate(image_format in valid_values, result)
        return result

    def offset_pwr_tracking2(self, SLC1: str, SLC2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, offs: str, ccp: str, OFF_par2: str, offs2: str, rwin, azwin, offsets: str, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, bw_frac, deramp, int_filt, pflag, pltflg, ccs: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_pwr_tracking2", supplied_args))

        if "offset_pwr_tracking2" in self.call_count:
            self.call_count["offset_pwr_tracking2"] += 1
        else:
            self.call_count["offset_pwr_tracking2"] = 1

        result = self._validate(Path(SLC1).exists(), result)
        result = self._validate(Path(SLC2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(offs).touch()
        Path(ccp).touch()
        result = self._validate(Path(OFF_par2).exists(), result)
        result = self._validate(Path(offs2).exists(), result)
        Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(deramp in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(int_filt in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        Path(ccs).touch()
        return result

    def offset_tracking(self, offs: str, ccp: str, SLC_par: str, OFF_par: str, disp_map: str, disp_val: str, mode, thres, poly_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_tracking", supplied_args))

        if "offset_tracking" in self.call_count:
            self.call_count["offset_tracking"] += 1
        else:
            self.call_count["offset_tracking"] = 1

        result = self._validate(Path(offs).exists(), result)
        result = self._validate(Path(ccp).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(disp_map).touch()
        Path(disp_val).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(poly_flag in valid_values, result)
        return result

    def SLC_interp_ScanSAR(self, SLC2_tab: str, SLC2_par: str, SLC1_tab: str, SLC1_par: str, OFF_par: str, SLC2R_tab, SLC_2R: str, SLC2R_par: str, mode, order, SLC2R_dir):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_interp_ScanSAR", supplied_args))

        if "SLC_interp_ScanSAR" in self.call_count:
            self.call_count["SLC_interp_ScanSAR"] += 1
        else:
            self.call_count["SLC_interp_ScanSAR"] = 1

        result = self._validate(Path(SLC2_tab).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(SLC1_tab).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(SLC_2R).touch()
        Path(SLC2R_par).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def error_stat(self, d1: str, d2: str, width, dtype, roff, loff, nr, nl, report):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "error_stat", supplied_args))

        if "error_stat" in self.call_count:
            self.call_count["error_stat"] += 1
        else:
            self.call_count["error_stat"] = 1

        result = self._validate(Path(d1).exists(), result)
        result = self._validate(Path(d2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def res_map(self, hgt: str, gr: str, data: str, SLC_par: str, OFF_par: str, res_hgt: str, res_data: str, nr, naz, azps_res, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "res_map", supplied_args))

        if "res_map" in self.call_count:
            self.call_count["res_map"] += 1
        else:
            self.call_count["res_map"] = 1

        result = self._validate(Path(hgt).exists(), result)
        result = self._validate(Path(gr).exists(), result)
        result = self._validate(Path(data).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(res_hgt).touch()
        Path(res_data).touch()
        return result

    def RSAT2_vec(self, SLC_par: str, RSAT2_orb, nstate):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "RSAT2_vec", supplied_args))

        if "RSAT2_vec" in self.call_count:
            self.call_count["RSAT2_vec"] += 1
        else:
            self.call_count["RSAT2_vec"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        return result

    def par_S1_GRD(self, GeoTIFF: str, annotation_XML: str, calibration_XML: str, noise_XML: str, MLI_par: str, MLI: str, GRD_par: str, GRD: str, eflg, rps, noise_pwr, edge_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_S1_GRD", supplied_args))

        if "par_S1_GRD" in self.call_count:
            self.call_count["par_S1_GRD"] += 1
        else:
            self.call_count["par_S1_GRD"] = 1

        result = self._validate(Path(GeoTIFF).exists(), result)
        result = self._validate(Path(annotation_XML).exists(), result)
        result = self._validate(Path(calibration_XML).exists(), result)
        result = self._validate(Path(noise_XML).exists(), result)
        Path(MLI_par).touch()
        Path(MLI).touch()
        Path(GRD_par).touch()
        Path(GRD).touch()
        return result

    def SLC_deramp_S1_TOPS(self, SLC1_tab: str, SLC2_tab: str, mode, phflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_deramp_S1_TOPS", supplied_args))

        if "SLC_deramp_S1_TOPS" in self.call_count:
            self.call_count["SLC_deramp_S1_TOPS"] += 1
        else:
            self.call_count["SLC_deramp_S1_TOPS"] = 1

        result = self._validate(Path(SLC1_tab).exists(), result)
        result = self._validate(Path(SLC2_tab).exists(), result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(phflg in valid_values, result)
        return result

    def par_RSAT_SGF(self, CEOS_leader: str, CEOS_data: str, GRD_par: str, GRD: str, sc_dB, dt):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RSAT_SGF", supplied_args))

        if "par_RSAT_SGF" in self.call_count:
            self.call_count["par_RSAT_SGF"] += 1
        else:
            self.call_count["par_RSAT_SGF"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        result = self._validate(Path(CEOS_data).exists(), result)
        Path(GRD_par).touch()
        Path(GRD).touch()
        return result

    def par_ICEYE_GRD(self, GeoTIFF: str, XML: str, MLI_par: str, MLI: str, GRD_par: str, GRD: str, rps):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ICEYE_GRD", supplied_args))

        if "par_ICEYE_GRD" in self.call_count:
            self.call_count["par_ICEYE_GRD"] += 1
        else:
            self.call_count["par_ICEYE_GRD"] = 1

        result = self._validate(Path(GeoTIFF).exists(), result)
        result = self._validate(Path(XML).exists(), result)
        Path(MLI_par).touch()
        Path(MLI).touch()
        Path(GRD_par).touch()
        Path(GRD).touch()
        return result

    def par_CS_SLC_TIF(self, GeoTIFF: str, XML: str, trunk: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_CS_SLC_TIF", supplied_args))

        if "par_CS_SLC_TIF" in self.call_count:
            self.call_count["par_CS_SLC_TIF"] += 1
        else:
            self.call_count["par_CS_SLC_TIF"] = 1

        result = self._validate(Path(GeoTIFF).exists(), result)
        result = self._validate(Path(XML).exists(), result)
        Path(trunk).touch()
        return result

    def DORIS_vec(self, SLC_par, DOR: str, nstate):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "DORIS_vec", supplied_args))

        if "DORIS_vec" in self.call_count:
            self.call_count["DORIS_vec"] += 1
        else:
            self.call_count["DORIS_vec"] = 1

        result = self._validate(Path(DOR).exists(), result)
        return result

    def bpf(self, data_in: str, data_out: str, width, fc_x, bw_x, fc_y, bw_y, roff, azoff, nr, naz, data_type, f_mode, beta, fir_len):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "bpf", supplied_args))

        if "bpf" in self.call_count:
            self.call_count["bpf"] += 1
        else:
            self.call_count["bpf"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        return result

    def image_stat(self, image: str, width, roff, loff, nr, nl, report):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "image_stat", supplied_args))

        if "image_stat" in self.call_count:
            self.call_count["image_stat"] += 1
        else:
            self.call_count["image_stat"] = 1

        result = self._validate(Path(image).exists(), result)
        return result

    def clear_flag(self, flag: str, width, flag_bits, Charges, BRANCH, xmin, xmax, ymin, ymax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "clear_flag", supplied_args))

        if "clear_flag" in self.call_count:
            self.call_count["clear_flag"] += 1
        else:
            self.call_count["clear_flag"] = 1

        result = self._validate(Path(flag).exists(), result)
        return result

    def par_NovaSAR_GRD(self, GeoTIFF: str, XML: str, polarization, MLI_par: str, MLI: str, GRD_par: str, GRD: str, rps):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_NovaSAR_GRD", supplied_args))

        if "par_NovaSAR_GRD" in self.call_count:
            self.call_count["par_NovaSAR_GRD"] += 1
        else:
            self.call_count["par_NovaSAR_GRD"] = 1

        result = self._validate(Path(GeoTIFF).exists(), result)
        result = self._validate(Path(XML).exists(), result)
        Path(MLI_par).touch()
        Path(MLI).touch()
        Path(GRD_par).touch()
        Path(GRD).touch()
        return result

    def slant_range(self, SLC_par: str, slr: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "slant_range", supplied_args))

        if "slant_range" in self.call_count:
            self.call_count["slant_range"] += 1
        else:
            self.call_count["slant_range"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        Path(slr).touch()
        return result

    def par_IECAS_SLC(self, aux_data: str, slc_Re: str, slc_Im: str, date: str, SLC_par: str, SLC: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_IECAS_SLC", supplied_args))

        if "par_IECAS_SLC" in self.call_count:
            self.call_count["par_IECAS_SLC"] += 1
        else:
            self.call_count["par_IECAS_SLC"] = 1

        result = self._validate(Path(aux_data).exists(), result)
        result = self._validate(Path(slc_Re).exists(), result)
        result = self._validate(Path(slc_Im).exists(), result)
        result = self._validate(Path(date).exists(), result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        return result

    def corr_flag(self, corr: str, flag, width, corr_thr, xmin, xmax, ymin, ymax, border):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "corr_flag", supplied_args))

        if "corr_flag" in self.call_count:
            self.call_count["corr_flag"] += 1
        else:
            self.call_count["corr_flag"] = 1

        result = self._validate(Path(corr).exists(), result)
        return result

    def ptarg_cal_MLI(self, MLI_par: str, MLI: str, r_samp, az_samp, psigma, c_r_samp, c_az_samp, ptr_image: str, r_plot: str, az_plot: str, pcal: str, osf, win, pltflg, psz, csz, theta_inc):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ptarg_cal_MLI", supplied_args))

        if "ptarg_cal_MLI" in self.call_count:
            self.call_count["ptarg_cal_MLI"] += 1
        else:
            self.call_count["ptarg_cal_MLI"] = 1

        result = self._validate(Path(MLI_par).exists(), result)
        result = self._validate(Path(MLI).exists(), result)
        Path(ptr_image).touch()
        Path(r_plot).touch()
        Path(az_plot).touch()
        Path(pcal).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        return result

    def ORB_filt(self, SLC_par_in: str, SLC_par_out: str, interval, extra):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ORB_filt", supplied_args))

        if "ORB_filt" in self.call_count:
            self.call_count["ORB_filt"] += 1
        else:
            self.call_count["ORB_filt"] = 1

        result = self._validate(Path(SLC_par_in).exists(), result)
        Path(SLC_par_out).touch()
        return result

    def offset_sub(self, offs: str, OFF_par: str, offs_sub: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_sub", supplied_args))

        if "offset_sub" in self.call_count:
            self.call_count["offset_sub"] += 1
        else:
            self.call_count["offset_sub"] = 1

        result = self._validate(Path(offs).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(offs_sub).touch()
        return result

    def radcal_MLI(self, MLI: str, MLI_PAR: str, OFF_par: str, CMLI: str, antenna: str, rloss_flag, ant_flag, refarea_flag, sc_dB, K_dB, pix_area: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "radcal_MLI", supplied_args))

        if "radcal_MLI" in self.call_count:
            self.call_count["radcal_MLI"] += 1
        else:
            self.call_count["radcal_MLI"] = 1

        result = self._validate(Path(MLI).exists(), result)
        result = self._validate(Path(MLI_PAR).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(CMLI).touch()
        result = self._validate(Path(antenna).exists(), result)
        Path(pix_area).touch()
        return result

    def S1_OPOD_vec(self, SLC_par, OPOD: str, nstate):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "S1_OPOD_vec", supplied_args))

        if "S1_OPOD_vec" in self.call_count:
            self.call_count["S1_OPOD_vec"] += 1
        else:
            self.call_count["S1_OPOD_vec"] = 1

        result = self._validate(Path(OPOD).exists(), result)
        return result

    def par_RCM_GRD(self, RCM_dir: str, polarization, radcal, noise, MLI_par: str, MLI: str, GRD_par: str, GRD: str, rps, noise_pwr: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RCM_GRD", supplied_args))

        if "par_RCM_GRD" in self.call_count:
            self.call_count["par_RCM_GRD"] += 1
        else:
            self.call_count["par_RCM_GRD"] = 1

        result = self._validate(Path(RCM_dir).exists(), result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(radcal in valid_values, result)
        Path(MLI_par).touch()
        Path(MLI).touch()
        Path(GRD_par).touch()
        Path(GRD).touch()
        Path(noise_pwr).touch()
        return result

    def offset_pwr_tracking(self, SLC1: str, SLC2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, offs: str, ccp: str, rwin, azwin, offsets: str, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, lanczos, bw_frac, deramp, int_filt, pflag, pltflg, ccs: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_pwr_tracking", supplied_args))

        if "offset_pwr_tracking" in self.call_count:
            self.call_count["offset_pwr_tracking"] += 1
        else:
            self.call_count["offset_pwr_tracking"] = 1

        result = self._validate(Path(SLC1).exists(), result)
        result = self._validate(Path(SLC2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(offs).touch()
        Path(ccp).touch()
        Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(deramp in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(int_filt in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        Path(ccs).touch()
        return result

    def SLC_phase_shift(self, SLC_1: str, SLC_par1: str, SLC_2: str, SLC_par2: str, ph_shift):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_phase_shift", supplied_args))

        if "SLC_phase_shift" in self.call_count:
            self.call_count["SLC_phase_shift"] += 1
        else:
            self.call_count["SLC_phase_shift"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_par1).exists(), result)
        Path(SLC_2).touch()
        Path(SLC_par2).touch()
        return result

    def unw_model(self, interf: str, unw_model: str, unw: str, width, xinit, yinit, ref_ph, width_model):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "unw_model", supplied_args))

        if "unw_model" in self.call_count:
            self.call_count["unw_model"] += 1
        else:
            self.call_count["unw_model"] = 1

        result = self._validate(Path(interf).exists(), result)
        result = self._validate(Path(unw_model).exists(), result)
        Path(unw).touch()
        return result

    def par_TX_SLC(self, annotation_XML: str, COSAR: str, SLC_par: str, SLC: str, pol):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_TX_SLC", supplied_args))

        if "par_TX_SLC" in self.call_count:
            self.call_count["par_TX_SLC"] += 1
        else:
            self.call_count["par_TX_SLC"] = 1

        result = self._validate(Path(annotation_XML).exists(), result)
        result = self._validate(Path(COSAR).exists(), result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        return result

    def ScanSAR_burst_MLI(self, SLC_tab: str, MLI_tab: str, rlks, azlks, bflg, SLCR_tab: str, MLI_dir):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_burst_MLI", supplied_args))

        if "ScanSAR_burst_MLI" in self.call_count:
            self.call_count["ScanSAR_burst_MLI"] += 1
        else:
            self.call_count["ScanSAR_burst_MLI"] = 1

        result = self._validate(Path(SLC_tab).exists(), result)
        Path(MLI_tab).touch()
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        result = self._validate(Path(SLCR_tab).exists(), result)
        return result

    def base_copy(self, SLC1_par: str, baseline_1: str, SLC2_par: str, baseline_2: str, time_rev):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_copy", supplied_args))

        if "base_copy" in self.call_count:
            self.call_count["base_copy"] += 1
        else:
            self.call_count["base_copy"] = 1

        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(baseline_1).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        Path(baseline_2).touch()
        return result

    def SLC_interp_map(self, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, SLC_2R: str, SLC2R_par: str, OFF_par2: str, coffs2_sm: str, loff, nlines, mode, order):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_interp_map", supplied_args))

        if "SLC_interp_map" in self.call_count:
            self.call_count["SLC_interp_map"] += 1
        else:
            self.call_count["SLC_interp_map"] = 1

        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(SLC_2R).touch()
        Path(SLC2R_par).touch()
        result = self._validate(Path(OFF_par2).exists(), result)
        result = self._validate(Path(coffs2_sm).exists(), result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def base_perp(self, baseline: str, SLC1_par: str, OFF_par: str, time_rev):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_perp", supplied_args))

        if "base_perp" in self.call_count:
            self.call_count["base_perp"] += 1
        else:
            self.call_count["base_perp"] = 1

        result = self._validate(Path(baseline).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        return result

    def SLC_adf(self, SLC: str, ref_SLC: str, ref_SLC_par: str, SLC_filt: str, mode, alpha, nfft_r, nfft_az, r_step, az_step, mwin_r, mwin_az):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_adf", supplied_args))

        if "SLC_adf" in self.call_count:
            self.call_count["SLC_adf"] += 1
        else:
            self.call_count["SLC_adf"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(ref_SLC).exists(), result)
        result = self._validate(Path(ref_SLC_par).exists(), result)
        Path(SLC_filt).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(mode in valid_values, result)
        return result

    def offset_SLC(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, offs: str, snr: str, rwin, azwin, offsets: str, n_ovr, nr, naz, thres, ISZ, pflag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_SLC", supplied_args))

        if "offset_SLC" in self.call_count:
            self.call_count["offset_SLC"] += 1
        else:
            self.call_count["offset_SLC"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(offs).touch()
        Path(snr).touch()
        Path(offsets).touch()
        return result

    def adf(self, interf: str, sm: str, cc: str, width, alpha, nfft, cc_win, step, loff, nlines, wfrac):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "adf", supplied_args))

        if "adf" in self.call_count:
            self.call_count["adf"] += 1
        else:
            self.call_count["adf"] = 1

        result = self._validate(Path(interf).exists(), result)
        Path(sm).touch()
        Path(cc).touch()
        return result

    def PRC_vec(self, SLC_par, PRC: str, nstate):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "PRC_vec", supplied_args))

        if "PRC_vec" in self.call_count:
            self.call_count["PRC_vec"] += 1
        else:
            self.call_count["PRC_vec"] = 1

        result = self._validate(Path(PRC).exists(), result)
        return result

    def adapt_filt(self, int: str, sm: str, width, low_snr_thr, filt_width, xmin, xmax, ymin, ymax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "adapt_filt", supplied_args))

        if "adapt_filt" in self.call_count:
            self.call_count["adapt_filt"] += 1
        else:
            self.call_count["adapt_filt"] = 1

        result = self._validate(Path(int).exists(), result)
        Path(sm).touch()
        return result

    def par_RSAT2_SG(self, product_XML: str, lut_XML: str, GeoTIFF: str, polarization: str, GRD_par: str, GRD: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RSAT2_SG", supplied_args))

        if "par_RSAT2_SG" in self.call_count:
            self.call_count["par_RSAT2_SG"] += 1
        else:
            self.call_count["par_RSAT2_SG"] = 1

        result = self._validate(Path(product_XML).exists(), result)
        result = self._validate(Path(lut_XML).exists(), result)
        result = self._validate(Path(GeoTIFF).exists(), result)
        result = self._validate(Path(polarization).exists(), result)
        Path(GRD_par).touch()
        Path(GRD).touch()
        return result

    def ScanSAR_burst_corners(self, SLC_par: str, TOPS_par: str, KML: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_burst_corners", supplied_args))

        if "ScanSAR_burst_corners" in self.call_count:
            self.call_count["ScanSAR_burst_corners"] += 1
        else:
            self.call_count["ScanSAR_burst_corners"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(TOPS_par).exists(), result)
        Path(KML).touch()
        return result

    def par_RCM_SLC_ScanSAR(self, RCM_dir: str, polarization, radcal, noise_in, root_name: str, SLC_tab: str, beam, noise_out):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RCM_SLC_ScanSAR", supplied_args))

        if "par_RCM_SLC_ScanSAR" in self.call_count:
            self.call_count["par_RCM_SLC_ScanSAR"] += 1
        else:
            self.call_count["par_RCM_SLC_ScanSAR"] = 1

        result = self._validate(Path(RCM_dir).exists(), result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(radcal in valid_values, result)
        Path(root_name).touch()
        Path(SLC_tab).touch()
        return result

    def offset_add(self, OFF_par1: str, OFF_par2: str, OFF_par3: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_add", supplied_args))

        if "offset_add" in self.call_count:
            self.call_count["offset_add"] += 1
        else:
            self.call_count["offset_add"] = 1

        result = self._validate(Path(OFF_par1).exists(), result)
        result = self._validate(Path(OFF_par2).exists(), result)
        Path(OFF_par3).touch()
        return result

    def multi_SLC_WSS(self, SLC: str, SLC_par: str, MLI: str, MLI_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_SLC_WSS", supplied_args))

        if "multi_SLC_WSS" in self.call_count:
            self.call_count["multi_SLC_WSS"] += 1
        else:
            self.call_count["multi_SLC_WSS"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        Path(MLI).touch()
        Path(MLI_par).touch()
        return result

    def SLC_cat_ScanSAR(self, SLC_tab1: str, SLC_tab2: str, SLC_tab3: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_cat_ScanSAR", supplied_args))

        if "SLC_cat_ScanSAR" in self.call_count:
            self.call_count["SLC_cat_ScanSAR"] += 1
        else:
            self.call_count["SLC_cat_ScanSAR"] = 1

        result = self._validate(Path(SLC_tab1).exists(), result)
        result = self._validate(Path(SLC_tab2).exists(), result)
        result = self._validate(Path(SLC_tab3).exists(), result)
        return result

    def par_RSAT_SLC(self, CEOS_leader: str, SLC_par: str, CEOS_data: str, SLC: str, sc_dB, dt):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RSAT_SLC", supplied_args))

        if "par_RSAT_SLC" in self.call_count:
            self.call_count["par_RSAT_SLC"] += 1
        else:
            self.call_count["par_RSAT_SLC"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        Path(SLC_par).touch()
        result = self._validate(Path(CEOS_data).exists(), result)
        Path(SLC).touch()
        return result

    def base_est_fft(self, interf: str, SLC1_par: str, OFF_par: str, baseline: str, nazfft, r_samp, az_line, nrfft):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_est_fft", supplied_args))

        if "base_est_fft" in self.call_count:
            self.call_count["base_est_fft"] += 1
        else:
            self.call_count["base_est_fft"] = 1

        result = self._validate(Path(interf).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(baseline).touch()
        return result

    def par_RCM_GRC(self, RCM_dir: str, polarization, radcal, noise, SLC_par: str, SLC: str, GRC_par: str, GRC: str, rps, noise_pwr: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RCM_GRC", supplied_args))

        if "par_RCM_GRC" in self.call_count:
            self.call_count["par_RCM_GRC"] += 1
        else:
            self.call_count["par_RCM_GRC"] = 1

        result = self._validate(Path(RCM_dir).exists(), result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(radcal in valid_values, result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        Path(GRC_par).touch()
        Path(GRC).touch()
        Path(noise_pwr).touch()
        return result

    def base_orbit(self, SLC1_par: str, SLC2_par: str, baseline: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_orbit", supplied_args))

        if "base_orbit" in self.call_count:
            self.call_count["base_orbit"] += 1
        else:
            self.call_count["base_orbit"] = 1

        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        Path(baseline).touch()
        return result

    def par_KS_DGM(self, HDF5: str, trunk: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_KS_DGM", supplied_args))

        if "par_KS_DGM" in self.call_count:
            self.call_count["par_KS_DGM"] += 1
        else:
            self.call_count["par_KS_DGM"] = 1

        result = self._validate(Path(HDF5).exists(), result)
        Path(trunk).touch()
        return result

    def par_RSI_ERS(self, CEOS_SAR_leader, SLC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RSI_ERS", supplied_args))

        if "par_RSI_ERS" in self.call_count:
            self.call_count["par_RSI_ERS"] += 1
        else:
            self.call_count["par_RSI_ERS"] = 1

        Path(SLC_par).touch()
        return result

    def af_SLC(self, SLC_par: str, SLC: str, rwin, azwin, dr, daz, thres, a1_flg, b0_flg, offsets: str, n_ovr, roff, azoff):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "af_SLC", supplied_args))

        if "af_SLC" in self.call_count:
            self.call_count["af_SLC"] += 1
        else:
            self.call_count["af_SLC"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(SLC).exists(), result)
        valid_values = [0, 1]
        result = self._validate(a1_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(b0_flg in valid_values, result)
        Path(offsets).touch()
        return result

    def par_EORC_JERS_SLC(self, CEOS_SAR_leader: str, SLC_par: str, CEOS_data: str, SLC: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_EORC_JERS_SLC", supplied_args))

        if "par_EORC_JERS_SLC" in self.call_count:
            self.call_count["par_EORC_JERS_SLC"] += 1
        else:
            self.call_count["par_EORC_JERS_SLC"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(SLC_par).touch()
        result = self._validate(Path(CEOS_data).exists(), result)
        Path(SLC).touch()
        return result

    def fill_gaps(self, data_in: str, width, data_out: str, dtype, method, max_dist, bp_flag, win, ds_method, ds_size, ds_data: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "fill_gaps", supplied_args))

        if "fill_gaps" in self.call_count:
            self.call_count["fill_gaps"] += 1
        else:
            self.call_count["fill_gaps"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        Path(ds_data).touch()
        return result

    def rascc_mask_thinning(self, ras_in: str, in_file: str, width, ras_out: str, nmax, thresh_1, thresh_nmax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "rascc_mask_thinning", supplied_args))

        if "rascc_mask_thinning" in self.call_count:
            self.call_count["rascc_mask_thinning"] += 1
        else:
            self.call_count["rascc_mask_thinning"] = 1

        result = self._validate(Path(ras_in).exists(), result)
        result = self._validate(Path(in_file).exists(), result)
        Path(ras_out).touch()
        return result

    def GRD_to_SR(self, GRD_par: str, MLI_par: str, OFF_par: str, in_file: str, out_file: str, rlks, azlks, interp_mode, sr_rsp, sr_azsp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "GRD_to_SR", supplied_args))

        if "GRD_to_SR" in self.call_count:
            self.call_count["GRD_to_SR"] += 1
        else:
            self.call_count["GRD_to_SR"] = 1

        result = self._validate(Path(GRD_par).exists(), result)
        Path(MLI_par).touch()
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(in_file).exists(), result)
        Path(out_file).touch()
        valid_values = [0, 1, 2]
        result = self._validate(interp_mode in valid_values, result)
        return result

    def multi_look2(self, SLC: str, SLC_par: str, MLI: str, MLI_par: str, r_dec, az_dec, rwin, azwin, wflg, lanczos, beta):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_look2", supplied_args))

        if "multi_look2" in self.call_count:
            self.call_count["multi_look2"] += 1
        else:
            self.call_count["multi_look2"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        Path(MLI).touch()
        Path(MLI_par).touch()
        valid_values = [0, 1]
        result = self._validate(wflg in valid_values, result)
        return result

    def dcomp_sirc(self, infile: str, outfile: str, samples, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "dcomp_sirc", supplied_args))

        if "dcomp_sirc" in self.call_count:
            self.call_count["dcomp_sirc"] += 1
        else:
            self.call_count["dcomp_sirc"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        return result

    def subtract_phase(self, interf_in: str, phase_file: str, interf_out: str, width, factor):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "subtract_phase", supplied_args))

        if "subtract_phase" in self.call_count:
            self.call_count["subtract_phase"] += 1
        else:
            self.call_count["subtract_phase"] = 1

        result = self._validate(Path(interf_in).exists(), result)
        result = self._validate(Path(phase_file).exists(), result)
        Path(interf_out).touch()
        return result

    def dcomp_sirc_quad(self, infile: str, outfile: str, samples, parameter, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "dcomp_sirc_quad", supplied_args))

        if "dcomp_sirc_quad" in self.call_count:
            self.call_count["dcomp_sirc_quad"] += 1
        else:
            self.call_count["dcomp_sirc_quad"] = 1

        result = self._validate(Path(infile).exists(), result)
        Path(outfile).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        result = self._validate(parameter in valid_values, result)
        return result

    def SLC_copy_ScanSAR(self, SLC1_tab: str, SLC2_tab, BURST_tab: str, dtype, SLC2_dir):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_copy_ScanSAR", supplied_args))

        if "SLC_copy_ScanSAR" in self.call_count:
            self.call_count["SLC_copy_ScanSAR"] += 1
        else:
            self.call_count["SLC_copy_ScanSAR"] = 1

        result = self._validate(Path(SLC1_tab).exists(), result)
        result = self._validate(Path(BURST_tab).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def par_RISAT_SLC(self, CEOS_leader: str, BAND_META: str, SLC_par: str, CEOS_image: str, SLC: str, line_dir, pix_dir, cal_flg, KdB):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RISAT_SLC", supplied_args))

        if "par_RISAT_SLC" in self.call_count:
            self.call_count["par_RISAT_SLC"] += 1
        else:
            self.call_count["par_RISAT_SLC"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        result = self._validate(Path(BAND_META).exists(), result)
        Path(SLC_par).touch()
        result = self._validate(Path(CEOS_image).exists(), result)
        Path(SLC).touch()
        valid_values = [0, 1]
        result = self._validate(cal_flg in valid_values, result)
        return result

    def SLC_freq_shift(self, SLC: str, SLC_par: str, SLC_shift: str, SLC_shift_par: str, freq_shift):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_freq_shift", supplied_args))

        if "SLC_freq_shift" in self.call_count:
            self.call_count["SLC_freq_shift"] += 1
        else:
            self.call_count["SLC_freq_shift"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        Path(SLC_shift).touch()
        Path(SLC_shift_par).touch()
        return result

    def DELFT_vec2(self, SLC_par: str, DELFT_dir, nstate, interval, ODR):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "DELFT_vec2", supplied_args))

        if "DELFT_vec2" in self.call_count:
            self.call_count["DELFT_vec2"] += 1
        else:
            self.call_count["DELFT_vec2"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        return result

    def SLC_intf(self, SLC_1: str, SLC_2R: str, SLC1_par: str, SLC2R_par: str, OFF_par: str, interf: str, rlks, azlks, loff, nlines, sps_flg, azf_flg, rp1_flg, rp2_flg, SLC_1s, SLC_2Rs, SLC1s_par, SLC2Rs_par, az_beta):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_intf", supplied_args))

        if "SLC_intf" in self.call_count:
            self.call_count["SLC_intf"] += 1
        else:
            self.call_count["SLC_intf"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2R).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2R_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(interf).touch()
        valid_values = [1, 0]
        result = self._validate(sps_flg in valid_values, result)
        valid_values = [1, 0]
        result = self._validate(azf_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(rp1_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(rp2_flg in valid_values, result)
        return result

    def ptarg_cal_SLC(self, SLC_par: str, SLC: str, r_samp, az_samp, psigma, c_r_samp, c_az_samp, ptr_image: str, r_plot: str, az_plot: str, pcal: str, osf, win, pltflg, psz, csz, c_image: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ptarg_cal_SLC", supplied_args))

        if "ptarg_cal_SLC" in self.call_count:
            self.call_count["ptarg_cal_SLC"] += 1
        else:
            self.call_count["ptarg_cal_SLC"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(SLC).exists(), result)
        Path(ptr_image).touch()
        Path(r_plot).touch()
        Path(az_plot).touch()
        Path(pcal).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        Path(c_image).touch()
        return result

    def sbi_filt(self, SLC_1: str, SLC1_par: str, SLC2R_par: str, SLCf: str, SLCf_par: str, SLCb: str, SLCb_par: str, norm_sq, iwflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "sbi_filt", supplied_args))

        if "sbi_filt" in self.call_count:
            self.call_count["sbi_filt"] += 1
        else:
            self.call_count["sbi_filt"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2R_par).exists(), result)
        Path(SLCf).touch()
        Path(SLCf_par).touch()
        Path(SLCb).touch()
        Path(SLCb_par).touch()
        valid_values = [0, 1]
        result = self._validate(iwflg in valid_values, result)
        return result

    def ph_slope_base(self, int_in: str, SLC_par: str, OFF_par: str, base: str, int_out: str, int_type, inverse):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ph_slope_base", supplied_args))

        if "ph_slope_base" in self.call_count:
            self.call_count["ph_slope_base"] += 1
        else:
            self.call_count["ph_slope_base"] = 1

        result = self._validate(Path(int_in).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(base).exists(), result)
        Path(int_out).touch()
        return result

    def multi_look_ScanSAR(self, SLC_tab: str, MLI: str, MLI_par: str, rlks, azlks, bflg, SLCR_tab: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_look_ScanSAR", supplied_args))

        if "multi_look_ScanSAR" in self.call_count:
            self.call_count["multi_look_ScanSAR"] += 1
        else:
            self.call_count["multi_look_ScanSAR"] = 1

        result = self._validate(Path(SLC_tab).exists(), result)
        Path(MLI).touch()
        Path(MLI_par).touch()
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        result = self._validate(Path(SLCR_tab).exists(), result)
        return result

    def par_ASNARO2(self, CEOS_data: str, CEOS_leader: str, SLC_par: str, SLC: str, reramp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ASNARO2", supplied_args))

        if "par_ASNARO2" in self.call_count:
            self.call_count["par_ASNARO2"] += 1
        else:
            self.call_count["par_ASNARO2"] = 1

        result = self._validate(Path(CEOS_data).exists(), result)
        result = self._validate(Path(CEOS_leader).exists(), result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        valid_values = [0, 1]
        result = self._validate(reramp in valid_values, result)
        return result

    def grasses(self, int: str, flag: str, unw: str, width, xmin, xmax, ymin, ymax, xinit, yinit, init_ph):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "grasses", supplied_args))

        if "grasses" in self.call_count:
            self.call_count["grasses"] += 1
        else:
            self.call_count["grasses"] = 1

        result = self._validate(Path(int).exists(), result)
        result = self._validate(Path(flag).exists(), result)
        Path(unw).touch()
        return result

    def mask_data(self, data_in: str, width, data_out: str, mask: str, dtype):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "mask_data", supplied_args))

        if "mask_data" in self.call_count:
            self.call_count["mask_data"] += 1
        else:
            self.call_count["mask_data"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        result = self._validate(Path(mask).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def par_ESA_JERS_SEASAT_SLC(self, CEOS_data: str, CEOS_leader: str, SLC_par: str, SLC: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ESA_JERS_SEASAT_SLC", supplied_args))

        if "par_ESA_JERS_SEASAT_SLC" in self.call_count:
            self.call_count["par_ESA_JERS_SEASAT_SLC"] += 1
        else:
            self.call_count["par_ESA_JERS_SEASAT_SLC"] = 1

        result = self._validate(Path(CEOS_data).exists(), result)
        result = self._validate(Path(CEOS_leader).exists(), result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        return result

    def ScanSAR_burst_overlap(self, SLC_tab: str, root_name: str, rlks, azlks, mode, bflg, SLCR_tab: str, dburst):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_burst_overlap", supplied_args))

        if "ScanSAR_burst_overlap" in self.call_count:
            self.call_count["ScanSAR_burst_overlap"] += 1
        else:
            self.call_count["ScanSAR_burst_overlap"] = 1

        result = self._validate(Path(SLC_tab).exists(), result)
        Path(root_name).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        result = self._validate(Path(SLCR_tab).exists(), result)
        return result

    def ORB_prop_SLC(self, SLC_par: str, nstate, interval, extra, mode):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ORB_prop_SLC", supplied_args))

        if "ORB_prop_SLC" in self.call_count:
            self.call_count["ORB_prop_SLC"] += 1
        else:
            self.call_count["ORB_prop_SLC"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        return result

    def interp_ad(self, data_in: str, data_out: str, width, r_max, np_min, np_max, w_mode, dtype, cp_data):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "interp_ad", supplied_args))

        if "interp_ad" in self.call_count:
            self.call_count["interp_ad"] += 1
        else:
            self.call_count["interp_ad"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(w_mode in valid_values, result)
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(cp_data in valid_values, result)
        return result

    def par_RISAT_GRD(self, CEOS_leader: str, BAND_META: str, GRD_par: str, CEOS_image: str, GRD: str, line_dir, pix_dir, cal_flg, KdB):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RISAT_GRD", supplied_args))

        if "par_RISAT_GRD" in self.call_count:
            self.call_count["par_RISAT_GRD"] += 1
        else:
            self.call_count["par_RISAT_GRD"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        result = self._validate(Path(BAND_META).exists(), result)
        Path(GRD_par).touch()
        result = self._validate(Path(CEOS_image).exists(), result)
        Path(GRD).touch()
        valid_values = [0, 1]
        result = self._validate(cal_flg in valid_values, result)
        return result

    def par_RSAT_SCW(self, CEOS_leader: str, CEOS_trailer: str, CEOS_data: str, GRD_par: str, GRD: str, sc_dB, dt):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RSAT_SCW", supplied_args))

        if "par_RSAT_SCW" in self.call_count:
            self.call_count["par_RSAT_SCW"] += 1
        else:
            self.call_count["par_RSAT_SCW"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        result = self._validate(Path(CEOS_trailer).exists(), result)
        result = self._validate(Path(CEOS_data).exists(), result)
        Path(GRD_par).touch()
        Path(GRD).touch()
        return result

    def neutron(self, intensity: str, flag: str, width, n_thres, ymin, ymax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "neutron", supplied_args))

        if "neutron" in self.call_count:
            self.call_count["neutron"] += 1
        else:
            self.call_count["neutron"] = 1

        result = self._validate(Path(intensity).exists(), result)
        result = self._validate(Path(flag).exists(), result)
        return result

    def SLC_mosaic_S1_TOPS(self, SLC_tab: str, SLC: str, SLC_par: str, rlks, azlks, bflg, SLCR_tab: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_mosaic_S1_TOPS", supplied_args))

        if "SLC_mosaic_S1_TOPS" in self.call_count:
            self.call_count["SLC_mosaic_S1_TOPS"] += 1
        else:
            self.call_count["SLC_mosaic_S1_TOPS"] = 1

        result = self._validate(Path(SLC_tab).exists(), result)
        Path(SLC).touch()
        Path(SLC_par).touch()
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        result = self._validate(Path(SLCR_tab).exists(), result)
        return result

    def multi_look(self, SLC: str, SLC_par: str, MLI: str, MLI_par: str, rlks, azlks, loff, nlines, scale, exp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_look", supplied_args))

        if "multi_look" in self.call_count:
            self.call_count["multi_look"] += 1
        else:
            self.call_count["multi_look"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        Path(MLI).touch()
        Path(MLI_par).touch()
        return result

    def mosaic_WB(self, data_tab: str, dtype: str, data_out: str, data_par_out: str, sc_flg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "mosaic_WB", supplied_args))

        if "mosaic_WB" in self.call_count:
            self.call_count["mosaic_WB"] += 1
        else:
            self.call_count["mosaic_WB"] = 1

        result = self._validate(Path(data_tab).exists(), result)
        result = self._validate(Path(dtype).exists(), result)
        Path(data_out).touch()
        Path(data_par_out).touch()
        valid_values = [0, 1]
        result = self._validate(sc_flg in valid_values, result)
        return result

    def ScanSAR_burst_to_mosaic(self, DATA_tab: str, mosaic: str, MLI_par: str, mflg, DATA_tab_ref: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_burst_to_mosaic", supplied_args))

        if "ScanSAR_burst_to_mosaic" in self.call_count:
            self.call_count["ScanSAR_burst_to_mosaic"] += 1
        else:
            self.call_count["ScanSAR_burst_to_mosaic"] = 1

        result = self._validate(Path(DATA_tab).exists(), result)
        Path(mosaic).touch()
        Path(MLI_par).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mflg in valid_values, result)
        result = self._validate(Path(DATA_tab_ref).exists(), result)
        return result

    def par_UAVSAR_SLC(self, ann: str, SLC_MLC_in: str, SLC_MLI_par: str, SLC_MLI_out: str, image_type, image_format, DOP: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_UAVSAR_SLC", supplied_args))

        if "par_UAVSAR_SLC" in self.call_count:
            self.call_count["par_UAVSAR_SLC"] += 1
        else:
            self.call_count["par_UAVSAR_SLC"] = 1

        result = self._validate(Path(ann).exists(), result)
        result = self._validate(Path(SLC_MLC_in).exists(), result)
        Path(SLC_MLI_par).touch()
        Path(SLC_MLI_out).touch()
        valid_values = [0, 2]
        result = self._validate(image_format in valid_values, result)
        result = self._validate(Path(DOP).exists(), result)
        return result

    def ScanSAR_burst_copy(self, SLC: str, SLC_par: str, TOPS_par: str, SLC_out: str, SLC_out_par: str, burst_num, drflg, SLC_par2: str, dtype):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_burst_copy", supplied_args))

        if "ScanSAR_burst_copy" in self.call_count:
            self.call_count["ScanSAR_burst_copy"] += 1
        else:
            self.call_count["ScanSAR_burst_copy"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(TOPS_par).exists(), result)
        Path(SLC_out).touch()
        Path(SLC_out_par).touch()
        valid_values = [0, 1]
        result = self._validate(drflg in valid_values, result)
        Path(SLC_par2).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def hgt_map(self, unw: str, SLC_par: str, OFF_par: str, baseline: str, hgt: str, gr: str, ph_flag, loff, nlines, SLC2R_par):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "hgt_map", supplied_args))

        if "hgt_map" in self.call_count:
            self.call_count["hgt_map"] += 1
        else:
            self.call_count["hgt_map"] = 1

        result = self._validate(Path(unw).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(baseline).exists(), result)
        Path(hgt).touch()
        Path(gr).touch()
        return result

    def par_CS_SLC(self, HDF5: str, trunk: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_CS_SLC", supplied_args))

        if "par_CS_SLC" in self.call_count:
            self.call_count["par_CS_SLC"] += 1
        else:
            self.call_count["par_CS_SLC"] = 1

        result = self._validate(Path(HDF5).exists(), result)
        Path(trunk).touch()
        return result

    def par_TX_GRD(self, annotation_XML: str, GeoTIFF: str, GRD_par, GRD: str, pol):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_TX_GRD", supplied_args))

        if "par_TX_GRD" in self.call_count:
            self.call_count["par_TX_GRD"] += 1
        else:
            self.call_count["par_TX_GRD"] = 1

        result = self._validate(Path(annotation_XML).exists(), result)
        result = self._validate(Path(GeoTIFF).exists(), result)
        Path(GRD).touch()
        return result

    def split_WB(self, data_in: str, data_par_in: str, data_tab: str, dtype: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "split_WB", supplied_args))

        if "split_WB" in self.call_count:
            self.call_count["split_WB"] += 1
        else:
            self.call_count["split_WB"] = 1

        result = self._validate(Path(data_in).exists(), result)
        result = self._validate(Path(data_par_in).exists(), result)
        result = self._validate(Path(data_tab).exists(), result)
        result = self._validate(Path(dtype).exists(), result)
        return result

    def base_init(self, SLC1_par: str, SLC2_par: str, OFF_par: str, interf: str, baseline: str, mflag, nrfft, nazfft, r_samp, az_line):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_init", supplied_args))

        if "base_init" in self.call_count:
            self.call_count["base_init"] += 1
        else:
            self.call_count["base_init"] = 1

        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(interf).exists(), result)
        Path(baseline).touch()
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(mflag in valid_values, result)
        return result

    def par_SIRC(self, CEOS_leader: str, SLC_par: str, UTC_MET):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_SIRC", supplied_args))

        if "par_SIRC" in self.call_count:
            self.call_count["par_SIRC"] += 1
        else:
            self.call_count["par_SIRC"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        Path(SLC_par).touch()
        return result

    def rascc_mask(self, cc: str, pwr: str, width, start_cc, start_pwr, nlines, pixavr, pixavaz, cc_thres, pwr_thres, cc_min, cc_max, scale, exp, LR, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "rascc_mask", supplied_args))

        if "rascc_mask" in self.call_count:
            self.call_count["rascc_mask"] += 1
        else:
            self.call_count["rascc_mask"] = 1

        result = self._validate(Path(cc).exists(), result)
        result = self._validate(Path(pwr).exists(), result)
        Path(rasf).touch()
        return result

    def init_offset_orbit(self, SLC1_par: str, SLC2_par: str, OFF_par, rpos, azpos, cflag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "init_offset_orbit", supplied_args))

        if "init_offset_orbit" in self.call_count:
            self.call_count["init_offset_orbit"] += 1
        else:
            self.call_count["init_offset_orbit"] = 1

        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(cflag in valid_values, result)
        return result

    def interf_SLC(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_pa: str, OFF_par: str, MLI_1: str, MLI_2: str, interf, nrlk, nazlk, loff, nltot, rfilt, azfilt, s_off):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "interf_SLC", supplied_args))

        if "interf_SLC" in self.call_count:
            self.call_count["interf_SLC"] += 1
        else:
            self.call_count["interf_SLC"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_pa).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(MLI_1).touch()
        Path(MLI_2).touch()
        valid_values = [0, 1]
        result = self._validate(rfilt in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(azfilt in valid_values, result)
        return result

    def MLI_cat(self, MLI_1: str, MLI_2: str, MLI1_par: str, MLI2_par: str, MLI_3: str, MLI3_par: str, dtype, mflg, overlap, interp_mode, degree, extrapol):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "MLI_cat", supplied_args))

        if "MLI_cat" in self.call_count:
            self.call_count["MLI_cat"] += 1
        else:
            self.call_count["MLI_cat"] = 1

        result = self._validate(Path(MLI_1).exists(), result)
        result = self._validate(Path(MLI_2).exists(), result)
        result = self._validate(Path(MLI1_par).exists(), result)
        result = self._validate(Path(MLI2_par).exists(), result)
        Path(MLI_3).touch()
        Path(MLI3_par).touch()
        valid_values = [0, 1]
        result = self._validate(mflg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(extrapol in valid_values, result)
        return result

    def par_RCM_SLC(self, RCM_dir: str, polarization, radcal, noise, SLC_par: str, SLC: str, noise_pwr: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RCM_SLC", supplied_args))

        if "par_RCM_SLC" in self.call_count:
            self.call_count["par_RCM_SLC"] += 1
        else:
            self.call_count["par_RCM_SLC"] = 1

        result = self._validate(Path(RCM_dir).exists(), result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(radcal in valid_values, result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        Path(noise_pwr).touch()
        return result

    def phase_slope(self, interf: str, slopes: str, width, win_sz, thres, xmin, xmax, ymin, ymax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "phase_slope", supplied_args))

        if "phase_slope" in self.call_count:
            self.call_count["phase_slope"] += 1
        else:
            self.call_count["phase_slope"] = 1

        result = self._validate(Path(interf).exists(), result)
        Path(slopes).touch()
        return result

    def par_TX_ScanSAR(self, annot_XML: str, swath, SLC_par: str, SLC: str, TOPS_par: str, bwflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_TX_ScanSAR", supplied_args))

        if "par_TX_ScanSAR" in self.call_count:
            self.call_count["par_TX_ScanSAR"] += 1
        else:
            self.call_count["par_TX_ScanSAR"] = 1

        result = self._validate(Path(annot_XML).exists(), result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        Path(TOPS_par).touch()
        valid_values = [0, 1]
        result = self._validate(bwflg in valid_values, result)
        return result

    def ave_image(self, im_list: str, width, ave: str, start, nlines, pixav_x, pixav_y, zflag, nmin):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ave_image", supplied_args))

        if "ave_image" in self.call_count:
            self.call_count["ave_image"] += 1
        else:
            self.call_count["ave_image"] = 1

        result = self._validate(Path(im_list).exists(), result)
        Path(ave).touch()
        valid_values = [0, 1]
        result = self._validate(zflag in valid_values, result)
        return result

    def multi_cpx(self, data_in: str, OFF_par_in: str, data_out: str, OFF_par_out, rlks, azlks, loff, nlines, roff, nsamp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_cpx", supplied_args))

        if "multi_cpx" in self.call_count:
            self.call_count["multi_cpx"] += 1
        else:
            self.call_count["multi_cpx"] = 1

        result = self._validate(Path(data_in).exists(), result)
        result = self._validate(Path(OFF_par_in).exists(), result)
        Path(data_out).touch()
        return result

    def ASAR_LO_phase_drift(self, SLC1_par: str, SLC2_par: str, OFF_par: str, ph_drift: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ASAR_LO_phase_drift", supplied_args))

        if "ASAR_LO_phase_drift" in self.call_count:
            self.call_count["ASAR_LO_phase_drift"] += 1
        else:
            self.call_count["ASAR_LO_phase_drift"] = 1

        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(ph_drift).touch()
        return result

    def radcal_pwr_stat(self, SLC_tab: str, SLC_tab_cal: str, plist: str, MSR_cal, PWR_cal, roff, loff, nr, nl, plist_out):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "radcal_pwr_stat", supplied_args))

        if "radcal_pwr_stat" in self.call_count:
            self.call_count["radcal_pwr_stat"] += 1
        else:
            self.call_count["radcal_pwr_stat"] = 1

        result = self._validate(Path(SLC_tab).exists(), result)
        result = self._validate(Path(SLC_tab_cal).exists(), result)
        result = self._validate(Path(plist).exists(), result)
        return result

    def par_ICEYE_SLC(self, HDF5: str, SLC_par: str, SLC: str, dtype):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ICEYE_SLC", supplied_args))

        if "par_ICEYE_SLC" in self.call_count:
            self.call_count["par_ICEYE_SLC"] += 1
        else:
            self.call_count["par_ICEYE_SLC"] = 1

        result = self._validate(Path(HDF5).exists(), result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def offset_SLC_tracking(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, offs: str, snr: str, rsw, azsw, offsets: str, n_ovr, thres, rstep, azstep, rstart, rstop, azstart, azstop, ISZ, pflag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_SLC_tracking", supplied_args))

        if "offset_SLC_tracking" in self.call_count:
            self.call_count["offset_SLC_tracking"] += 1
        else:
            self.call_count["offset_SLC_tracking"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(offs).touch()
        Path(snr).touch()
        Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        return result

    def tree_cc(self, flag: str, width, mbl, xmin, xmax, ymin, ymax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "tree_cc", supplied_args))

        if "tree_cc" in self.call_count:
            self.call_count["tree_cc"] += 1
        else:
            self.call_count["tree_cc"] = 1

        result = self._validate(Path(flag).exists(), result)
        return result

    def MLI_copy(self, MLI_in: str, MLI_in_par: str, MLI_out: str, MLI_out_par: str, roff, nr, loff, nl):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "MLI_copy", supplied_args))

        if "MLI_copy" in self.call_count:
            self.call_count["MLI_copy"] += 1
        else:
            self.call_count["MLI_copy"] = 1

        result = self._validate(Path(MLI_in).exists(), result)
        result = self._validate(Path(MLI_in_par).exists(), result)
        Path(MLI_out).touch()
        Path(MLI_out_par).touch()
        return result

    def ORRM_vec(self, SLC_par, ORRM: str, nstate):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ORRM_vec", supplied_args))

        if "ORRM_vec" in self.call_count:
            self.call_count["ORRM_vec"] += 1
        else:
            self.call_count["ORRM_vec"] = 1

        result = self._validate(Path(ORRM).exists(), result)
        return result

    def SLC_ovr(self, SLC: str, SLC_par: str, SLC_ovr: str, SLC_ovr_par: str, r_ovr, az_ovr, mode, order):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_ovr", supplied_args))

        if "SLC_ovr" in self.call_count:
            self.call_count["SLC_ovr"] += 1
        else:
            self.call_count["SLC_ovr"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        Path(SLC_ovr).touch()
        Path(SLC_ovr_par).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def tree_gzw(self, flag: str, width, mbl, xmin, xmax, ymin, ymax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "tree_gzw", supplied_args))

        if "tree_gzw" in self.call_count:
            self.call_count["tree_gzw"] += 1
        else:
            self.call_count["tree_gzw"] = 1

        result = self._validate(Path(flag).exists(), result)
        return result

    def mcf(self, interf: str, wgt: str, mask: str, unw: str, width, tri_mode, roff, loff, nr, nlines, npat_r, npat_az, ovrlap, r_init, az_init, init_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "mcf", supplied_args))

        if "mcf" in self.call_count:
            self.call_count["mcf"] += 1
        else:
            self.call_count["mcf"] = 1

        result = self._validate(Path(interf).exists(), result)
        result = self._validate(Path(wgt).exists(), result)
        result = self._validate(Path(mask).exists(), result)
        Path(unw).touch()
        valid_values = [0, 1]
        result = self._validate(tri_mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(init_flag in valid_values, result)
        return result

    def par_ESA_ERS(self, CEOS_SAR_leader: str, SLC_par: str, CEOS_DAT: str, SLC: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ESA_ERS", supplied_args))

        if "par_ESA_ERS" in self.call_count:
            self.call_count["par_ESA_ERS"] += 1
        else:
            self.call_count["par_ESA_ERS"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(SLC_par).touch()
        result = self._validate(Path(CEOS_DAT).exists(), result)
        Path(SLC).touch()
        return result

    def SLC_interp(self, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, SLC_2R: str, SLC2R_par: str, loff, nlines, mode, order):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_interp", supplied_args))

        if "SLC_interp" in self.call_count:
            self.call_count["SLC_interp"] += 1
        else:
            self.call_count["SLC_interp"] = 1

        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(SLC_2R).touch()
        Path(SLC2R_par).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def par_S1_SLC(self, GeoTIFF: str, annotation_XML: str, calibration_XML: str, noise_XML: str, SLC_par: str, SLC: str, TOPS_par: str, dtype, sc_dB, noise_pwr):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_S1_SLC", supplied_args))

        if "par_S1_SLC" in self.call_count:
            self.call_count["par_S1_SLC"] += 1
        else:
            self.call_count["par_S1_SLC"] = 1

        result = self._validate(Path(GeoTIFF).exists(), result)
        result = self._validate(Path(annotation_XML).exists(), result)
        result = self._validate(Path(calibration_XML).exists(), result)
        result = self._validate(Path(noise_XML).exists(), result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        Path(TOPS_par).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def par_ASAR(self, ASAR_ERS_file: str, output_name: str, K_dB, to):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ASAR", supplied_args))

        if "par_ASAR" in self.call_count:
            self.call_count["par_ASAR"] += 1
        else:
            self.call_count["par_ASAR"] = 1

        result = self._validate(Path(ASAR_ERS_file).exists(), result)
        Path(output_name).touch()
        return result

    def par_ASF_96(self, CEOS_SAR_leader: str, SLC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ASF_96", supplied_args))

        if "par_ASF_96" in self.call_count:
            self.call_count["par_ASF_96"] += 1
        else:
            self.call_count["par_ASF_96"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(SLC_par).touch()
        return result

    def ScanSAR_mosaic_to_burst(self, DATA: str, MLI_par, DATA_tab):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_mosaic_to_burst", supplied_args))

        if "ScanSAR_mosaic_to_burst" in self.call_count:
            self.call_count["ScanSAR_mosaic_to_burst"] += 1
        else:
            self.call_count["ScanSAR_mosaic_to_burst"] = 1

        result = self._validate(Path(DATA).exists(), result)
        return result

    def base_ls(self, SLC_par: str, OFF_par: str, gcp_ph: str, baseline: str, ph_flag, bc_flag, bn_flag, bcdot_flag, bndot_flag, bperp_min, SLC2R_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_ls", supplied_args))

        if "base_ls" in self.call_count:
            self.call_count["base_ls"] += 1
        else:
            self.call_count["base_ls"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        result = self._validate(Path(gcp_ph).exists(), result)
        result = self._validate(Path(baseline).exists(), result)
        result = self._validate(Path(SLC2R_par).exists(), result)
        return result

    def az_spec_SLC(self, SLC: str, SLC_par: str, spectrum: str, roff, namb, pltflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "az_spec_SLC", supplied_args))

        if "az_spec_SLC" in self.call_count:
            self.call_count["az_spec_SLC"] += 1
        else:
            self.call_count["az_spec_SLC"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        Path(spectrum).touch()
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def SLC_copy(self, SLC_in: str, SLC_par_in: str, SLC_out: str, SLC_par_out: str, fcase, sc, roff, nr, loff, nl, swap, header_lines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_copy", supplied_args))

        if "SLC_copy" in self.call_count:
            self.call_count["SLC_copy"] += 1
        else:
            self.call_count["SLC_copy"] = 1

        result = self._validate(Path(SLC_in).exists(), result)
        result = self._validate(Path(SLC_par_in).exists(), result)
        Path(SLC_out).touch()
        Path(SLC_par_out).touch()
        valid_values = [1, 2, 3, 4]
        result = self._validate(fcase in valid_values, result)
        valid_values = [0, 1, 2]
        result = self._validate(swap in valid_values, result)
        return result

    def az_integrate(self, data: str, width: str, azi: str, cflg, scale, lz):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "az_integrate", supplied_args))

        if "az_integrate" in self.call_count:
            self.call_count["az_integrate"] += 1
        else:
            self.call_count["az_integrate"] = 1

        result = self._validate(Path(data).exists(), result)
        result = self._validate(Path(width).exists(), result)
        Path(azi).touch()
        valid_values = [0, 1]
        result = self._validate(cflg in valid_values, result)
        return result

    def SLC_cat(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, SLC_3: str, SLC3_par: str, dopflg, iflg, phflg, gainflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_cat", supplied_args))

        if "SLC_cat" in self.call_count:
            self.call_count["SLC_cat"] += 1
        else:
            self.call_count["SLC_cat"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(SLC_3).touch()
        Path(SLC3_par).touch()
        valid_values = [0, 1]
        result = self._validate(dopflg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(iflg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(phflg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(gainflg in valid_values, result)
        return result

    def par_NovaSAR_SLC(self, GeoTIFF: str, XML: str, polarization, SLC_par: str, SLC: str, dtype):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_NovaSAR_SLC", supplied_args))

        if "par_NovaSAR_SLC" in self.call_count:
            self.call_count["par_NovaSAR_SLC"] += 1
        else:
            self.call_count["par_NovaSAR_SLC"] = 1

        result = self._validate(Path(GeoTIFF).exists(), result)
        result = self._validate(Path(XML).exists(), result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def SLC_corners(self, SLC_par: str, terra_alt: str, kml: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_corners", supplied_args))

        if "SLC_corners" in self.call_count:
            self.call_count["SLC_corners"] += 1
        else:
            self.call_count["SLC_corners"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(terra_alt).exists(), result)
        Path(kml).touch()
        return result

    def SLC_deramp(self, SLC_1: str, SLC_par1: str, SLC_2: str, SLC_par2: str, mode, dop_ph: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_deramp", supplied_args))

        if "SLC_deramp" in self.call_count:
            self.call_count["SLC_deramp"] += 1
        else:
            self.call_count["SLC_deramp"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_par1).exists(), result)
        Path(SLC_2).touch()
        Path(SLC_par2).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        Path(dop_ph).touch()
        return result

    def residue(self, int: str, flag: str, width, xmin, xmax, ymin, ymax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "residue", supplied_args))

        if "residue" in self.call_count:
            self.call_count["residue"] += 1
        else:
            self.call_count["residue"] = 1

        result = self._validate(Path(int).exists(), result)
        result = self._validate(Path(flag).exists(), result)
        return result

    def par_PRI(self, CEOS_SAR_leader: str, PRI_par: str, CEOS_DAT: str, PRI: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_PRI", supplied_args))

        if "par_PRI" in self.call_count:
            self.call_count["par_PRI"] += 1
        else:
            self.call_count["par_PRI"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PRI_par).touch()
        result = self._validate(Path(CEOS_DAT).exists(), result)
        Path(PRI).touch()
        return result

    def create_offset(self, SLC1_par: str, SLC2_par: str, OFF_par, algorithm, rlks, azlks, iflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "create_offset", supplied_args))

        if "create_offset" in self.call_count:
            self.call_count["create_offset"] += 1
        else:
            self.call_count["create_offset"] = 1

        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        valid_values = [1, 2]
        result = self._validate(algorithm in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(iflg in valid_values, result)
        return result

    def multi_look_MLI(self, MLI_in: str, MLI_in_par: str, MLI_out: str, MLI_out_par: str, rlks, azlks, loff, nlines, scale, e_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_look_MLI", supplied_args))

        if "multi_look_MLI" in self.call_count:
            self.call_count["multi_look_MLI"] += 1
        else:
            self.call_count["multi_look_MLI"] = 1

        result = self._validate(Path(MLI_in).exists(), result)
        result = self._validate(Path(MLI_in_par).exists(), result)
        Path(MLI_out).touch()
        Path(MLI_out_par).touch()
        valid_values = [0, 1]
        result = self._validate(e_flag in valid_values, result)
        return result

    def multi_real(self, data_in: str, OFF_par_in: str, data_out: str, OFF_par_out, rlks, azlks, loff, nlines, roff, nsamp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_real", supplied_args))

        if "multi_real" in self.call_count:
            self.call_count["multi_real"] += 1
        else:
            self.call_count["multi_real"] = 1

        result = self._validate(Path(data_in).exists(), result)
        result = self._validate(Path(OFF_par_in).exists(), result)
        Path(data_out).touch()
        return result

    def SLC_intf2(self, SLC_1: str, SLC_2R: str, SLC1_par: str, SLC2R_par: str, MLI_1: str, MLI_2R: str, MLI1_par: str, MLI2R_par: str, interf: str, cc: str, r_dec, az_dec, rwin, azwin, wflg, lanczos, beta):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_intf2", supplied_args))

        if "SLC_intf2" in self.call_count:
            self.call_count["SLC_intf2"] += 1
        else:
            self.call_count["SLC_intf2"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2R).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2R_par).exists(), result)
        Path(MLI_1).touch()
        Path(MLI_2R).touch()
        Path(MLI1_par).touch()
        Path(MLI2R_par).touch()
        Path(interf).touch()
        Path(cc).touch()
        valid_values = [0, 1]
        result = self._validate(wflg in valid_values, result)
        return result

    def par_ASF_PRI(self, CEOS_leader: str, CEOS_data: str, GRD_par: str, GRD: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ASF_PRI", supplied_args))

        if "par_ASF_PRI" in self.call_count:
            self.call_count["par_ASF_PRI"] += 1
        else:
            self.call_count["par_ASF_PRI"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        result = self._validate(Path(CEOS_data).exists(), result)
        Path(GRD_par).touch()
        Path(GRD).touch()
        return result

    def offset_pwr(self, SLC1: str, SLC2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, offs: str, ccp: str, rwin, azwin, offsets: str, n_ovr, nr, naz, thres, lanczos, bw_frac, deramp, int_filt, pflag, pltflg, ccs: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_pwr", supplied_args))

        if "offset_pwr" in self.call_count:
            self.call_count["offset_pwr"] += 1
        else:
            self.call_count["offset_pwr"] = 1

        result = self._validate(Path(SLC1).exists(), result)
        result = self._validate(Path(SLC2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(OFF_par).exists(), result)
        Path(offs).touch()
        Path(ccp).touch()
        Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(deramp in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(int_filt in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        Path(ccs).touch()
        return result

    def par_ATLSCI_ERS(self, CEOS_SAR_leader, CEOS_Image: str, SLC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ATLSCI_ERS", supplied_args))

        if "par_ATLSCI_ERS" in self.call_count:
            self.call_count["par_ATLSCI_ERS"] += 1
        else:
            self.call_count["par_ATLSCI_ERS"] = 1

        result = self._validate(Path(CEOS_Image).exists(), result)
        Path(SLC_par).touch()
        return result

    def par_PRI_ESRIN_JERS(self, CEOS_SAR_leader: str, PRI_par: str, CEOS_DAT: str, PRI: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_PRI_ESRIN_JERS", supplied_args))

        if "par_PRI_ESRIN_JERS" in self.call_count:
            self.call_count["par_PRI_ESRIN_JERS"] += 1
        else:
            self.call_count["par_PRI_ESRIN_JERS"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(PRI_par).touch()
        result = self._validate(Path(CEOS_DAT).exists(), result)
        Path(PRI).touch()
        return result

    def par_KC_PALSAR_slr(self, facter_m: str, CEOS_leader: str, SLC_par: str, pol, pls_mode, KC_data: str, pwr: str, fdtab: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_KC_PALSAR_slr", supplied_args))

        if "par_KC_PALSAR_slr" in self.call_count:
            self.call_count["par_KC_PALSAR_slr"] += 1
        else:
            self.call_count["par_KC_PALSAR_slr"] = 1

        result = self._validate(Path(facter_m).exists(), result)
        result = self._validate(Path(CEOS_leader).exists(), result)
        Path(SLC_par).touch()
        valid_values = [1, 2, 3]
        result = self._validate(pls_mode in valid_values, result)
        result = self._validate(Path(KC_data).exists(), result)
        Path(pwr).touch()
        Path(fdtab).touch()
        return result

    def ptarg_SLC(self, SLC_par: str, SLC: str, r_samp, az_samp, ptr_image: str, r_plot: str, az_plot: str, ptr_par: str, osf, win, pltflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ptarg_SLC", supplied_args))

        if "ptarg_SLC" in self.call_count:
            self.call_count["ptarg_SLC"] += 1
        else:
            self.call_count["ptarg_SLC"] = 1

        result = self._validate(Path(SLC_par).exists(), result)
        result = self._validate(Path(SLC).exists(), result)
        Path(ptr_image).touch()
        Path(r_plot).touch()
        Path(az_plot).touch()
        Path(ptr_par).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        return result

    def par_EORC_PALSAR(self, CEOS_leader: str, SLC_par: str, CEOS_data: str, SLC: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_EORC_PALSAR", supplied_args))

        if "par_EORC_PALSAR" in self.call_count:
            self.call_count["par_EORC_PALSAR"] += 1
        else:
            self.call_count["par_EORC_PALSAR"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        Path(SLC_par).touch()
        result = self._validate(Path(CEOS_data).exists(), result)
        Path(SLC).touch()
        return result

    def S1_burstloc(self, annotation_XML: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "S1_burstloc", supplied_args))

        if "S1_burstloc" in self.call_count:
            self.call_count["S1_burstloc"] += 1
        else:
            self.call_count["S1_burstloc"] = 1

        result = self._validate(Path(annotation_XML).exists(), result)
        return result

    def par_GF3_SLC(self, GeoTIFF: str, annotation_XML: str, SLC_par: str, SLC: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_GF3_SLC", supplied_args))

        if "par_GF3_SLC" in self.call_count:
            self.call_count["par_GF3_SLC"] += 1
        else:
            self.call_count["par_GF3_SLC"] = 1

        result = self._validate(Path(GeoTIFF).exists(), result)
        result = self._validate(Path(annotation_XML).exists(), result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        return result

    def cc_wave(self, interf: str, MLI_1: str, MLI_2: str, cc: str, width, bx, by, wflg, xmin, xmax, ymin, ymax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "cc_wave", supplied_args))

        if "cc_wave" in self.call_count:
            self.call_count["cc_wave"] += 1
        else:
            self.call_count["cc_wave"] = 1

        result = self._validate(Path(interf).exists(), result)
        result = self._validate(Path(MLI_1).exists(), result)
        result = self._validate(Path(MLI_2).exists(), result)
        Path(cc).touch()
        return result

    def par_RSAT2_SLC(self, product_XML: str, lut_XML: str, GeoTIFF: str, polarization: str, SLC_par: str, SLC: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RSAT2_SLC", supplied_args))

        if "par_RSAT2_SLC" in self.call_count:
            self.call_count["par_RSAT2_SLC"] += 1
        else:
            self.call_count["par_RSAT2_SLC"] = 1

        result = self._validate(Path(product_XML).exists(), result)
        result = self._validate(Path(lut_XML).exists(), result)
        result = self._validate(Path(GeoTIFF).exists(), result)
        result = self._validate(Path(polarization).exists(), result)
        Path(SLC_par).touch()
        Path(SLC).touch()
        return result

    def residue_cc(self, int: str, flag: str, width, xmin, xmax, ymin, ymax):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "residue_cc", supplied_args))

        if "residue_cc" in self.call_count:
            self.call_count["residue_cc"] += 1
        else:
            self.call_count["residue_cc"] = 1

        result = self._validate(Path(int).exists(), result)
        result = self._validate(Path(flag).exists(), result)
        return result

    def par_PulSAR(self, CEOS_SAR_leader: str, SLC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_PulSAR", supplied_args))

        if "par_PulSAR" in self.call_count:
            self.call_count["par_PulSAR"] += 1
        else:
            self.call_count["par_PulSAR"] = 1

        result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        Path(SLC_par).touch()
        return result

    def par_ASF_91(self, CEOS_leader: str, CEOS_trailer: str, SLC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ASF_91", supplied_args))

        if "par_ASF_91" in self.call_count:
            self.call_count["par_ASF_91"] += 1
        else:
            self.call_count["par_ASF_91"] = 1

        result = self._validate(Path(CEOS_leader).exists(), result)
        result = self._validate(Path(CEOS_trailer).exists(), result)
        Path(SLC_par).touch()
        return result

    def fspf(self, data_in: str, data_out: str, width, dtype, r_max, spf_type, MLI_par):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "fspf", supplied_args))

        if "fspf" in self.call_count:
            self.call_count["fspf"] += 1
        else:
            self.call_count["fspf"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1, 2]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(spf_type in valid_values, result)
        return result

    def radcal_SLC(self, SLC: str, SLC_PAR: str, CSLC: str, CSLC_PAR: str, fcase, antenna, rloss_flag, ant_flag, refarea_flag, sc_dB, K_dB, pix_area: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "radcal_SLC", supplied_args))

        if "radcal_SLC" in self.call_count:
            self.call_count["radcal_SLC"] += 1
        else:
            self.call_count["radcal_SLC"] = 1

        result = self._validate(Path(SLC).exists(), result)
        result = self._validate(Path(SLC_PAR).exists(), result)
        Path(CSLC).touch()
        Path(CSLC_PAR).touch()
        valid_values = [1, 2, 3, 4]
        result = self._validate(fcase in valid_values, result)
        Path(pix_area).touch()
        return result

    def line_interp(self, input, output, width):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "line_interp", supplied_args))

        if "line_interp" in self.call_count:
            self.call_count["line_interp"] += 1
        else:
            self.call_count["line_interp"] = 1

        return result

    def product_cpx(self, f1: str, f2: str, f_out: str, width, start, nlines, conjg_flg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "product_cpx", supplied_args))

        if "product_cpx" in self.call_count:
            self.call_count["product_cpx"] += 1
        else:
            self.call_count["product_cpx"] = 1

        result = self._validate(Path(f1).exists(), result)
        result = self._validate(Path(f2).exists(), result)
        Path(f_out).touch()
        valid_values = [0, 1]
        result = self._validate(conjg_flg in valid_values, result)
        return result

    def ras_majority(self, ras_in: str, ras_out: str, filter_width, LR):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_majority", supplied_args))

        if "ras_majority" in self.call_count:
            self.call_count["ras_majority"] += 1
        else:
            self.call_count["ras_majority"] = 1

        result = self._validate(Path(ras_in).exists(), result)
        Path(ras_out).touch()
        return result

    def average_filter(self, din: str, dout: str, width, bx, by, wflg, min_pt, zflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "average_filter", supplied_args))

        if "average_filter" in self.call_count:
            self.call_count["average_filter"] += 1
        else:
            self.call_count["average_filter"] = 1

        result = self._validate(Path(din).exists(), result)
        Path(dout).touch()
        valid_values = [0, 1, 2]
        result = self._validate(wflg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def mask_class(self, class_map: str, file_in: str, file_out: str, format_flag, LR, selection_flag, n_class, class_1, class_n, null_value):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "mask_class", supplied_args))

        if "mask_class" in self.call_count:
            self.call_count["mask_class"] += 1
        else:
            self.call_count["mask_class"] = 1

        result = self._validate(Path(class_map).exists(), result)
        result = self._validate(Path(file_in).exists(), result)
        Path(file_out).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(format_flag in valid_values, result)
        return result

    def ras_ratio_dB(self, pwr1: str, pwr2: str, width, start_pwr1, start_pwr2, nlines, pixavr, pixavaz, min_value, max_value, dB_offset, LR, abs_flag, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_ratio_dB", supplied_args))

        if "ras_ratio_dB" in self.call_count:
            self.call_count["ras_ratio_dB"] += 1
        else:
            self.call_count["ras_ratio_dB"] = 1

        result = self._validate(Path(pwr1).exists(), result)
        result = self._validate(Path(pwr2).exists(), result)
        Path(rasf).touch()
        return result

    def linear_to_dB(self, data_in: str, data_out: str, width, inverse_flag, null_value):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "linear_to_dB", supplied_args))

        if "linear_to_dB" in self.call_count:
            self.call_count["linear_to_dB"] += 1
        else:
            self.call_count["linear_to_dB"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1]
        result = self._validate(inverse_flag in valid_values, result)
        return result

    def histogram(self, data_in: str, width, polygon: str, hist: str, stat: str, min, max, nbins, mode, lin_log):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "histogram", supplied_args))

        if "histogram" in self.call_count:
            self.call_count["histogram"] += 1
        else:
            self.call_count["histogram"] = 1

        result = self._validate(Path(data_in).exists(), result)
        result = self._validate(Path(polygon).exists(), result)
        Path(hist).touch()
        Path(stat).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(lin_log in valid_values, result)
        return result

    def gamma_map(self, input_data: str, output_data: str, width, nlooks, bx, by):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "gamma_map", supplied_args))

        if "gamma_map" in self.call_count:
            self.call_count["gamma_map"] += 1
        else:
            self.call_count["gamma_map"] = 1

        result = self._validate(Path(input_data).exists(), result)
        Path(output_data).touch()
        return result

    def m_chi(self, s0: str, m: str, s2chi: str, S_par: str, c1: str, c2: str, c3: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "m-chi", supplied_args))

        if "m-chi" in self.call_count:
            self.call_count["m-chi"] += 1
        else:
            self.call_count["m-chi"] = 1

        result = self._validate(Path(s0).exists(), result)
        result = self._validate(Path(m).exists(), result)
        result = self._validate(Path(s2chi).exists(), result)
        result = self._validate(Path(S_par).exists(), result)
        Path(c1).touch()
        Path(c2).touch()
        Path(c3).touch()
        return result

    def reallks(self, image: str, ML_image: str, width, rlks, azlks, start, nlines, r_start, nsamp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "reallks", supplied_args))

        if "reallks" in self.call_count:
            self.call_count["reallks"] += 1
        else:
            self.call_count["reallks"] = 1

        result = self._validate(Path(image).exists(), result)
        Path(ML_image).touch()
        return result

    def frost(self, pwr1: str, pwr1_frost: str, width, fx, sx, power):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "frost", supplied_args))

        if "frost" in self.call_count:
            self.call_count["frost"] += 1
        else:
            self.call_count["frost"] = 1

        result = self._validate(Path(pwr1).exists(), result)
        Path(pwr1_frost).touch()
        return result

    def mt_lee_filt(self, im_list: str, ref_image: str, width, winsz, L_ref, L, cthres, out_list: str, ref_out: str, b_coeff: str, filt_num: str, msr: str, ctr: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "mt_lee_filt", supplied_args))

        if "mt_lee_filt" in self.call_count:
            self.call_count["mt_lee_filt"] += 1
        else:
            self.call_count["mt_lee_filt"] = 1

        result = self._validate(Path(im_list).exists(), result)
        result = self._validate(Path(ref_image).exists(), result)
        result = self._validate(Path(out_list).exists(), result)
        Path(ref_out).touch()
        Path(b_coeff).touch()
        Path(filt_num).touch()
        Path(msr).touch()
        Path(ctr).touch()
        return result

    def multi_stat(self, im_list: str, width, im_out: str, mode, rank, nmin):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "multi_stat", supplied_args))

        if "multi_stat" in self.call_count:
            self.call_count["multi_stat"] += 1
        else:
            self.call_count["multi_stat"] = 1

        result = self._validate(Path(im_list).exists(), result)
        Path(im_out).touch()
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(mode in valid_values, result)
        return result

    def trigo(self, data1: str, func, data2: str, width):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "trigo", supplied_args))

        if "trigo" in self.call_count:
            self.call_count["trigo"] += 1
        else:
            self.call_count["trigo"] = 1

        result = self._validate(Path(data1).exists(), result)
        valid_values = [1]
        result = self._validate(func in valid_values, result)
        Path(data2).touch()
        return result

    def temp_filt(self, data_tab: str, width, waz, wr, wt_flag, zero_flag, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "temp_filt", supplied_args))

        if "temp_filt" in self.call_count:
            self.call_count["temp_filt"] += 1
        else:
            self.call_count["temp_filt"] = 1

        result = self._validate(Path(data_tab).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(wt_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def cpxlks(self, CMPLX: str, ML_CMPLX: str, width, rlks, azlks, start, nlines, r_start, nsamp):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "cpxlks", supplied_args))

        if "cpxlks" in self.call_count:
            self.call_count["cpxlks"] += 1
        else:
            self.call_count["cpxlks"] = 1

        result = self._validate(Path(CMPLX).exists(), result)
        Path(ML_CMPLX).touch()
        return result

    def ras_to_rgb(self, red, green, blue, ras_out: str, LR, null_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_to_rgb", supplied_args))

        if "ras_to_rgb" in self.call_count:
            self.call_count["ras_to_rgb"] += 1
        else:
            self.call_count["ras_to_rgb"] = 1

        Path(ras_out).touch()
        valid_values = [0, 1]
        result = self._validate(null_flag in valid_values, result)
        return result

    def ave2pwr(self, pwr1: str, pwr2: str, pwr_out: str, width, scale_factor):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ave2pwr", supplied_args))

        if "ave2pwr" in self.call_count:
            self.call_count["ave2pwr"] += 1
        else:
            self.call_count["ave2pwr"] = 1

        result = self._validate(Path(pwr1).exists(), result)
        result = self._validate(Path(pwr2).exists(), result)
        Path(pwr_out).touch()
        return result

    def ras_m_chi(self, s1: str, c1: str, c2: str, c3: str, width, start, nlines, pixavr, pixavaz, scale, exp, rasf: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_m-chi", supplied_args))

        if "ras_m-chi" in self.call_count:
            self.call_count["ras_m-chi"] += 1
        else:
            self.call_count["ras_m-chi"] = 1

        result = self._validate(Path(s1).exists(), result)
        result = self._validate(Path(c1).exists(), result)
        result = self._validate(Path(c2).exists(), result)
        result = self._validate(Path(c3).exists(), result)
        Path(rasf).touch()
        return result

    def sigma2gamma(self, pwr1: str, inc: str, gamma: str, width):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "sigma2gamma", supplied_args))

        if "sigma2gamma" in self.call_count:
            self.call_count["sigma2gamma"] += 1
        else:
            self.call_count["sigma2gamma"] = 1

        result = self._validate(Path(pwr1).exists(), result)
        result = self._validate(Path(inc).exists(), result)
        Path(gamma).touch()
        return result

    def hsi_color_scale(self, file_out: str, nval, chip_width, gap, height):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "hsi_color_scale", supplied_args))

        if "hsi_color_scale" in self.call_count:
            self.call_count["hsi_color_scale"] += 1
        else:
            self.call_count["hsi_color_scale"] = 1

        Path(file_out).touch()
        valid_values = [0]
        result = self._validate(nval in valid_values, result)
        return result

    def unw_to_cpx(self, unw: str, cpx: str, width):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "unw_to_cpx", supplied_args))

        if "unw_to_cpx" in self.call_count:
            self.call_count["unw_to_cpx"] += 1
        else:
            self.call_count["unw_to_cpx"] = 1

        result = self._validate(Path(unw).exists(), result)
        Path(cpx).touch()
        return result

    def polyx(self, ):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "polyx", supplied_args))

        if "polyx" in self.call_count:
            self.call_count["polyx"] += 1
        else:
            self.call_count["polyx"] = 1

        return result

    def pauli(self, SLC_HH: str, SLC_VV: str, SLC_HV: str, SLC_HH_par: str, SLC_VV_par: str, SLC_HV_par: str, P: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "pauli", supplied_args))

        if "pauli" in self.call_count:
            self.call_count["pauli"] += 1
        else:
            self.call_count["pauli"] = 1

        result = self._validate(Path(SLC_HH).exists(), result)
        result = self._validate(Path(SLC_VV).exists(), result)
        result = self._validate(Path(SLC_HV).exists(), result)
        result = self._validate(Path(SLC_HH_par).exists(), result)
        result = self._validate(Path(SLC_VV_par).exists(), result)
        result = self._validate(Path(SLC_HV_par).exists(), result)
        Path(P).touch()
        return result

    def m_alpha(self, s0: str, m: str, alpha: str, S_par: str, c1: str, c2: str, c3: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "m-alpha", supplied_args))

        if "m-alpha" in self.call_count:
            self.call_count["m-alpha"] += 1
        else:
            self.call_count["m-alpha"] = 1

        result = self._validate(Path(s0).exists(), result)
        result = self._validate(Path(m).exists(), result)
        result = self._validate(Path(alpha).exists(), result)
        result = self._validate(Path(S_par).exists(), result)
        Path(c1).touch()
        Path(c2).touch()
        Path(c3).touch()
        return result

    def stokes(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, S: str, S_par: str, rlks, azlks, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "stokes", supplied_args))

        if "stokes" in self.call_count:
            self.call_count["stokes"] += 1
        else:
            self.call_count["stokes"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        Path(S).touch()
        Path(S_par).touch()
        return result

    def histogram_ras(self, ras_in: str, polygon, histograms, mean_stdev, percent, lr_flag, start, stop):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "histogram_ras", supplied_args))

        if "histogram_ras" in self.call_count:
            self.call_count["histogram_ras"] += 1
        else:
            self.call_count["histogram_ras"] = 1

        result = self._validate(Path(ras_in).exists(), result)
        return result

    def product(self, data_1: str, data_2: str, product: str, width, bx, by, wgt_flg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "product", supplied_args))

        if "product" in self.call_count:
            self.call_count["product"] += 1
        else:
            self.call_count["product"] = 1

        result = self._validate(Path(data_1).exists(), result)
        result = self._validate(Path(data_2).exists(), result)
        Path(product).touch()
        valid_values = [0, 1, 2]
        result = self._validate(wgt_flg in valid_values, result)
        return result

    def mt_lee_filt_cpx(self, cpx_list: str, ref_image: str, width, winsz, L_ref, cthres, out_list: str, ref_out: str, b_coeff: str, filt_num: str, msr: str, ctr: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "mt_lee_filt_cpx", supplied_args))

        if "mt_lee_filt_cpx" in self.call_count:
            self.call_count["mt_lee_filt_cpx"] += 1
        else:
            self.call_count["mt_lee_filt_cpx"] = 1

        result = self._validate(Path(cpx_list).exists(), result)
        result = self._validate(Path(ref_image).exists(), result)
        result = self._validate(Path(out_list).exists(), result)
        Path(ref_out).touch()
        Path(b_coeff).touch()
        Path(filt_num).touch()
        Path(msr).touch()
        Path(ctr).touch()
        return result

    def polcovar(self, SLC_1: str, SLC_2: str, SLC_3: str, SLC1_par: str, SLC2_par: str, SLC3_par: str, C: str, C_par: str, rlks, azlks, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "polcovar", supplied_args))

        if "polcovar" in self.call_count:
            self.call_count["polcovar"] += 1
        else:
            self.call_count["polcovar"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC_3).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(SLC3_par).exists(), result)
        Path(C).touch()
        Path(C_par).touch()
        return result

    def stokes_qm(self, S: str, S_par: str, m: str, s2chi: str, s2psi: str, m_l: str, m_c: str, lp_ratio: str, cp_ratio: str, mu: str, delta: str, alpha: str, phi: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "stokes_qm", supplied_args))

        if "stokes_qm" in self.call_count:
            self.call_count["stokes_qm"] += 1
        else:
            self.call_count["stokes_qm"] = 1

        result = self._validate(Path(S).exists(), result)
        result = self._validate(Path(S_par).exists(), result)
        Path(m).touch()
        Path(s2chi).touch()
        Path(s2psi).touch()
        Path(m_l).touch()
        Path(m_c).touch()
        Path(lp_ratio).touch()
        Path(cp_ratio).touch()
        Path(mu).touch()
        Path(delta).touch()
        Path(alpha).touch()
        Path(phi).touch()
        return result

    def frame(self, data_in: str, data_out: str, width, dtype, dx1, dx2, dy1, dy2, null_flag, all_flag, null_value, frame_value):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "frame", supplied_args))

        if "frame" in self.call_count:
            self.call_count["frame"] += 1
        else:
            self.call_count["frame"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(null_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(all_flag in valid_values, result)
        return result

    def looks(self, ):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "looks", supplied_args))

        if "looks" in self.call_count:
            self.call_count["looks"] += 1
        else:
            self.call_count["looks"] = 1

        return result

    def polcoh(self, SLC_1: str, SLC_2: str, SLC_3: str, SLC1_par: str, SLC2_par: str, SLC3_par: str, T: str, T_par: str, rlks, azlks, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "polcoh", supplied_args))

        if "polcoh" in self.call_count:
            self.call_count["polcoh"] += 1
        else:
            self.call_count["polcoh"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC_3).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        result = self._validate(Path(SLC3_par).exists(), result)
        Path(T).touch()
        Path(T_par).touch()
        return result

    def lin_comb(self, nfiles, f1: str, f2: str, constant, factor1, factor2, f_out: str, width, start, nlines, pixav_x, pixav_y, zero_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "lin_comb", supplied_args))

        if "lin_comb" in self.call_count:
            self.call_count["lin_comb"] += 1
        else:
            self.call_count["lin_comb"] = 1

        result = self._validate(Path(f1).exists(), result)
        result = self._validate(Path(f2).exists(), result)
        Path(f_out).touch()
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def multi_class_mapping(self, nfiles, f1: str, f2: str, fn: str, classf: str, ras_out: str, width, start, nlines, pixav_x, pixav_y, LR, color_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "multi_class_mapping", supplied_args))

        if "multi_class_mapping" in self.call_count:
            self.call_count["multi_class_mapping"] += 1
        else:
            self.call_count["multi_class_mapping"] = 1

        result = self._validate(Path(f1).exists(), result)
        result = self._validate(Path(f2).exists(), result)
        result = self._validate(Path(fn).exists(), result)
        result = self._validate(Path(classf).exists(), result)
        Path(ras_out).touch()
        valid_values = [0, 1]
        result = self._validate(color_flag in valid_values, result)
        return result

    def takecut(self, data_in: str, width, report: str, mode, pos, pr_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "takecut", supplied_args))

        if "takecut" in self.call_count:
            self.call_count["takecut"] += 1
        else:
            self.call_count["takecut"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(report).touch()
        valid_values = [0, 1]
        result = self._validate(pr_flag in valid_values, result)
        return result

    def polyx_phase(self, data: str, width, polygon: str, report: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "polyx_phase", supplied_args))

        if "polyx_phase" in self.call_count:
            self.call_count["polyx_phase"] += 1
        else:
            self.call_count["polyx_phase"] = 1

        result = self._validate(Path(data).exists(), result)
        result = self._validate(Path(polygon).exists(), result)
        Path(report).touch()
        return result

    def temp_lin_var(self, data_tab: str, mean: str, stdev: str, width, waz, wr, wt_flag, zero_flag, loff, nlines, norm_pow):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "temp_lin_var", supplied_args))

        if "temp_lin_var" in self.call_count:
            self.call_count["temp_lin_var"] += 1
        else:
            self.call_count["temp_lin_var"] = 1

        result = self._validate(Path(data_tab).exists(), result)
        Path(mean).touch()
        Path(stdev).touch()
        valid_values = [0, 1, 2]
        result = self._validate(wt_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def ave_cpx(self, cpx_list: str, width, ave: str, start, nlines, zflag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ave_cpx", supplied_args))

        if "ave_cpx" in self.call_count:
            self.call_count["ave_cpx"] += 1
        else:
            self.call_count["ave_cpx"] = 1

        result = self._validate(Path(cpx_list).exists(), result)
        Path(ave).touch()
        valid_values = [0, 1]
        result = self._validate(zflag in valid_values, result)
        return result

    def mask_op(self, mask_1: str, mask_2: str, mask_out: str, mode):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "mask_op", supplied_args))

        if "mask_op" in self.call_count:
            self.call_count["mask_op"] += 1
        else:
            self.call_count["mask_op"] = 1

        result = self._validate(Path(mask_1).exists(), result)
        result = self._validate(Path(mask_2).exists(), result)
        Path(mask_out).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        return result

    def quad2cp(self, SLC_HH: str, SLC_HV: str, SLC_VH: str, SLC_VV: str, SLC_HH_par: str, SLC_HV_par: str, SLC_VH_par: str, SLC_VV_par: str, CP: str, TX_pol):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "quad2cp", supplied_args))

        if "quad2cp" in self.call_count:
            self.call_count["quad2cp"] += 1
        else:
            self.call_count["quad2cp"] = 1

        result = self._validate(Path(SLC_HH).exists(), result)
        result = self._validate(Path(SLC_HV).exists(), result)
        result = self._validate(Path(SLC_VH).exists(), result)
        result = self._validate(Path(SLC_VV).exists(), result)
        result = self._validate(Path(SLC_HH_par).exists(), result)
        result = self._validate(Path(SLC_HV_par).exists(), result)
        result = self._validate(Path(SLC_VH_par).exists(), result)
        result = self._validate(Path(SLC_VV_par).exists(), result)
        Path(CP).touch()
        valid_values = [0, 1]
        result = self._validate(TX_pol in valid_values, result)
        return result

    def poly_math(self, ):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "poly_math", supplied_args))

        if "poly_math" in self.call_count:
            self.call_count["poly_math"] += 1
        else:
            self.call_count["poly_math"] = 1

        return result

    def single_class_mapping(self, nfiles, f1: str, lt1, ut1, fn: str, ltn, utn, ras_out: str, width, start, nlines, pixav_x, pixav_y, LR):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "single_class_mapping", supplied_args))

        if "single_class_mapping" in self.call_count:
            self.call_count["single_class_mapping"] += 1
        else:
            self.call_count["single_class_mapping"] = 1

        result = self._validate(Path(f1).exists(), result)
        result = self._validate(Path(fn).exists(), result)
        Path(ras_out).touch()
        return result

    def drawthat(self, ras_in: str, ras_out: str, pt_list: str, mode, r, g, b, xs, zflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "drawthat", supplied_args))

        if "drawthat" in self.call_count:
            self.call_count["drawthat"] += 1
        else:
            self.call_count["drawthat"] = 1

        result = self._validate(Path(ras_in).exists(), result)
        Path(ras_out).touch()
        result = self._validate(Path(pt_list).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def m_delta(self, s0: str, m: str, delta: str, S_par: str, c1: str, c2: str, c3: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "m-delta", supplied_args))

        if "m-delta" in self.call_count:
            self.call_count["m-delta"] += 1
        else:
            self.call_count["m-delta"] = 1

        result = self._validate(Path(s0).exists(), result)
        result = self._validate(Path(m).exists(), result)
        result = self._validate(Path(delta).exists(), result)
        result = self._validate(Path(S_par).exists(), result)
        Path(c1).touch()
        Path(c2).touch()
        Path(c3).touch()
        return result

    def temp_filt_ad(self, data_tab: str, width, zero_flag, loffset, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "temp_filt_ad", supplied_args))

        if "temp_filt_ad" in self.call_count:
            self.call_count["temp_filt_ad"] += 1
        else:
            self.call_count["temp_filt_ad"] = 1

        result = self._validate(Path(data_tab).exists(), result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def bm3d(self, data_in: str, width, data_out: str, dtype, profile, looks, sigma, block_size, s_dist, step, d_max, t1d):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "bm3d", supplied_args))

        if "bm3d" in self.call_count:
            self.call_count["bm3d"] += 1
        else:
            self.call_count["bm3d"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1, 2]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1, 2, 3, 4, 5, 6]
        result = self._validate(profile in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(d_max in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(t1d in valid_values, result)
        return result

    def restore_float(self, input, output, width, interp_limit):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "restore_float", supplied_args))

        if "restore_float" in self.call_count:
            self.call_count["restore_float"] += 1
        else:
            self.call_count["restore_float"] = 1

        return result

    def haalpha(self, alpha: str, beta: str, gamma: str, SLC_par: str, anisotropy: str, entropy: str, lambda1: str, lambda2: str, lambda3: str, MLI_par: str, rlks, azlks, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "haalpha", supplied_args))

        if "haalpha" in self.call_count:
            self.call_count["haalpha"] += 1
        else:
            self.call_count["haalpha"] = 1

        Path(alpha).touch()
        result = self._validate(Path(beta).exists(), result)
        result = self._validate(Path(gamma).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        Path(anisotropy).touch()
        Path(entropy).touch()
        Path(lambda1).touch()
        Path(lambda2).touch()
        Path(lambda3).touch()
        Path(MLI_par).touch()
        return result

    def wolf(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, J: str, J_par: str, rlks, azlks, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "wolf", supplied_args))

        if "wolf" in self.call_count:
            self.call_count["wolf"] += 1
        else:
            self.call_count["wolf"] = 1

        result = self._validate(Path(SLC_1).exists(), result)
        result = self._validate(Path(SLC_2).exists(), result)
        result = self._validate(Path(SLC1_par).exists(), result)
        result = self._validate(Path(SLC2_par).exists(), result)
        Path(J).touch()
        Path(J_par).touch()
        return result

    def takethat_dem_par(self, data_in: str, width, positions: str, DEM_par: str, report: str, mode, zero_flag, nn_flag, print_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "takethat_dem_par", supplied_args))

        if "takethat_dem_par" in self.call_count:
            self.call_count["takethat_dem_par"] += 1
        else:
            self.call_count["takethat_dem_par"] = 1

        result = self._validate(Path(data_in).exists(), result)
        result = self._validate(Path(positions).exists(), result)
        result = self._validate(Path(DEM_par).exists(), result)
        Path(report).touch()
        return result

    def lee(self, input_data: str, output_data: str, width, nlooks, bx, by):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "lee", supplied_args))

        if "lee" in self.call_count:
            self.call_count["lee"] += 1
        else:
            self.call_count["lee"] = 1

        result = self._validate(Path(input_data).exists(), result)
        Path(output_data).touch()
        return result

    def enh_lee(self, input_data: str, output_data: str, width, nlooks, damp, bx, by):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "enh_lee", supplied_args))

        if "enh_lee" in self.call_count:
            self.call_count["enh_lee"] += 1
        else:
            self.call_count["enh_lee"] = 1

        result = self._validate(Path(input_data).exists(), result)
        Path(output_data).touch()
        return result

    def ratio(self, d1: str, d2: str, ratio: str, width, bx, by, wgt_flg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ratio", supplied_args))

        if "ratio" in self.call_count:
            self.call_count["ratio"] += 1
        else:
            self.call_count["ratio"] = 1

        result = self._validate(Path(d1).exists(), result)
        result = self._validate(Path(d2).exists(), result)
        Path(ratio).touch()
        valid_values = [0, 1, 2]
        result = self._validate(wgt_flg in valid_values, result)
        return result

    def cc_ad(self, interf: str, pwr1: str, pwr2: str, slope: str, texture: str, cc_ad: str, width, box_min, box_max, wgt_flag, loff, nl):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "cc_ad", supplied_args))

        if "cc_ad" in self.call_count:
            self.call_count["cc_ad"] += 1
        else:
            self.call_count["cc_ad"] = 1

        result = self._validate(Path(interf).exists(), result)
        result = self._validate(Path(pwr1).exists(), result)
        result = self._validate(Path(pwr2).exists(), result)
        result = self._validate(Path(slope).exists(), result)
        result = self._validate(Path(texture).exists(), result)
        Path(cc_ad).touch()
        valid_values = [0, 1]
        result = self._validate(wgt_flag in valid_values, result)
        return result

    def ras_ras(self, ras_in: str, ras_out: str, col_looks, row_looks, LR, r_lin_log, g_lin_log, b_lin_log, force24):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_ras", supplied_args))

        if "ras_ras" in self.call_count:
            self.call_count["ras_ras"] += 1
        else:
            self.call_count["ras_ras"] = 1

        result = self._validate(Path(ras_in).exists(), result)
        Path(ras_out).touch()
        valid_values = [0, 1]
        result = self._validate(force24 in valid_values, result)
        return result

    def lin_comb_cpx(self, nfiles, f1: str, f2: str, constant_r, constant_i, factor1_r, factor1_i, factor2_r, factor2_i, f_out: str, width, start, nlines, pixav_x, pixav_y, zero_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "lin_comb_cpx", supplied_args))

        if "lin_comb_cpx" in self.call_count:
            self.call_count["lin_comb_cpx"] += 1
        else:
            self.call_count["lin_comb_cpx"] = 1

        result = self._validate(Path(f1).exists(), result)
        result = self._validate(Path(f2).exists(), result)
        Path(f_out).touch()
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def median_filter(self, din: str, dout: str, width, bx, by, min_pt, zflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "median_filter", supplied_args))

        if "median_filter" in self.call_count:
            self.call_count["median_filter"] += 1
        else:
            self.call_count["median_filter"] = 1

        result = self._validate(Path(din).exists(), result)
        Path(dout).touch()
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def takethat(self, data_in: str, width, positions: str, report: str, mode, zero_flag, nn_flag, print_flag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "takethat", supplied_args))

        if "takethat" in self.call_count:
            self.call_count["takethat"] += 1
        else:
            self.call_count["takethat"] = 1

        result = self._validate(Path(data_in).exists(), result)
        result = self._validate(Path(positions).exists(), result)
        Path(report).touch()
        return result

    def cc_monitoring(self, nfiles, f1: str, f2: str, ras_out: str, width, cc_thresh, start, nlines, pixav_x, pixav_y, LR):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "cc_monitoring", supplied_args))

        if "cc_monitoring" in self.call_count:
            self.call_count["cc_monitoring"] += 1
        else:
            self.call_count["cc_monitoring"] = 1

        result = self._validate(Path(f1).exists(), result)
        result = self._validate(Path(f2).exists(), result)
        Path(ras_out).touch()
        return result

    def temp_log_var(self, data_tab: str, mean: str, stdev: str, width, waz, wr, wt_flag, zero_flag, loff, nlines):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "temp_log_var", supplied_args))

        if "temp_log_var" in self.call_count:
            self.call_count["temp_log_var"] += 1
        else:
            self.call_count["temp_log_var"] = 1

        result = self._validate(Path(data_tab).exists(), result)
        Path(mean).touch()
        Path(stdev).touch()
        valid_values = [0, 1, 2]
        result = self._validate(wt_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def ras_to_hsi(self, HUE: str, SATURATION: str, INTENSITY: str, ras_out: str, LR, cflg):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_to_hsi", supplied_args))

        if "ras_to_hsi" in self.call_count:
            self.call_count["ras_to_hsi"] += 1
        else:
            self.call_count["ras_to_hsi"] = 1

        result = self._validate(Path(HUE).exists(), result)
        result = self._validate(Path(SATURATION).exists(), result)
        result = self._validate(Path(INTENSITY).exists(), result)
        Path(ras_out).touch()
        valid_values = [0, 1]
        result = self._validate(cflg in valid_values, result)
        return result

    def edge_detection(self, data_in: str, width, data_out: str, dtype, op_flg, sigma_x, sigma_y, T1, T2, min_seg_size, max_reg_len, max_reg_std, max_reg_dist, seg_out: str, line_filt, max_line_std):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "edge_detection", supplied_args))

        if "edge_detection" in self.call_count:
            self.call_count["edge_detection"] += 1
        else:
            self.call_count["edge_detection"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(data_out).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1, 2]
        result = self._validate(op_flg in valid_values, result)
        Path(seg_out).touch()
        valid_values = [0, 1]
        result = self._validate(line_filt in valid_values, result)
        return result

    def texture(self, data_in: str, format_flag, texture: str, width, type, bx, by, r_looks, az_looks, wgt_flag, data_in_mean: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "texture", supplied_args))

        if "texture" in self.call_count:
            self.call_count["texture"] += 1
        else:
            self.call_count["texture"] = 1

        result = self._validate(Path(data_in).exists(), result)
        Path(texture).touch()
        valid_values = [0, 1]
        result = self._validate(type in valid_values, result)
        valid_values = [0, 1, 2]
        result = self._validate(wgt_flag in valid_values, result)
        result = self._validate(Path(data_in_mean).exists(), result)
        return result

    def diplane_helix(self, LL: str, RR: str, SLC_par: str, diplane: str, helix: str, MLI_par: str, rlks, azlks, loff, nlines, scale):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "diplane_helix", supplied_args))

        if "diplane_helix" in self.call_count:
            self.call_count["diplane_helix"] += 1
        else:
            self.call_count["diplane_helix"] = 1

        result = self._validate(Path(LL).exists(), result)
        result = self._validate(Path(RR).exists(), result)
        result = self._validate(Path(SLC_par).exists(), result)
        Path(diplane).touch()
        Path(helix).touch()
        Path(MLI_par).touch()
        return result

    def lin_comb_ref(self, f1: str, f2: str, constant, factor1, factor2, f_out: str, width, roff, loff, nr, nl, zflag):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "lin_comb_ref", supplied_args))

        if "lin_comb_ref" in self.call_count:
            self.call_count["lin_comb_ref"] += 1
        else:
            self.call_count["lin_comb_ref"] = 1

        result = self._validate(Path(f1).exists(), result)
        result = self._validate(Path(f2).exists(), result)
        Path(f_out).touch()
        valid_values = [0, 1]
        result = self._validate(zflag in valid_values, result)
        return result
