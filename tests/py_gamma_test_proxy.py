from pathlib import Path
from typing import Sequence, NamedTuple, Dict, Union

PyGammaCall = NamedTuple("PyGammaCall", [("module", str), ("program", str), ("parameters", Dict[str, object])])


class SimpleParFile(object):
    values = {}

    def __init__(self, path):
        with open(path, 'r') as file:
            lines = file.read().splitlines()[2:]  # Skip header lines

            for line in lines:
                value_id = line.split(':')[0]
                value_data = line[len(value_id)+2:].strip()

                self.values[value_id] = value_data

    def get_value(self, value_id: str, dtype = str, index: int = 0):
        if dtype == str:
            return self.values[value_id]

        return dtype(self.values[value_id].split()[index])


class PyGammaTestProxy(object):
    ParFile = SimpleParFile

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

    def gc_map_fine(self, gc_in: str, width, DIFF_par, gc_out: str, ref_flg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map_fine", supplied_args))

        if "gc_map_fine" in self.call_count:
            self.call_count["gc_map_fine"] += 1
        else:
            self.call_count["gc_map_fine"] = 1

        if gc_in is not None:
            result = self._validate(Path(gc_in).exists(), result)
        if gc_out is not None:
            Path(gc_out).touch()
        valid_values = [0, 1]
        result = self._validate(ref_flg in valid_values, result)
        return result

    def diff_ls_fit(self, unw_1: str, unw_2: str, DIFF_par: str, nr = None, naz = None, mask = None, plot_data: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "diff_ls_fit", supplied_args))

        if "diff_ls_fit" in self.call_count:
            self.call_count["diff_ls_fit"] += 1
        else:
            self.call_count["diff_ls_fit"] = 1

        if unw_1 is not None:
            result = self._validate(Path(unw_1).exists(), result)
        if unw_2 is not None:
            result = self._validate(Path(unw_2).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if plot_data is not None:
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

        if WSS_tab is not None:
            result = self._validate(Path(WSS_tab).exists(), result)
        if MLI_par is not None:
            Path(MLI_par).touch()
        if WSS_data is not None:
            Path(WSS_data).touch()
        valid_values = [0, 1]
        result = self._validate(type in valid_values, result)
        return result

    def dispmap_vec_offset(self, DEM_par: str, DEM: str, dispmap_r: str, dispmap_az: str, lv_theta: str, lv_phi: str, dv_norm: str, dv_theta: str = None, dv_phi: str = None, dv_x: str = None, dv_y: str = None, dv_z: str = None, mask_angle = None, mode = None, ax_north = None, ax_east = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_vec_offset", supplied_args))

        if "dispmap_vec_offset" in self.call_count:
            self.call_count["dispmap_vec_offset"] += 1
        else:
            self.call_count["dispmap_vec_offset"] = 1

        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if dispmap_r is not None:
            result = self._validate(Path(dispmap_r).exists(), result)
        if dispmap_az is not None:
            result = self._validate(Path(dispmap_az).exists(), result)
        if lv_theta is not None:
            result = self._validate(Path(lv_theta).exists(), result)
        if lv_phi is not None:
            result = self._validate(Path(lv_phi).exists(), result)
        if dv_norm is not None:
            Path(dv_norm).touch()
        if dv_theta is not None:
            Path(dv_theta).touch()
        if dv_phi is not None:
            Path(dv_phi).touch()
        if dv_x is not None:
            Path(dv_x).touch()
        if dv_y is not None:
            Path(dv_y).touch()
        if dv_z is not None:
            Path(dv_z).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def create_diff_par(self, PAR_1: str, PAR_2: str, DIFF_par: str, PAR_type = None, iflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "create_diff_par", supplied_args))

        if "create_diff_par" in self.call_count:
            self.call_count["create_diff_par"] += 1
        else:
            self.call_count["create_diff_par"] = 1

        if PAR_1 is not None:
            result = self._validate(Path(PAR_1).exists(), result)
        if PAR_2 is not None:
            result = self._validate(Path(PAR_2).exists(), result)
        if DIFF_par is not None and not Path(DIFF_par).exists():
            Path(DIFF_par).touch()
        valid_values = [0, 1, 2]
        result = self._validate(PAR_type in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(iflg in valid_values, result)
        return result

    def gc_map_grd(self, MLI_par: str, DEM_par: str, DEM: str, DEM_seg_par: str, DEM_seg: str, lookup_table: str, lat_ovr = None, lon_ovr = None, sim_sar: str = None, u: str = None, v: str = None, inc: str = None, psi: str = None, pix: str = None, ls_map: str = None, frame = None, ls_mode = None, r_ovr = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map_grd", supplied_args))

        if "gc_map_grd" in self.call_count:
            self.call_count["gc_map_grd"] += 1
        else:
            self.call_count["gc_map_grd"] = 1

        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if DEM_seg_par is not None and not Path(DEM_seg_par).exists():
            Path(DEM_seg_par).touch()
        if DEM_seg is not None:
            Path(DEM_seg).touch()
        if lookup_table is not None:
            Path(lookup_table).touch()
        if sim_sar is not None:
            Path(sim_sar).touch()
        if u is not None:
            Path(u).touch()
        if v is not None:
            Path(v).touch()
        if inc is not None:
            Path(inc).touch()
        if psi is not None:
            Path(psi).touch()
        if pix is not None:
            Path(pix).touch()
        if ls_map is not None:
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

        if data_tab is not None:
            result = self._validate(Path(data_tab).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        if DEM_par_out is not None:
            Path(DEM_par_out).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(format_flag in valid_values, result)
        return result

    def gc_map1(self, MLI_par: str, OFF_par: str, DEM_par: str, DEM: str, DEM_seg_par: str, DEM_seg: str, lookup_table: str, lat_ovr = None, lon_ovr = None, sim_sar: str = None, u: str = None, v: str = None, inc: str = None, psi: str = None, pix: str = None, ls_map: str = None, frame = None, ls_mode = None, r_ovr = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map1", supplied_args))

        if "gc_map1" in self.call_count:
            self.call_count["gc_map1"] += 1
        else:
            self.call_count["gc_map1"] = 1

        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if DEM_seg_par is not None and not Path(DEM_seg_par).exists():
            Path(DEM_seg_par).touch()
        if DEM_seg is not None:
            Path(DEM_seg).touch()
        if lookup_table is not None:
            Path(lookup_table).touch()
        if sim_sar is not None:
            Path(sim_sar).touch()
        if u is not None:
            Path(u).touch()
        if v is not None:
            Path(v).touch()
        if inc is not None:
            Path(inc).touch()
        if psi is not None:
            Path(psi).touch()
        if pix is not None:
            Path(pix).touch()
        if ls_map is not None:
            Path(ls_map).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(ls_mode in valid_values, result)
        return result

    def SLC_diff_intf(self, SLC_1: str, SLC_2R: str, SLC1_par: str, SLC2R_par: str, OFF_par: str, sim_unw: str, diff_int: str, rlks, azlks, sps_flg = None, azf_flg = None, rbw_min = None, rp1_flg = None, rp2_flg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "SLC_diff_intf", supplied_args))

        if "SLC_diff_intf" in self.call_count:
            self.call_count["SLC_diff_intf"] += 1
        else:
            self.call_count["SLC_diff_intf"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2R is not None:
            result = self._validate(Path(SLC_2R).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2R_par is not None:
            result = self._validate(Path(SLC2R_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if sim_unw is not None:
            result = self._validate(Path(sim_unw).exists(), result)
        if diff_int is not None:
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

    def pixel_area(self, MLI_par: str, DEM_par: str, DEM: str, lookup_table: str, ls_map: str, inc_map: str, pix_sigma0: str, pix_gamma0: str = None, nstep = None, area_fact = None, sigma0_ratio: str = None, gamma0_ratio: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "pixel_area", supplied_args))

        if "pixel_area" in self.call_count:
            self.call_count["pixel_area"] += 1
        else:
            self.call_count["pixel_area"] = 1

        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if lookup_table is not None:
            result = self._validate(Path(lookup_table).exists(), result)
        if ls_map is not None:
            result = self._validate(Path(ls_map).exists(), result)
        if inc_map is not None:
            result = self._validate(Path(inc_map).exists(), result)
        if pix_sigma0 is not None:
            Path(pix_sigma0).touch()
        if pix_gamma0 is not None:
            Path(pix_gamma0).touch()
        if sigma0_ratio is not None:
            Path(sigma0_ratio).touch()
        if gamma0_ratio is not None:
            Path(gamma0_ratio).touch()
        return result

    def extract_gcp(self, DEM_rdc: str, OFF_par: str, GCP: str, nr, naz, mask: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "extract_gcp", supplied_args))

        if "extract_gcp" in self.call_count:
            self.call_count["extract_gcp"] += 1
        else:
            self.call_count["extract_gcp"] = 1

        if DEM_rdc is not None:
            result = self._validate(Path(DEM_rdc).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if GCP is not None:
            Path(GCP).touch()
        if mask is not None:
            result = self._validate(Path(mask).exists(), result)
        return result

    def gc_map2(self, MLI_par: str, DEM_par: str, DEM: str, DEM_seg_par: str = None, DEM_seg: str = None, lookup_table: str = None, lat_ovr = None, lon_ovr = None, ls_map: str = None, ls_map_rdc: str = None, inc: str = None, res: str = None, offnadir: str = None, sim_sar: str = None, u: str = None, v: str = None, psi: str = None, pix: str = None, r_ovr = None, az_dec = None, mask = None, frame = None, ls_scaling = None, DIFF_par: str = None, ref_flg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map2", supplied_args))

        if "gc_map2" in self.call_count:
            self.call_count["gc_map2"] += 1
        else:
            self.call_count["gc_map2"] = 1

        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if DEM_seg_par is not None:
            Path(DEM_seg_par).touch()
        if DEM_seg is not None:
            Path(DEM_seg).touch()
        if lookup_table is not None:
            Path(lookup_table).touch()
        if ls_map is not None:
            Path(ls_map).touch()
        if ls_map_rdc is not None:
            Path(ls_map_rdc).touch()
        if inc is not None:
            Path(inc).touch()
        if res is not None:
            Path(res).touch()
        if offnadir is not None:
            Path(offnadir).touch()
        if sim_sar is not None:
            Path(sim_sar).touch()
        if u is not None:
            Path(u).touch()
        if v is not None:
            Path(v).touch()
        if psi is not None:
            Path(psi).touch()
        if pix is not None:
            Path(pix).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(mask in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(ls_scaling in valid_values, result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(ref_flg in valid_values, result)
        return result

    def dem_trans(self, DEM1_par: str, DEM1: str, DEM2_par: str, DEM2: str, lat_ovr = None, lon_ovr = None, datum_shift = None, bflg = None, lookup_table: str = None, interp_mode = None, order = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_trans", supplied_args))

        if "dem_trans" in self.call_count:
            self.call_count["dem_trans"] += 1
        else:
            self.call_count["dem_trans"] = 1

        if DEM1_par is not None:
            result = self._validate(Path(DEM1_par).exists(), result)
        if DEM1 is not None:
            result = self._validate(Path(DEM1).exists(), result)
        if DEM2_par is not None and not Path(DEM2_par).exists():
            Path(DEM2_par).touch()
        if DEM2 is not None:
            Path(DEM2).touch()
        valid_values = [0, 1]
        result = self._validate(datum_shift in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        if lookup_table is not None:
            Path(lookup_table).touch()
        valid_values = [0, 1, 2]
        result = self._validate(interp_mode in valid_values, result)
        return result

    def quad_sub(self, int_1: str, DIFF_par: str, int_2: str, int_type, mode = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "quad_sub", supplied_args))

        if "quad_sub" in self.call_count:
            self.call_count["quad_sub"] += 1
        else:
            self.call_count["quad_sub"] = 1

        if int_1 is not None:
            result = self._validate(Path(int_1).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if int_2 is not None:
            Path(int_2).touch()
        valid_values = [0, 1]
        result = self._validate(int_type in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def phase_sum(self, im_list: str, width, sum: str, start = None, nlines = None, pixav_x = None, pixav_y = None, zflag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "phase_sum", supplied_args))

        if "phase_sum" in self.call_count:
            self.call_count["phase_sum"] += 1
        else:
            self.call_count["phase_sum"] = 1

        if im_list is not None:
            result = self._validate(Path(im_list).exists(), result)
        if sum is not None:
            Path(sum).touch()
        valid_values = [0, 1]
        result = self._validate(zflag in valid_values, result)
        return result

    def base_add(self, base_1: str, base_2: str, base_out: str, mode = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "base_add", supplied_args))

        if "base_add" in self.call_count:
            self.call_count["base_add"] += 1
        else:
            self.call_count["base_add"] = 1

        if base_1 is not None:
            result = self._validate(Path(base_1).exists(), result)
        if base_2 is not None:
            result = self._validate(Path(base_2).exists(), result)
        if base_out is not None:
            Path(base_out).touch()
        return result

    def gec_map_grd(self, GRD_par: str, DEM_par: str, href: str, DEM_seg_par: str, lookup_table: str, lat_ovr = None, lon_ovr = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gec_map_grd", supplied_args))

        if "gec_map_grd" in self.call_count:
            self.call_count["gec_map_grd"] += 1
        else:
            self.call_count["gec_map_grd"] = 1

        if GRD_par is not None:
            result = self._validate(Path(GRD_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if href is not None:
            result = self._validate(Path(href).exists(), result)
        if DEM_seg_par is not None and not Path(DEM_seg_par).exists():
            Path(DEM_seg_par).touch()
        if lookup_table is not None:
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

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if SAR_coord is not None:
            result = self._validate(Path(SAR_coord).exists(), result)
        if MAP_coord is not None:
            Path(MAP_coord).touch()
        return result

    def dispmap_vec(self, DEM_par: str, dispmap: str, lv_theta: str, lv_phi: str, fv_theta: str, fv_phi: str, dv_norm: str, dv_theta: str = None, dv_phi: str = None, dv_x: str = None, dv_y: str = None, dv_z: str = None, mask_angle = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_vec", supplied_args))

        if "dispmap_vec" in self.call_count:
            self.call_count["dispmap_vec"] += 1
        else:
            self.call_count["dispmap_vec"] = 1

        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if dispmap is not None:
            result = self._validate(Path(dispmap).exists(), result)
        if lv_theta is not None:
            result = self._validate(Path(lv_theta).exists(), result)
        if lv_phi is not None:
            result = self._validate(Path(lv_phi).exists(), result)
        if fv_theta is not None:
            result = self._validate(Path(fv_theta).exists(), result)
        if fv_phi is not None:
            result = self._validate(Path(fv_phi).exists(), result)
        if dv_norm is not None:
            Path(dv_norm).touch()
        if dv_theta is not None:
            Path(dv_theta).touch()
        if dv_phi is not None:
            Path(dv_phi).touch()
        if dv_x is not None:
            Path(dv_x).touch()
        if dv_y is not None:
            Path(dv_y).touch()
        if dv_z is not None:
            Path(dv_z).touch()
        return result

    def scale_base(self, unw_2: str, scaled_unw_2: str, baseline_1: str, SLC1_par_1: str, OFF_par_1: str, baseline_2: str, SLC_1_par_2: str, OFF_par_2: str, int_type):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "scale_base", supplied_args))

        if "scale_base" in self.call_count:
            self.call_count["scale_base"] += 1
        else:
            self.call_count["scale_base"] = 1

        if unw_2 is not None:
            result = self._validate(Path(unw_2).exists(), result)
        if scaled_unw_2 is not None:
            Path(scaled_unw_2).touch()
        if baseline_1 is not None:
            result = self._validate(Path(baseline_1).exists(), result)
        if SLC1_par_1 is not None:
            result = self._validate(Path(SLC1_par_1).exists(), result)
        if OFF_par_1 is not None:
            result = self._validate(Path(OFF_par_1).exists(), result)
        if baseline_2 is not None:
            result = self._validate(Path(baseline_2).exists(), result)
        if SLC_1_par_2 is not None:
            result = self._validate(Path(SLC_1_par_2).exists(), result)
        if OFF_par_2 is not None:
            result = self._validate(Path(OFF_par_2).exists(), result)
        return result

    def par_EORC_PALSAR_geo(self, CEOS_leader: str, MLI_par: str, DEM_par: str, CEOS_data: str, MLI: str = None, cal = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_EORC_PALSAR_geo", supplied_args))

        if "par_EORC_PALSAR_geo" in self.call_count:
            self.call_count["par_EORC_PALSAR_geo"] += 1
        else:
            self.call_count["par_EORC_PALSAR_geo"] = 1

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if MLI_par is not None:
            Path(MLI_par).touch()
        if DEM_par is not None:
            Path(DEM_par).touch()
        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if MLI is not None:
            Path(MLI).touch()
        return result

    def SLC_intf_geo(self, SLC_1: str, SLC_2: str, DEM_par: str, interf: str, DEM_par2: str, e_lks, n_lks, MLI_1: str = None, MLI_2: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "SLC_intf_geo", supplied_args))

        if "SLC_intf_geo" in self.call_count:
            self.call_count["SLC_intf_geo"] += 1
        else:
            self.call_count["SLC_intf_geo"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if interf is not None:
            Path(interf).touch()
        if DEM_par2 is not None:
            Path(DEM_par2).touch()
        if MLI_1 is not None:
            Path(MLI_1).touch()
        if MLI_2 is not None:
            Path(MLI_2).touch()
        return result

    def map_trans(self, DEM1_par: str, data1: str, DEM2_par: str, data2: str, lat_ovr = None, lon_ovr = None, interp_mode = None, dtype = None, bflg = None, lookup_table: str = None, order = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "map_trans", supplied_args))

        if "map_trans" in self.call_count:
            self.call_count["map_trans"] += 1
        else:
            self.call_count["map_trans"] = 1

        if DEM1_par is not None:
            result = self._validate(Path(DEM1_par).exists(), result)
        if data1 is not None:
            result = self._validate(Path(data1).exists(), result)
        if DEM2_par is not None and not Path(DEM2_par).exists():
            Path(DEM2_par).touch()
        if data2 is not None:
            Path(data2).touch()
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        if lookup_table is not None:
            Path(lookup_table).touch()
        return result

    def gc_map_inversion(self, gc_map, width_in, gc_map_out, width_out, nlines_out = None, interp_mode = None, n_ovr = None, rad_max = None, nintr = None):
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

    def par_JERS_geo(self, CEOS_leader: str, CEOS_data: str, MLI_par: str, DEM_par: str, GEO: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_JERS_geo", supplied_args))

        if "par_JERS_geo" in self.call_count:
            self.call_count["par_JERS_geo"] += 1
        else:
            self.call_count["par_JERS_geo"] = 1

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if MLI_par is not None:
            Path(MLI_par).touch()
        if DEM_par is not None:
            Path(DEM_par).touch()
        if GEO is not None:
            Path(GEO).touch()
        return result

    def interp_data(self, data2: str, DIFF_par: str, data2_out: str, interp_mode = None, dtype = None, order = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "interp_data", supplied_args))

        if "interp_data" in self.call_count:
            self.call_count["interp_data"] += 1
        else:
            self.call_count["interp_data"] = 1

        if data2 is not None:
            result = self._validate(Path(data2).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if data2_out is not None:
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

        if MLI1_par is not None:
            result = self._validate(Path(MLI1_par).exists(), result)
        if DEM_RDC is not None:
            result = self._validate(Path(DEM_RDC).exists(), result)
        if MLI2_par is not None:
            result = self._validate(Path(MLI2_par).exists(), result)
        if lt is not None:
            Path(lt).touch()
        return result

    def SLC_interp_lt_ScanSAR(self, SLC2_tab: str, SLC2_par: str, SLC1_tab: str, SLC1_par: str, lookup_table: str, MLI1_par: str, MLI2_par: str, OFF_par: str, SLC2R_tab: str, SLC_2R: str = None, SLC2R_par: str = None, mode = None, order = None, SLC2R_dir = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "SLC_interp_lt_ScanSAR", supplied_args))

        if "SLC_interp_lt_ScanSAR" in self.call_count:
            self.call_count["SLC_interp_lt_ScanSAR"] += 1
        else:
            self.call_count["SLC_interp_lt_ScanSAR"] = 1

        if SLC2_tab is not None:
            result = self._validate(Path(SLC2_tab).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if SLC1_tab is not None:
            result = self._validate(Path(SLC1_tab).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if lookup_table is not None:
            result = self._validate(Path(lookup_table).exists(), result)
        if MLI1_par is not None:
            result = self._validate(Path(MLI1_par).exists(), result)
        if MLI2_par is not None:
            result = self._validate(Path(MLI2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if SLC2R_tab is not None and not Path(SLC2R_tab).exists():
            Path(SLC2R_tab).touch()
        if SLC_2R is not None:
            Path(SLC_2R).touch()
        if SLC2R_par is not None:
            Path(SLC2R_par).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def par_CS_geo(self, HDF5: str, MLI_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_CS_geo", supplied_args))

        if "par_CS_geo" in self.call_count:
            self.call_count["par_CS_geo"] += 1
        else:
            self.call_count["par_CS_geo"] = 1

        if HDF5 is not None:
            result = self._validate(Path(HDF5).exists(), result)
        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        return result

    def par_RCM_geo(self, RCM_dir: str, polarization, MLI_par: str, DEM_par: str, GEO: str, dtype = None, ps = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_RCM_geo", supplied_args))

        if "par_RCM_geo" in self.call_count:
            self.call_count["par_RCM_geo"] += 1
        else:
            self.call_count["par_RCM_geo"] = 1

        if RCM_dir is not None:
            result = self._validate(Path(RCM_dir).exists(), result)
        if MLI_par is not None:
            Path(MLI_par).touch()
        if DEM_par is not None:
            Path(DEM_par).touch()
        if GEO is not None:
            Path(GEO).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def offset_pwr_trackingm2(self, MLI_1: str, MLI_2: str, DIFF_par: str, offs: str, ccp: str, DIFF_par2: str = None, offs2: str = None, rwin = None, azwin = None, offsets: str = None, n_ovr = None, thres = None, rstep = None, azstep = None, rstart = None, rstop = None, azstart = None, azstop = None, bw_frac = None, pflag = None, pltflg = None, ccs: str = None, std_mean = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_pwr_trackingm2", supplied_args))

        if "offset_pwr_trackingm2" in self.call_count:
            self.call_count["offset_pwr_trackingm2"] += 1
        else:
            self.call_count["offset_pwr_trackingm2"] = 1

        if MLI_1 is not None:
            result = self._validate(Path(MLI_1).exists(), result)
        if MLI_2 is not None:
            result = self._validate(Path(MLI_2).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if offs is not None:
            Path(offs).touch()
        if ccp is not None:
            Path(ccp).touch()
        if DIFF_par2 is not None:
            result = self._validate(Path(DIFF_par2).exists(), result)
        if offs2 is not None:
            result = self._validate(Path(offs2).exists(), result)
        if offsets is not None:
            Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        if ccs is not None:
            Path(ccs).touch()
        return result

    def init_offset_orbitm(self, MLI1_par: str, MLI2_par: str, DIFF_par: str, rpos = None, azpos = None, cflag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "init_offset_orbitm", supplied_args))

        if "init_offset_orbitm" in self.call_count:
            self.call_count["init_offset_orbitm"] += 1
        else:
            self.call_count["init_offset_orbitm"] = 1

        if MLI1_par is not None:
            result = self._validate(Path(MLI1_par).exists(), result)
        if MLI2_par is not None:
            result = self._validate(Path(MLI2_par).exists(), result)
        if DIFF_par is not None and not Path(DIFF_par).exists():
            Path(DIFF_par).touch()
        valid_values = [0, 1]
        result = self._validate(cflag in valid_values, result)
        return result

    def quad_fit(self, unw: str, DIFF_par: str, dr = None, daz = None, mask = None, plot_data: str = None, model = None, pmodel: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "quad_fit", supplied_args))

        if "quad_fit" in self.call_count:
            self.call_count["quad_fit"] += 1
        else:
            self.call_count["quad_fit"] = 1

        if unw is not None:
            result = self._validate(Path(unw).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if plot_data is not None:
            Path(plot_data).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6, 7]
        result = self._validate(model in valid_values, result)
        if pmodel is not None:
            Path(pmodel).touch()
        return result

    def stacking(self, DIFF_tab: str, width, ph_rate: str, sig_ph_rate: str, sig_ph: str, roff, loff, nr = None, nl = None, np_min = None, tscale = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "stacking", supplied_args))

        if "stacking" in self.call_count:
            self.call_count["stacking"] += 1
        else:
            self.call_count["stacking"] = 1

        if DIFF_tab is not None:
            result = self._validate(Path(DIFF_tab).exists(), result)
        if ph_rate is not None:
            Path(ph_rate).touch()
        if sig_ph_rate is not None:
            Path(sig_ph_rate).touch()
        if sig_ph is not None:
            Path(sig_ph).touch()
        valid_values = [0, 1]
        result = self._validate(tscale in valid_values, result)
        return result

    def atm_mod2(self, diff_unw: str, hgt: str, MLI_par: str, model: str, dr = None, daz = None, mask: str = None, mode = None, roff = None, loff = None, report: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "atm_mod2", supplied_args))

        if "atm_mod2" in self.call_count:
            self.call_count["atm_mod2"] += 1
        else:
            self.call_count["atm_mod2"] = 1

        if diff_unw is not None:
            result = self._validate(Path(diff_unw).exists(), result)
        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if model is not None:
            Path(model).touch()
        if mask is not None:
            result = self._validate(Path(mask).exists(), result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(mode in valid_values, result)
        if report is not None:
            Path(report).touch()
        return result

    def par_UAVSAR_geo(self, ann: str, SLC_MLI_par: str = None, DEM_par: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_UAVSAR_geo", supplied_args))

        if "par_UAVSAR_geo" in self.call_count:
            self.call_count["par_UAVSAR_geo"] += 1
        else:
            self.call_count["par_UAVSAR_geo"] = 1

        if ann is not None:
            result = self._validate(Path(ann).exists(), result)
        if SLC_MLI_par is not None:
            Path(SLC_MLI_par).touch()
        if DEM_par is not None:
            Path(DEM_par).touch()
        return result

    def dem_import(self, input_DEM: str, DEM: str, DEM_par: str, input_type = None, priority = None, geoid: str = None, geoid_par: str = None, geoid_type = None, latN_shift = None, lonE_shift = None, zflg = None, no_data = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_import", supplied_args))

        if "dem_import" in self.call_count:
            self.call_count["dem_import"] += 1
        else:
            self.call_count["dem_import"] = 1

        if input_DEM is not None:
            result = self._validate(Path(input_DEM).exists(), result)
        if DEM is not None:
            Path(DEM).touch()
        if DEM_par is not None and not Path(DEM_par).exists():
            Path(DEM_par).touch()
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(input_type in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(priority in valid_values, result)
        if geoid is not None:
            result = self._validate(Path(geoid).exists(), result)
        if geoid_par is not None:
            result = self._validate(Path(geoid_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(geoid_type in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def SLC_interp_lt(self, SLC_2: str, SLC1_par: str, SLC2_par: str, lookup_table: str, MLI1_par: str, MLI2_par: str, OFF_par: str, SLC_2R: str, SLC2R_par: str, blk_size = None, mode = None, order = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "SLC_interp_lt", supplied_args))

        if "SLC_interp_lt" in self.call_count:
            self.call_count["SLC_interp_lt"] += 1
        else:
            self.call_count["SLC_interp_lt"] = 1

        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if lookup_table is not None:
            result = self._validate(Path(lookup_table).exists(), result)
        if MLI1_par is not None:
            result = self._validate(Path(MLI1_par).exists(), result)
        if MLI2_par is not None:
            result = self._validate(Path(MLI2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if SLC_2R is not None:
            Path(SLC_2R).touch()
        if SLC2R_par is not None:
            Path(SLC2R_par).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def dispmap(self, unw: str, hgt: str, MLI_par: str, OFF_par: str, disp_map: str, mode = None, sflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap", supplied_args))

        if "dispmap" in self.call_count:
            self.call_count["dispmap"] += 1
        else:
            self.call_count["dispmap"] = 1

        if unw is not None:
            result = self._validate(Path(unw).exists(), result)
        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if disp_map is not None:
            Path(disp_map).touch()
        return result

    def atm_mod(self, diff_unw: str, hgt: str, DIFF_par: str, model: str, dr = None, daz = None, mask: str = None, mode = None, roff = None, loff = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "atm_mod", supplied_args))

        if "atm_mod" in self.call_count:
            self.call_count["atm_mod"] += 1
        else:
            self.call_count["atm_mod"] = 1

        if diff_unw is not None:
            result = self._validate(Path(diff_unw).exists(), result)
        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if model is not None:
            Path(model).touch()
        if mask is not None:
            result = self._validate(Path(mask).exists(), result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def dispmap_LOS(self, unw: str, width, freq, disp_map: str, sflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_LOS", supplied_args))

        if "dispmap_LOS" in self.call_count:
            self.call_count["dispmap_LOS"] += 1
        else:
            self.call_count["dispmap_LOS"] = 1

        if unw is not None:
            result = self._validate(Path(unw).exists(), result)
        if disp_map is not None:
            Path(disp_map).touch()
        return result

    def sub_phase(self, int_1: str, unw_2: str, DIFF_par: str, diff_int: str, dtype, mode = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "sub_phase", supplied_args))

        if "sub_phase" in self.call_count:
            self.call_count["sub_phase"] += 1
        else:
            self.call_count["sub_phase"] = 1

        if int_1 is not None:
            result = self._validate(Path(int_1).exists(), result)
        if unw_2 is not None:
            result = self._validate(Path(unw_2).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if diff_int is not None:
            Path(diff_int).touch()
        valid_values = [0, 1, 2]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def phase_sim(self, SLC1_par: str, OFF_par: str, baseline: str, hgt: str, sim_unw: str, ph_flag = None, bflag = None, definition: str = None, delta_t: str = None, int_mode: str = None, SLC2R_par: str = None, ph_mode = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "phase_sim", supplied_args))

        if "phase_sim" in self.call_count:
            self.call_count["phase_sim"] += 1
        else:
            self.call_count["phase_sim"] = 1

        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if baseline is not None:
            result = self._validate(Path(baseline).exists(), result)
        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if sim_unw is not None:
            Path(sim_unw).touch()
        valid_values = [0, 1]
        result = self._validate(ph_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(bflag in valid_values, result)
        if definition is not None:
            result = self._validate(Path(definition).exists(), result)
        if delta_t is not None:
            result = self._validate(Path(delta_t).exists(), result)
        if int_mode is not None:
            result = self._validate(Path(int_mode).exists(), result)
        if SLC2R_par is not None:
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

        if annotation_XML is not None:
            result = self._validate(Path(annotation_XML).exists(), result)
        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if polarization is not None:
            result = self._validate(Path(polarization).exists(), result)
        if DEM_par is not None:
            Path(DEM_par).touch()
        if MLI_par is not None:
            Path(MLI_par).touch()
        if MLI is not None:
            Path(MLI).touch()
        return result

    def MLI_interp_lt(self, MLI_2: str, MLI1_par: str, MLI2_par: str, lookup_table: str, MLI3_par: str, MLI4_par: str, DIFF_par: str, MLI_2R: str, MLI2R_par: str, order = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "MLI_interp_lt", supplied_args))

        if "MLI_interp_lt" in self.call_count:
            self.call_count["MLI_interp_lt"] += 1
        else:
            self.call_count["MLI_interp_lt"] = 1

        if MLI_2 is not None:
            result = self._validate(Path(MLI_2).exists(), result)
        if MLI1_par is not None:
            result = self._validate(Path(MLI1_par).exists(), result)
        if MLI2_par is not None:
            result = self._validate(Path(MLI2_par).exists(), result)
        if lookup_table is not None:
            result = self._validate(Path(lookup_table).exists(), result)
        if MLI3_par is not None:
            result = self._validate(Path(MLI3_par).exists(), result)
        if MLI4_par is not None:
            result = self._validate(Path(MLI4_par).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if MLI_2R is not None:
            Path(MLI_2R).touch()
        if MLI2R_par is not None:
            Path(MLI2R_par).touch()
        return result

    def lk_vec_lt(self, MLI_par: str, DEM_par: str, DEM: str, lt: str, lv_theta: str, lv_phi: str, lv_ENU: str = None, azv_ENU: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "lk_vec_lt", supplied_args))

        if "lk_vec_lt" in self.call_count:
            self.call_count["lk_vec_lt"] += 1
        else:
            self.call_count["lk_vec_lt"] = 1

        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if lt is not None:
            result = self._validate(Path(lt).exists(), result)
        if lv_theta is not None:
            Path(lv_theta).touch()
        if lv_phi is not None:
            Path(lv_phi).touch()
        if lv_ENU is not None:
            Path(lv_ENU).touch()
        if azv_ENU is not None:
            Path(azv_ENU).touch()
        return result

    def coord_to_sarpix(self, SLC_par, OFF_par: str, DEM_par: str, north_lat: str = None, east_lon: str = None, hgt: str = None, DIFF_par = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "coord_to_sarpix", supplied_args))

        if "coord_to_sarpix" in self.call_count:
            self.call_count["coord_to_sarpix"] += 1
        else:
            self.call_count["coord_to_sarpix"] = 1

        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if north_lat is not None:
            result = self._validate(Path(north_lat).exists(), result)
        if east_lon is not None:
            result = self._validate(Path(east_lon).exists(), result)
        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        return result

    def WSS_intf(self, SLC_1: str, SLC_2R: str, SLC1_par: str, SLC2R_par: str, OFF_par: str, interf: str, rlks = None, sps_flg = None, azf_flg = None, m_flg = None, boff = None, bstep = None, bmax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "WSS_intf", supplied_args))

        if "WSS_intf" in self.call_count:
            self.call_count["WSS_intf"] += 1
        else:
            self.call_count["WSS_intf"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2R is not None:
            result = self._validate(Path(SLC_2R).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2R_par is not None:
            result = self._validate(Path(SLC2R_par).exists(), result)
        if OFF_par is not None:
            Path(OFF_par).touch()
        if interf is not None:
            Path(interf).touch()
        valid_values = [1, 0]
        result = self._validate(sps_flg in valid_values, result)
        valid_values = [1, 0]
        result = self._validate(azf_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(m_flg in valid_values, result)
        return result

    def map_section(self, DEM_par: str, n1, e1, n2, e2, post_north, post_east, DEM_par2: str, lt: str = None, MLI_par1: str = None, MLI_par2: str = None, cflg = None, lt2: str = None, MLI_coord: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "map_section", supplied_args))

        if "map_section" in self.call_count:
            self.call_count["map_section"] += 1
        else:
            self.call_count["map_section"] = 1

        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM_par2 is not None:
            Path(DEM_par2).touch()
        if lt is not None:
            result = self._validate(Path(lt).exists(), result)
        if MLI_par1 is not None:
            result = self._validate(Path(MLI_par1).exists(), result)
        if MLI_par2 is not None:
            result = self._validate(Path(MLI_par2).exists(), result)
        if lt2 is not None:
            Path(lt2).touch()
        if MLI_coord is not None:
            Path(MLI_coord).touch()
        return result

    def offset_list_fitm(self, cp_list: str, DIFF_par: str, DEM_par: str, lookup_table: str = None, lt_type = None, type1 = None, type2 = None, coffsets: str = None, poly_order = None, interact_flag = None, trans_list = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_list_fitm", supplied_args))

        if "offset_list_fitm" in self.call_count:
            self.call_count["offset_list_fitm"] += 1
        else:
            self.call_count["offset_list_fitm"] = 1

        if cp_list is not None:
            result = self._validate(Path(cp_list).exists(), result)
        if DIFF_par is not None and not Path(DIFF_par).exists():
            Path(DIFF_par).touch()
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if lookup_table is not None:
            result = self._validate(Path(lookup_table).exists(), result)
        valid_values = [1, 2]
        result = self._validate(lt_type in valid_values, result)
        valid_values = [1, 2, 3]
        result = self._validate(type1 in valid_values, result)
        valid_values = [1, 2, 3]
        result = self._validate(type2 in valid_values, result)
        if coffsets is not None:
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

        if DEM_par1 is not None:
            result = self._validate(Path(DEM_par1).exists(), result)
        if gc_map is not None:
            result = self._validate(Path(gc_map).exists(), result)
        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if mask is not None:
            result = self._validate(Path(mask).exists(), result)
        if clist_RDC is not None:
            Path(clist_RDC).touch()
        if clist_MAP is not None:
            Path(clist_MAP).touch()
        if DEM_par2 is not None:
            Path(DEM_par2).touch()
        return result

    def multi_look_geo(self, geo_SLC: str, SLC_DEM_par: str, MLI: str, MLI_DEM_par: str, e_lks, n_lks, dtype = None, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "multi_look_geo", supplied_args))

        if "multi_look_geo" in self.call_count:
            self.call_count["multi_look_geo"] += 1
        else:
            self.call_count["multi_look_geo"] = 1

        if geo_SLC is not None:
            result = self._validate(Path(geo_SLC).exists(), result)
        if SLC_DEM_par is not None:
            result = self._validate(Path(SLC_DEM_par).exists(), result)
        if MLI is not None:
            Path(MLI).touch()
        if MLI_DEM_par is not None:
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

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if lookup_table is not None:
            result = self._validate(Path(lookup_table).exists(), result)
        if MLI1_par is not None:
            result = self._validate(Path(MLI1_par).exists(), result)
        if MLI2_par is not None:
            result = self._validate(Path(MLI2_par).exists(), result)
        if DIFF_par1 is not None:
            result = self._validate(Path(DIFF_par1).exists(), result)
        if SLC_2R is not None:
            Path(SLC_2R).touch()
        if SLC2R_par is not None:
            Path(SLC2R_par).touch()
        if DIFF_par2 is not None:
            Path(DIFF_par2).touch()
        return result

    def dispmap_vec2(self, DEM_par: str, DEM: str, dispmap1: str, lv1_theta: str, lv1_phi: str, dispmap2: str, lv2_theta: str, lv2_phi: str, dv_norm: str, dv_theta: str = None, dv_phi: str = None, dv_x: str = None, dv_y: str = None, dv_z: str = None, mask_angle = None, mode = None, ax_north = None, ax_east = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_vec2", supplied_args))

        if "dispmap_vec2" in self.call_count:
            self.call_count["dispmap_vec2"] += 1
        else:
            self.call_count["dispmap_vec2"] = 1

        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if dispmap1 is not None:
            result = self._validate(Path(dispmap1).exists(), result)
        if lv1_theta is not None:
            result = self._validate(Path(lv1_theta).exists(), result)
        if lv1_phi is not None:
            result = self._validate(Path(lv1_phi).exists(), result)
        if dispmap2 is not None:
            result = self._validate(Path(dispmap2).exists(), result)
        if lv2_theta is not None:
            result = self._validate(Path(lv2_theta).exists(), result)
        if lv2_phi is not None:
            result = self._validate(Path(lv2_phi).exists(), result)
        if dv_norm is not None:
            Path(dv_norm).touch()
        if dv_theta is not None:
            Path(dv_theta).touch()
        if dv_phi is not None:
            Path(dv_phi).touch()
        if dv_x is not None:
            Path(dv_x).touch()
        if dv_y is not None:
            Path(dv_y).touch()
        if dv_z is not None:
            Path(dv_z).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        return result

    def diff_ls_unw(self, int_1: str, unw_2: str, DIFF_par: str, diff_int: str, int_type = None, ph_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "diff_ls_unw", supplied_args))

        if "diff_ls_unw" in self.call_count:
            self.call_count["diff_ls_unw"] += 1
        else:
            self.call_count["diff_ls_unw"] = 1

        if int_1 is not None:
            result = self._validate(Path(int_1).exists(), result)
        if unw_2 is not None:
            result = self._validate(Path(unw_2).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if diff_int is not None:
            Path(diff_int).touch()
        return result

    def offset_pwr_list(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, clist_RDC: str, clist_MAP: str, offs: str, ccp: str, nx, ny, rwin = None, azwin = None, offsets: str = None, n_ovr = None, thres = None, bw_frac = None, deramp = None, int_filt = None, pflag = None, pltflg = None, ccs: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_pwr_list", supplied_args))

        if "offset_pwr_list" in self.call_count:
            self.call_count["offset_pwr_list"] += 1
        else:
            self.call_count["offset_pwr_list"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if clist_RDC is not None:
            result = self._validate(Path(clist_RDC).exists(), result)
        if clist_MAP is not None:
            result = self._validate(Path(clist_MAP).exists(), result)
        if offs is not None:
            Path(offs).touch()
        if ccp is not None:
            Path(ccp).touch()
        if offsets is not None:
            Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(int_filt in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        if ccs is not None:
            Path(ccs).touch()
        return result

    def geocode(self, lookup_table: str, data_in: str, width_in, data_out: str, width_out, nlines_out = None, interp_mode = None, dtype = None, lr_in = None, lr_out = None, n_ovr = None, rad_max = None, nintr = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "geocode", supplied_args))

        if "geocode" in self.call_count:
            self.call_count["geocode"] += 1
        else:
            self.call_count["geocode"] = 1

        if lookup_table is not None:
            result = self._validate(Path(lookup_table).exists(), result)
        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(interp_mode in valid_values, result)
        valid_values = [0, 1, 2, 3, 4, 5, 6]
        result = self._validate(dtype in valid_values, result)
        return result

    def par_TX_geo(self, annotation_XML: str, GeoTIFF: str, MLI_par: str, DEM_par: str, GEO: str, pol = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_TX_geo", supplied_args))

        if "par_TX_geo" in self.call_count:
            self.call_count["par_TX_geo"] += 1
        else:
            self.call_count["par_TX_geo"] = 1

        if annotation_XML is not None:
            result = self._validate(Path(annotation_XML).exists(), result)
        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if MLI_par is not None:
            Path(MLI_par).touch()
        if DEM_par is not None:
            Path(DEM_par).touch()
        if GEO is not None:
            Path(GEO).touch()
        return result

    def data2xyz(self, DEM_par: str, data: str, data_xyz: str, dflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "data2xyz", supplied_args))

        if "data2xyz" in self.call_count:
            self.call_count["data2xyz"] += 1
        else:
            self.call_count["data2xyz"] = 1

        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if data_xyz is not None:
            Path(data_xyz).touch()
        valid_values = [0, 1]
        result = self._validate(dflg in valid_values, result)
        return result

    def ScanSAR_burst_diff_intf(self, SLC1_tab: str, SLC2R_tab: str, SIM_tab: str, DIFF_tab: str, SLCR_tab: str = None, DIFF_dir = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "ScanSAR_burst_diff_intf", supplied_args))

        if "ScanSAR_burst_diff_intf" in self.call_count:
            self.call_count["ScanSAR_burst_diff_intf"] += 1
        else:
            self.call_count["ScanSAR_burst_diff_intf"] = 1

        if SLC1_tab is not None:
            result = self._validate(Path(SLC1_tab).exists(), result)
        if SLC2R_tab is not None:
            result = self._validate(Path(SLC2R_tab).exists(), result)
        if SIM_tab is not None:
            result = self._validate(Path(SIM_tab).exists(), result)
        if DIFF_tab is not None and not Path(DIFF_tab).exists():
            Path(DIFF_tab).touch()
        if SLCR_tab is not None:
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

        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if DEM_XYZ is not None:
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

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if lv_theta is not None:
            Path(lv_theta).touch()
        if lv_phi is not None:
            Path(lv_phi).touch()
        return result

    def create_dem_par(self, DEM_par: str, SLC_par: str = None, terra_alt = None, delta_y = None, delta_x = None, EPSG = None, iflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "create_dem_par", supplied_args))

        if "create_dem_par" in self.call_count:
            self.call_count["create_dem_par"] += 1
        else:
            self.call_count["create_dem_par"] = 1

        if DEM_par is not None and not Path(DEM_par).exists():
            Path(DEM_par).touch()
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(iflg in valid_values, result)
        return result

    def dem_x_y_z(self, DEM_par: str, DEM: str, DEM_X: str, DEM_Y: str, DEM_Z: str, format_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_x_y_z", supplied_args))

        if "dem_x_y_z" in self.call_count:
            self.call_count["dem_x_y_z"] += 1
        else:
            self.call_count["dem_x_y_z"] = 1

        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if DEM_X is not None:
            Path(DEM_X).touch()
        if DEM_Y is not None:
            Path(DEM_Y).touch()
        if DEM_Z is not None:
            Path(DEM_Z).touch()
        valid_values = [0, 1]
        result = self._validate(format_flag in valid_values, result)
        return result

    def ras_clist(self, clist: str, ras_in: str, ras_out: str, xsf = None, ysf = None, r = None, g = None, b = None, xs = None, zflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "ras_clist", supplied_args))

        if "ras_clist" in self.call_count:
            self.call_count["ras_clist"] += 1
        else:
            self.call_count["ras_clist"] = 1

        if clist is not None:
            result = self._validate(Path(clist).exists(), result)
        if ras_in is not None:
            result = self._validate(Path(ras_in).exists(), result)
        if ras_out is not None:
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

        if LV is not None:
            result = self._validate(Path(LV).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if disp_east is not None:
            result = self._validate(Path(disp_east).exists(), result)
        if disp_north is not None:
            result = self._validate(Path(disp_north).exists(), result)
        if disp_up is not None:
            result = self._validate(Path(disp_up).exists(), result)
        if disp_LOS is not None:
            Path(disp_LOS).touch()
        return result

    def dem_gradient(self, DEM_par: str, DEM: str, theta: str, phi: str, mag: str, type = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_gradient", supplied_args))

        if "dem_gradient" in self.call_count:
            self.call_count["dem_gradient"] += 1
        else:
            self.call_count["dem_gradient"] = 1

        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if theta is not None:
            Path(theta).touch()
        if phi is not None:
            Path(phi).touch()
        if mag is not None:
            Path(mag).touch()
        valid_values = [0, 1]
        result = self._validate(type in valid_values, result)
        return result

    def offset_pwr_trackingm(self, MLI_1: str, MLI_2: str, DIFF_par: str, offs: str, ccp: str, rwin = None, azwin = None, offsets: str = None, n_ovr = None, thres = None, rstep = None, azstep = None, rstart = None, rstop = None, azstart = None, azstop = None, lanczos = None, bw_frac = None, pflag = None, pltflg = None, ccs: str = None, std_mean = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_pwr_trackingm", supplied_args))

        if "offset_pwr_trackingm" in self.call_count:
            self.call_count["offset_pwr_trackingm"] += 1
        else:
            self.call_count["offset_pwr_trackingm"] = 1

        if MLI_1 is not None:
            result = self._validate(Path(MLI_1).exists(), result)
        if MLI_2 is not None:
            result = self._validate(Path(MLI_2).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if offs is not None:
            Path(offs).touch()
        if ccp is not None:
            Path(ccp).touch()
        if offsets is not None:
            Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        if ccs is not None:
            Path(ccs).touch()
        return result

    def init_offsetm(self, MLI_1: str, MLI_2: str, DIFF_par, rlks = None, azlks = None, rpos = None, azpos = None, offr = None, offaz = None, thres = None, patch = None, cflag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "init_offsetm", supplied_args))

        if "init_offsetm" in self.call_count:
            self.call_count["init_offsetm"] += 1
        else:
            self.call_count["init_offsetm"] = 1

        if MLI_1 is not None:
            result = self._validate(Path(MLI_1).exists(), result)
        if MLI_2 is not None:
            result = self._validate(Path(MLI_2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(cflag in valid_values, result)
        return result

    def gc_map_fd(self, MLI_par: str, fd_tab: str, DEM_par: str, DEM: str, DEM_seg_par: str, DEM_seg: str, lookup_table: str, lat_ovr = None, lon_ovr = None, sim_sar: str = None, u: str = None, v: str = None, inc: str = None, psi: str = None, pix: str = None, ls_map: str = None, frame = None, ls_mode = None, r_ovr = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_map_fd", supplied_args))

        if "gc_map_fd" in self.call_count:
            self.call_count["gc_map_fd"] += 1
        else:
            self.call_count["gc_map_fd"] = 1

        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if fd_tab is not None:
            result = self._validate(Path(fd_tab).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if DEM_seg_par is not None and not Path(DEM_seg_par).exists():
            Path(DEM_seg_par).touch()
        if DEM_seg is not None:
            Path(DEM_seg).touch()
        if lookup_table is not None:
            Path(lookup_table).touch()
        if sim_sar is not None:
            Path(sim_sar).touch()
        if u is not None:
            Path(u).touch()
        if v is not None:
            Path(v).touch()
        if inc is not None:
            Path(inc).touch()
        if psi is not None:
            Path(psi).touch()
        if pix is not None:
            Path(pix).touch()
        if ls_map is not None:
            Path(ls_map).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(ls_mode in valid_values, result)
        return result

    def phase_sim_orb(self, SLC1_par: str, SLC2R_par: str, OFF_par: str, hgt: str, sim_unw: str, SLC_ref_par: str = None, definition: str = None, delta_t: str = None, int_mode: str = None, ph_mode = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "phase_sim_orb", supplied_args))

        if "phase_sim_orb" in self.call_count:
            self.call_count["phase_sim_orb"] += 1
        else:
            self.call_count["phase_sim_orb"] = 1

        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2R_par is not None:
            result = self._validate(Path(SLC2R_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if sim_unw is not None:
            Path(sim_unw).touch()
        if SLC_ref_par is not None:
            result = self._validate(Path(SLC_ref_par).exists(), result)
        if definition is not None:
            result = self._validate(Path(definition).exists(), result)
        if delta_t is not None:
            result = self._validate(Path(delta_t).exists(), result)
        if int_mode is not None:
            result = self._validate(Path(int_mode).exists(), result)
        valid_values = [0, 1]
        result = self._validate(ph_mode in valid_values, result)
        return result

    def gec_map(self, SLC_par: str, OFF_par: str, DEM_par: str, href: str, DEM_seg_par: str, lookup_table: str, lat_ovr = None, lon_ovr = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gec_map", supplied_args))

        if "gec_map" in self.call_count:
            self.call_count["gec_map"] += 1
        else:
            self.call_count["gec_map"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if href is not None:
            result = self._validate(Path(href).exists(), result)
        if DEM_seg_par is not None and not Path(DEM_seg_par).exists():
            Path(DEM_seg_par).touch()
        if lookup_table is not None:
            Path(lookup_table).touch()
        return result

    def dh_map_orb(self, SLC1_par: str, SLC2R_par: str, OFF_par: str, hgt: str, dp: str, dpdh: str, dh: str, SLC_ref_par: str = None, int_mode: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dh_map_orb", supplied_args))

        if "dh_map_orb" in self.call_count:
            self.call_count["dh_map_orb"] += 1
        else:
            self.call_count["dh_map_orb"] = 1

        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2R_par is not None:
            result = self._validate(Path(SLC2R_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if dp is not None:
            result = self._validate(Path(dp).exists(), result)
        if dpdh is not None:
            Path(dpdh).touch()
        if dh is not None:
            Path(dh).touch()
        if SLC_ref_par is not None:
            result = self._validate(Path(SLC_ref_par).exists(), result)
        if int_mode is not None:
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

        if offs is not None:
            result = self._validate(Path(offs).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if offs_sub is not None:
            Path(offs_sub).touch()
        return result

    def offset_trackingm(self, offs: str, snr: str, MLI_par: str, DIFF_par: str, coffs_map: str, coffsets: str = None, mode = None, thres = None, poly_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_trackingm", supplied_args))

        if "offset_trackingm" in self.call_count:
            self.call_count["offset_trackingm"] += 1
        else:
            self.call_count["offset_trackingm"] = 1

        if offs is not None:
            result = self._validate(Path(offs).exists(), result)
        if snr is not None:
            result = self._validate(Path(snr).exists(), result)
        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if coffs_map is not None:
            Path(coffs_map).touch()
        if coffsets is not None:
            Path(coffsets).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(poly_flag in valid_values, result)
        return result

    def comb_interfs(self, int_1, int_2, base_1, base_2, factor_1, factor_2, width, combi_out, combi_base, sm = None, Only = None, The = None):
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

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if lookup_table is not None:
            Path(lookup_table).touch()
        return result

    def par_KS_geo(self, HDF5: str, MLI_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "par_KS_geo", supplied_args))

        if "par_KS_geo" in self.call_count:
            self.call_count["par_KS_geo"] += 1
        else:
            self.call_count["par_KS_geo"] = 1

        if HDF5 is not None:
            result = self._validate(Path(HDF5).exists(), result)
        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        return result

    def offset_pwrm(self, MLI_1: str, MLI_2: str, DIFF_par: str, offs: str, ccp: str, rwin = None, azwin = None, offsets: str = None, n_ovr = None, nr = None, naz = None, thres = None, lanczos = None, bw_frac = None, pflag = None, pltflg = None, ccs: str = None, std_mean = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_pwrm", supplied_args))

        if "offset_pwrm" in self.call_count:
            self.call_count["offset_pwrm"] += 1
        else:
            self.call_count["offset_pwrm"] = 1

        if MLI_1 is not None:
            result = self._validate(Path(MLI_1).exists(), result)
        if MLI_2 is not None:
            result = self._validate(Path(MLI_2).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if offs is not None:
            Path(offs).touch()
        if ccp is not None:
            Path(ccp).touch()
        if offsets is not None:
            Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        if ccs is not None:
            Path(ccs).touch()
        return result

    def coord_to_sarpix_list(self, SLC_par: str, OFF_par: str, DEM_par: str, MAP_coord: str, SAR_coord: str, DIFF_par: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "coord_to_sarpix_list", supplied_args))

        if "coord_to_sarpix_list" in self.call_count:
            self.call_count["coord_to_sarpix_list"] += 1
        else:
            self.call_count["coord_to_sarpix_list"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if MAP_coord is not None:
            result = self._validate(Path(MAP_coord).exists(), result)
        if SAR_coord is not None:
            Path(SAR_coord).touch()
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        return result

    def offset_fitm(self, offs: str, ccp: str, DIFF_par: str, coffs: str = None, coffsets: str = None, thres = None, npoly = None, interact_mode = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "offset_fitm", supplied_args))

        if "offset_fitm" in self.call_count:
            self.call_count["offset_fitm"] += 1
        else:
            self.call_count["offset_fitm"] = 1

        if offs is not None:
            result = self._validate(Path(offs).exists(), result)
        if ccp is not None:
            result = self._validate(Path(ccp).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if coffs is not None:
            Path(coffs).touch()
        if coffsets is not None:
            Path(coffsets).touch()
        valid_values = [0, 1]
        result = self._validate(interact_mode in valid_values, result)
        return result

    def dispmap_ENU(self, LV_tab: str, DISP_tab: str, SIGMA_tab: str, DEM_par: str, disp_east: str, disp_north: str, disp_up: str, sigma_east: str = None, sigma_north: str = None, sigma_up: str = None, chi2: str = None, min_obs = None, tol = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dispmap_ENU", supplied_args))

        if "dispmap_ENU" in self.call_count:
            self.call_count["dispmap_ENU"] += 1
        else:
            self.call_count["dispmap_ENU"] = 1

        if LV_tab is not None:
            result = self._validate(Path(LV_tab).exists(), result)
        if DISP_tab is not None:
            result = self._validate(Path(DISP_tab).exists(), result)
        if SIGMA_tab is not None:
            result = self._validate(Path(SIGMA_tab).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if disp_east is not None:
            Path(disp_east).touch()
        if disp_north is not None:
            Path(disp_north).touch()
        if disp_up is not None:
            Path(disp_up).touch()
        if sigma_east is not None:
            Path(sigma_east).touch()
        if sigma_north is not None:
            Path(sigma_north).touch()
        if sigma_up is not None:
            Path(sigma_up).touch()
        if chi2 is not None:
            Path(chi2).touch()
        return result

    def dem_coord(self, DEM_par: str, east: str, north: str, dtype = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "dem_coord", supplied_args))

        if "dem_coord" in self.call_count:
            self.call_count["dem_coord"] += 1
        else:
            self.call_count["dem_coord"] = 1

        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if east is not None:
            Path(east).touch()
        if north is not None:
            Path(north).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def rotate_image(self, data_in: str, width_in, angle, data_out: str, width_out, nlines_out, interp_mode = None, dtype = None, order = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "rotate_image", supplied_args))

        if "rotate_image" in self.call_count:
            self.call_count["rotate_image"] += 1
        else:
            self.call_count["rotate_image"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6, 7]
        result = self._validate(interp_mode in valid_values, result)
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        return result

    def geocode_back(self, data_in: str, width_in, lookup_table: str, data_out: str, width_out, nlines_out = None, interp_mode = None, dtype = None, lr_in = None, lr_out = None, order = None, e_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "geocode_back", supplied_args))

        if "geocode_back" in self.call_count:
            self.call_count["geocode_back"] += 1
        else:
            self.call_count["geocode_back"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if lookup_table is not None:
            result = self._validate(Path(lookup_table).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(e_flag in valid_values, result)
        return result

    def sarpix_coord(self, SLC_par: str, OFF_par: str = None, DEM_par: str = None, azlin = None, rpix = None, ref_hgt = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "sarpix_coord", supplied_args))

        if "sarpix_coord" in self.call_count:
            self.call_count["sarpix_coord"] += 1
        else:
            self.call_count["sarpix_coord"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        return result

    def resamp_image(self, data_in: str, width_in, xscale, yscale, data_out: str, width_out, nlines_out, interp_mode = None, dtype = None, order = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "resamp_image", supplied_args))

        if "resamp_image" in self.call_count:
            self.call_count["resamp_image"] += 1
        else:
            self.call_count["resamp_image"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        return result

    def pol2rec(self, data1: str, SLC_par1: str, data2: str, SLC_par2: str, pix_size: str, dtype, mode = None, xmin = None, nx = None, ymin = None, ny = None, rmax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "pol2rec", supplied_args))

        if "pol2rec" in self.call_count:
            self.call_count["pol2rec"] += 1
        else:
            self.call_count["pol2rec"] = 1

        if data1 is not None:
            result = self._validate(Path(data1).exists(), result)
        if SLC_par1 is not None:
            result = self._validate(Path(SLC_par1).exists(), result)
        if data2 is not None:
            Path(data2).touch()
        if SLC_par2 is not None:
            Path(SLC_par2).touch()
        if pix_size is not None:
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

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if DIFF_par is not None:
            result = self._validate(Path(DIFF_par).exists(), result)
        if SLC_2R is not None:
            Path(SLC_2R).touch()
        if SLC2R_par is not None:
            Path(SLC2R_par).touch()
        return result

    def gc_GPRI_map(self, MLI_par: str, DEM_par: str, DEM: str, DEM_seg_par: str, DEM_seg: str, lookup_table: str, lat_ovr = None, lon_ovr = None, sim_sar: str = None, lv_theta: str = None, lv_phi: str = None, u: str = None, v: str = None, inc: str = None, psi: str = None, pix: str = None, ls_map: str = None, frame = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DIFF", "gc_GPRI_map", supplied_args))

        if "gc_GPRI_map" in self.call_count:
            self.call_count["gc_GPRI_map"] += 1
        else:
            self.call_count["gc_GPRI_map"] = 1

        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if DEM_seg_par is not None and not Path(DEM_seg_par).exists():
            Path(DEM_seg_par).touch()
        if DEM_seg is not None:
            Path(DEM_seg).touch()
        if lookup_table is not None:
            Path(lookup_table).touch()
        if sim_sar is not None:
            Path(sim_sar).touch()
        if lv_theta is not None:
            Path(lv_theta).touch()
        if lv_phi is not None:
            Path(lv_phi).touch()
        if u is not None:
            Path(u).touch()
        if v is not None:
            Path(v).touch()
        if inc is not None:
            Path(inc).touch()
        if psi is not None:
            Path(psi).touch()
        if pix is not None:
            Path(pix).touch()
        if ls_map is not None:
            Path(ls_map).touch()
        return result

    def dop_mlcc(self, SAR_par: str, PROC_par: str, signal_data: str, plot_data: str = None, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "dop_mlcc", supplied_args))

        if "dop_mlcc" in self.call_count:
            self.call_count["dop_mlcc"] += 1
        else:
            self.call_count["dop_mlcc"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if plot_data is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
            Path(PROC_par).touch()
        return result

    def doppler_2d(self, SAR_par: str, PROC_par: str, signal_data: str, dop2d: str, loff = None, blsz = None, nbl = None, a2_flg = None, b0_flg = None, b1_flg = None, c0_flg = None, ambig_flag = None, namb = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "doppler_2d", supplied_args))

        if "doppler_2d" in self.call_count:
            self.call_count["doppler_2d"] += 1
        else:
            self.call_count["doppler_2d"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if dop2d is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
            Path(PROC_par).touch()
        return result

    def az_proc(self, SAR_par: str, PROC_par: str, rc_data: str, SLC: str, az_patch = None, SLC_format = None, cal_fact = None, SLC_type = None, kaiser = None, npatch = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "az_proc", supplied_args))

        if "az_proc" in self.call_count:
            self.call_count["az_proc"] += 1
        else:
            self.call_count["az_proc"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if rc_data is not None:
            result = self._validate(Path(rc_data).exists(), result)
        if SLC is not None:
            Path(SLC).touch()
        valid_values = [0, 1]
        result = self._validate(SLC_format in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(SLC_type in valid_values, result)
        return result

    def dop_interf(self, SAR_par1, PROC_par1, PROC_par2, PROC_par1_out: str, PROC_par2_out: str, dop: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "dop_interf", supplied_args))

        if "dop_interf" in self.call_count:
            self.call_count["dop_interf"] += 1
        else:
            self.call_count["dop_interf"] = 1

        if PROC_par1_out is not None:
            Path(PROC_par1_out).touch()
        if PROC_par2_out is not None:
            Path(PROC_par2_out).touch()
        if dop is not None:
            Path(dop).touch()
        return result

    def CS_proc(self, HDF5: str, SAR_par: str, PROC_par: str, raw_out: str, loff = None, nl = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "CS_proc", supplied_args))

        if "CS_proc" in self.call_count:
            self.call_count["CS_proc"] += 1
        else:
            self.call_count["CS_proc"] = 1

        if HDF5 is not None:
            result = self._validate(Path(HDF5).exists(), result)
        if SAR_par is not None:
            Path(SAR_par).touch()
        if PROC_par is not None:
            Path(PROC_par).touch()
        if raw_out is not None:
            Path(raw_out).touch()
        return result

    def azsp_SLC(self, SAR_par: str, PROC_par: str, SAR_data: str, spectrum, loff = None, roff = None, nsub = None, data_format = None, pltflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "azsp_SLC", supplied_args))

        if "azsp_SLC" in self.call_count:
            self.call_count["azsp_SLC"] += 1
        else:
            self.call_count["azsp_SLC"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if SAR_data is not None:
            result = self._validate(Path(SAR_data).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_format in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def pre_rc_JERS(self, SAR_par: str, PROC_par: str, rspec: str, signal_data: str, rc_data: str, prefilt_dec = None, kaiser = None, filt_lm = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "pre_rc_JERS", supplied_args))

        if "pre_rc_JERS" in self.call_count:
            self.call_count["pre_rc_JERS"] += 1
        else:
            self.call_count["pre_rc_JERS"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if rspec is not None:
            result = self._validate(Path(rspec).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if rc_data is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
            Path(PROC_par).touch()
        return result

    def PALSAR_proc(self, CEOS_SAR_leader: str, SAR_par: str, PROC_par: str, CEOS_raw_data: str, raw_out: str, TX_POL = None, RX_POL = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "PALSAR_proc", supplied_args))

        if "PALSAR_proc" in self.call_count:
            self.call_count["PALSAR_proc"] += 1
        else:
            self.call_count["PALSAR_proc"] = 1

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if SAR_par is not None:
            Path(SAR_par).touch()
        if PROC_par is not None:
            Path(PROC_par).touch()
        if CEOS_raw_data is not None:
            result = self._validate(Path(CEOS_raw_data).exists(), result)
        if raw_out is not None:
            Path(raw_out).touch()
        valid_values = [0, 1]
        result = self._validate(TX_POL in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(RX_POL in valid_values, result)
        return result

    def rspec_real(self, SAR_par: str, PROC_par: str, signal_data: str, range_spec: str, loff = None, nlspec = None, nrfft = None, pltflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "rspec_real", supplied_args))

        if "rspec_real" in self.call_count:
            self.call_count["rspec_real"] += 1
        else:
            self.call_count["rspec_real"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if range_spec is not None:
            Path(range_spec).touch()
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def rspec_IQ(self, SAR_par: str, PROC_par: str, signal_data: str, range_spec: str, loff = None, nlspec = None, nrfft = None, roff = None, pltflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "rspec_IQ", supplied_args))

        if "rspec_IQ" in self.call_count:
            self.call_count["rspec_IQ"] += 1
        else:
            self.call_count["rspec_IQ"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if range_spec is not None:
            Path(range_spec).touch()
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def af(self, SAR_par: str, PROC_par: str, SLC: str, rwin = None, azwin = None, dr = None, daz = None, thres = None, update_flg = None, a1_flg = None, b0_flg = None, offsets: str = None, dac_flg = None, n_ovr = None, roff = None, azoff = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "af", supplied_args))

        if "af" in self.call_count:
            self.call_count["af"] += 1
        else:
            self.call_count["af"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        valid_values = [0, 1]
        result = self._validate(update_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(a1_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(b0_flg in valid_values, result)
        if offsets is not None:
            Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(dac_flg in valid_values, result)
        return result

    def hist_IQ(self, SAR_par: str, PROC_par: str, signal_data: str, historgram: str, loff = None, nl = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "hist_IQ", supplied_args))

        if "hist_IQ" in self.call_count:
            self.call_count["hist_IQ"] += 1
        else:
            self.call_count["hist_IQ"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if historgram is not None:
            Path(historgram).touch()
        return result

    def PALSAR_burst_sync(self, SAR_par1: str, PROC_par1: str, raw1: str, SAR_par2: str, PROC_par2: str, raw2: str, PROC_par1_out: str, raw1_out: str, PROC_par2_out: str, raw2_out: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "PALSAR_burst_sync", supplied_args))

        if "PALSAR_burst_sync" in self.call_count:
            self.call_count["PALSAR_burst_sync"] += 1
        else:
            self.call_count["PALSAR_burst_sync"] = 1

        if SAR_par1 is not None:
            result = self._validate(Path(SAR_par1).exists(), result)
        if PROC_par1 is not None:
            result = self._validate(Path(PROC_par1).exists(), result)
        if raw1 is not None:
            result = self._validate(Path(raw1).exists(), result)
        if SAR_par2 is not None:
            result = self._validate(Path(SAR_par2).exists(), result)
        if PROC_par2 is not None:
            result = self._validate(Path(PROC_par2).exists(), result)
        if raw2 is not None:
            result = self._validate(Path(raw2).exists(), result)
        if PROC_par1_out is not None:
            Path(PROC_par1_out).touch()
        if raw1_out is not None:
            Path(raw1_out).touch()
        if PROC_par2_out is not None:
            Path(PROC_par2_out).touch()
        if raw2_out is not None:
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

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if raw_IQ is not None:
            result = self._validate(Path(raw_IQ).exists(), result)
        if raw_IQ_swap is not None:
            Path(raw_IQ_swap).touch()
        return result

    def dop_ambig(self, SAR_par: str, PROC_par: str, signal_data: str, algorithm = None, loff = None, output_plot: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "dop_ambig", supplied_args))

        if "dop_ambig" in self.call_count:
            self.call_count["dop_ambig"] += 1
        else:
            self.call_count["dop_ambig"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        valid_values = [1, 2]
        result = self._validate(algorithm in valid_values, result)
        if output_plot is not None:
            Path(output_plot).touch()
        return result

    def PRC_proc(self, PROC_par: str, PRC, nstate = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "PRC_proc", supplied_args))

        if "PRC_proc" in self.call_count:
            self.call_count["PRC_proc"] += 1
        else:
            self.call_count["PRC_proc"] = 1

        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        return result

    def rc_fmcw(self, SAR_par: str, PROC_par: str, signal_data: str, rc_data: str, nrc_off, nrc_samp = None, loff = None, nl = None, kaiser = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "rc_fmcw", supplied_args))

        if "rc_fmcw" in self.call_count:
            self.call_count["rc_fmcw"] += 1
        else:
            self.call_count["rc_fmcw"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if rc_data is not None:
            Path(rc_data).touch()
        return result

    def RSAT_lks(self, SLC_PROC_par: str, MLI_PROC_par: str, SLC_image: str, ML_image, kaiser = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "RSAT_lks", supplied_args))

        if "RSAT_lks" in self.call_count:
            self.call_count["RSAT_lks"] += 1
        else:
            self.call_count["RSAT_lks"] = 1

        if SLC_PROC_par is not None:
            result = self._validate(Path(SLC_PROC_par).exists(), result)
        if MLI_PROC_par is not None:
            Path(MLI_PROC_par).touch()
        if SLC_image is not None:
            result = self._validate(Path(SLC_image).exists(), result)
        return result

    def extract_psd(self, spectra, num, output_spectrum):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "extract_psd", supplied_args))

        if "extract_psd" in self.call_count:
            self.call_count["extract_psd"] += 1
        else:
            self.call_count["extract_psd"] = 1

        return result

    def doppler_real(self, SAR_par: str, PROC_par: str, signal_data: str, doppler: str, loff = None, nsub = None, ambig_flag = None, namb = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "doppler_real", supplied_args))

        if "doppler_real" in self.call_count:
            self.call_count["doppler_real"] += 1
        else:
            self.call_count["doppler_real"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if doppler is not None:
            Path(doppler).touch()
        return result

    def pre_rc(self, SAR_par: str, PROC_par: str, signal_data: str, rc_data: str, prefilt_dec = None, loff = None, nl = None, nr_samp = None, kaiser = None, filt_lm = None, nr_ext = None, fr_ext = None, pre_ext = None, post_ext = None, RFI_filt = None, RFI_thres = None, fc_offset = None, win_bw = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "pre_rc", supplied_args))

        if "pre_rc" in self.call_count:
            self.call_count["pre_rc"] += 1
        else:
            self.call_count["pre_rc"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if rc_data is not None:
            Path(rc_data).touch()
        valid_values = [0, 1]
        result = self._validate(RFI_filt in valid_values, result)
        return result

    def PALSAR_proc_WB(self, CEOS_SAR_leader: str, SAR_par: str, PROC_par: str, CEOS_raw_data: str, beam: str, raw_out: str, prf: str = None, wflg: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "PALSAR_proc_WB", supplied_args))

        if "PALSAR_proc_WB" in self.call_count:
            self.call_count["PALSAR_proc_WB"] += 1
        else:
            self.call_count["PALSAR_proc_WB"] = 1

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if SAR_par is not None:
            Path(SAR_par).touch()
        if PROC_par is not None:
            Path(PROC_par).touch()
        if CEOS_raw_data is not None:
            result = self._validate(Path(CEOS_raw_data).exists(), result)
        if beam is not None:
            result = self._validate(Path(beam).exists(), result)
        if raw_out is not None:
            Path(raw_out).touch()
        if prf is not None:
            result = self._validate(Path(prf).exists(), result)
        if wflg is not None:
            result = self._validate(Path(wflg).exists(), result)
        return result

    def ORRM_proc(self, PROC_par: str, ORRM: str, nstate = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ORRM_proc", supplied_args))

        if "ORRM_proc" in self.call_count:
            self.call_count["ORRM_proc"] += 1
        else:
            self.call_count["ORRM_proc"] = 1

        if PROC_par is not None and not Path(PROC_par).exists():
            Path(PROC_par).touch()
        if ORRM is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
            Path(PROC_par).touch()
        return result

    def RSAT_raw(self, CEOS_leader: str, SAR_par: str, PROC_par: str, raw_data_files: str = None, raw_out: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "RSAT_raw", supplied_args))

        if "RSAT_raw" in self.call_count:
            self.call_count["RSAT_raw"] += 1
        else:
            self.call_count["RSAT_raw"] = 1

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if SAR_par is not None:
            Path(SAR_par).touch()
        if PROC_par is not None:
            Path(PROC_par).touch()
        if raw_data_files is not None:
            result = self._validate(Path(raw_data_files).exists(), result)
        if raw_out is not None:
            Path(raw_out).touch()
        return result

    def ERS_ENVISAT_proc(self, L0: str, SAR_par: str, PROC_par: str, raw: str, loff = None, nl = None, swst_flg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_ENVISAT_proc", supplied_args))

        if "ERS_ENVISAT_proc" in self.call_count:
            self.call_count["ERS_ENVISAT_proc"] += 1
        else:
            self.call_count["ERS_ENVISAT_proc"] = 1

        if L0 is not None:
            result = self._validate(Path(L0).exists(), result)
        if SAR_par is not None and not Path(SAR_par).exists():
            Path(SAR_par).touch()
        if PROC_par is not None:
            Path(PROC_par).touch()
        if raw is not None:
            Path(raw).touch()
        valid_values = [0, 1]
        result = self._validate(swst_flg in valid_values, result)
        return result

    def multi_SLC(self, SLC_PROC_par: str, MLI_PROC_par: str, SLC: str, MLI: str, rlks, azlks, slc_format = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "multi_SLC", supplied_args))

        if "multi_SLC" in self.call_count:
            self.call_count["multi_SLC"] += 1
        else:
            self.call_count["multi_SLC"] = 1

        if SLC_PROC_par is not None:
            result = self._validate(Path(SLC_PROC_par).exists(), result)
        if MLI_PROC_par is not None:
            Path(MLI_PROC_par).touch()
        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if MLI is not None:
            Path(MLI).touch()
        valid_values = [0, 1]
        result = self._validate(slc_format in valid_values, result)
        return result

    def JERS_acs(self, USER_HEADER: str, SEG_DESCR: str, ORBIT_DATA: str, SENSOR_DATA: str, track: str, SAR_par: str, PROC_par: str, raw_out: str, loff, nl, nsx = None, fsx = None, deskew = None, terra_alt = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "JERS_acs", supplied_args))

        if "JERS_acs" in self.call_count:
            self.call_count["JERS_acs"] += 1
        else:
            self.call_count["JERS_acs"] = 1

        if USER_HEADER is not None:
            result = self._validate(Path(USER_HEADER).exists(), result)
        if SEG_DESCR is not None:
            result = self._validate(Path(SEG_DESCR).exists(), result)
        if ORBIT_DATA is not None:
            result = self._validate(Path(ORBIT_DATA).exists(), result)
        if SENSOR_DATA is not None:
            result = self._validate(Path(SENSOR_DATA).exists(), result)
        if track is not None:
            Path(track).touch()
        if SAR_par is not None:
            Path(SAR_par).touch()
        if PROC_par is not None:
            Path(PROC_par).touch()
        if raw_out is not None:
            Path(raw_out).touch()
        return result

    def DORIS_proc(self, PROC_par: str, DOR: str, nstate = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "DORIS_proc", supplied_args))

        if "DORIS_proc" in self.call_count:
            self.call_count["DORIS_proc"] += 1
        else:
            self.call_count["DORIS_proc"] = 1

        if PROC_par is not None and not Path(PROC_par).exists():
            Path(PROC_par).touch()
        if DOR is not None:
            result = self._validate(Path(DOR).exists(), result)
        return result

    def ASAR_AP_proc(self, L0: str, INS: str, SAR_par1: str, SAR_par2: str, PROC_par1: str, PROC_par2: str, raw1: str, raw2: str, ant_gain1: str, ant_gain2: str, loff = None, nl = None, roff = None, nr = None, refer = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ASAR_AP_proc", supplied_args))

        if "ASAR_AP_proc" in self.call_count:
            self.call_count["ASAR_AP_proc"] += 1
        else:
            self.call_count["ASAR_AP_proc"] = 1

        if L0 is not None:
            result = self._validate(Path(L0).exists(), result)
        if INS is not None:
            result = self._validate(Path(INS).exists(), result)
        if SAR_par1 is not None:
            Path(SAR_par1).touch()
        if SAR_par2 is not None:
            Path(SAR_par2).touch()
        if PROC_par1 is not None:
            Path(PROC_par1).touch()
        if PROC_par2 is not None:
            Path(PROC_par2).touch()
        if raw1 is not None:
            Path(raw1).touch()
        if raw2 is not None:
            Path(raw2).touch()
        if ant_gain1 is not None:
            result = self._validate(Path(ant_gain1).exists(), result)
        if ant_gain2 is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
            Path(PROC_par).touch()
        return result

    def rc_real(self, SAR_par: str, PROC_par: str, signal_data: str, rc_data: str, loff = None, nl = None, kaiser = None, nr_ext = None, fr_ext = None, r_chirp: str = None, rfi_filt = None, rfi_thres = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "rc_real", supplied_args))

        if "rc_real" in self.call_count:
            self.call_count["rc_real"] += 1
        else:
            self.call_count["rc_real"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if rc_data is not None:
            Path(rc_data).touch()
        if r_chirp is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
            Path(PROC_par).touch()
        return result

    def prefilt(self, SAR_par: str, PROC_par: str, rc_data: str, prefilt_out: str, prefilt_dec, filt_lm = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "prefilt", supplied_args))

        if "prefilt" in self.call_count:
            self.call_count["prefilt"] += 1
        else:
            self.call_count["prefilt"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if rc_data is not None:
            result = self._validate(Path(rc_data).exists(), result)
        if prefilt_out is not None:
            Path(prefilt_out).touch()
        return result

    def doppler(self, SAR_par: str, PROC_par: str, signal_data: str, doppler: str, loff = None, nsub = None, ambig_flag = None, namb = None, order = None, pltflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "doppler", supplied_args))

        if "doppler" in self.call_count:
            self.call_count["doppler"] += 1
        else:
            self.call_count["doppler"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if doppler is not None:
            Path(doppler).touch()
        valid_values = [0, 1, 2]
        result = self._validate(ambig_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def DELFT_proc2(self, PROC_par: str, DELFT_dir, nstate = None, interval = None, ODR = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "DELFT_proc2", supplied_args))

        if "DELFT_proc2" in self.call_count:
            self.call_count["DELFT_proc2"] += 1
        else:
            self.call_count["DELFT_proc2"] = 1

        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        return result

    def copy(self, infile: str, outfile: str, lbytes, start = None, nlines = None, offset = None, file_ldr = None, offb = None, nbyte = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "copy", supplied_args))

        if "copy" in self.call_count:
            self.call_count["copy"] += 1
        else:
            self.call_count["copy"] = 1

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if CEOS_trailer is not None:
            result = self._validate(Path(CEOS_trailer).exists(), result)
        if PROC_par is not None:
            Path(PROC_par).touch()
        return result

    def ORB_prop(self, PROC_par: str, nstate = None, interval = None, extra = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ORB_prop", supplied_args))

        if "ORB_prop" in self.call_count:
            self.call_count["ORB_prop"] += 1
        else:
            self.call_count["ORB_prop"] = 1

        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        return result

    def ERS_fix(self, ERS_PAF, SAR_par: str, PROC_par: str, cc_flag, raw = None, output_file: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_fix", supplied_args))

        if "ERS_fix" in self.call_count:
            self.call_count["ERS_fix"] += 1
        else:
            self.call_count["ERS_fix"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if output_file is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
            Path(PROC_par).touch()
        return result

    def ptarg(self, SLC: str, width: str, r_samp: str, az_samp: str, ptr_image: str, r_plot: str, az_plot: str, data_format, win = None, pltflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ptarg", supplied_args))

        if "ptarg" in self.call_count:
            self.call_count["ptarg"] += 1
        else:
            self.call_count["ptarg"] = 1

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if width is not None:
            result = self._validate(Path(width).exists(), result)
        if r_samp is not None:
            result = self._validate(Path(r_samp).exists(), result)
        if az_samp is not None:
            result = self._validate(Path(az_samp).exists(), result)
        if ptr_image is not None:
            Path(ptr_image).touch()
        if r_plot is not None:
            Path(r_plot).touch()
        if az_plot is not None:
            Path(az_plot).touch()
        valid_values = [0, 1]
        result = self._validate(data_format in valid_values, result)
        valid_values = [0, 1, 2]
        result = self._validate(pltflg in valid_values, result)
        return result

    def multi_GRD_SLC(self, SLC_PROC_par: str, GRD_PROC_par: str, SLC_image: str, GRD_image: str, rlks, azlks, interp_mode = None, sample_spacing = None, gr_start = None, t_start = None, t_end = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "multi_GRD_SLC", supplied_args))

        if "multi_GRD_SLC" in self.call_count:
            self.call_count["multi_GRD_SLC"] += 1
        else:
            self.call_count["multi_GRD_SLC"] = 1

        if SLC_PROC_par is not None:
            result = self._validate(Path(SLC_PROC_par).exists(), result)
        if GRD_PROC_par is not None:
            Path(GRD_PROC_par).touch()
        if SLC_image is not None:
            result = self._validate(Path(SLC_image).exists(), result)
        if GRD_image is not None:
            Path(GRD_image).touch()
        valid_values = [0, 1]
        result = self._validate(interp_mode in valid_values, result)
        return result

    def create_sar_par(self, SAR_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "create_sar_par", supplied_args))

        if "create_sar_par" in self.call_count:
            self.call_count["create_sar_par"] += 1
        else:
            self.call_count["create_sar_par"] = 1

        if SAR_par is not None and not Path(SAR_par).exists():
            Path(SAR_par).touch()
        return result

    def rspec_JERS(self, SAR_par: str, PROC_par: str, signal_data: str, range_spec: str, nr_samp = None, nl_spec = None, loff = None, nlines = None, nr_ext = None, fr_ext = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "rspec_JERS", supplied_args))

        if "rspec_JERS" in self.call_count:
            self.call_count["rspec_JERS"] += 1
        else:
            self.call_count["rspec_JERS"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if range_spec is not None:
            Path(range_spec).touch()
        return result

    def azsp_IQ(self, SAR_par: str, PROC_par: str, signal_data: str, spectrum: str, loff = None, roff = None, nsub = None, ambig_flg = None, namb = None, pltflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "azsp_IQ", supplied_args))

        if "azsp_IQ" in self.call_count:
            self.call_count["azsp_IQ"] += 1
        else:
            self.call_count["azsp_IQ"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if spectrum is not None:
            Path(spectrum).touch()
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def ERS_proc_ACRES(self, CEOS_SAR_leader: str, PROC_par: str, type = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ERS_proc_ACRES", supplied_args))

        if "ERS_proc_ACRES" in self.call_count:
            self.call_count["ERS_proc_ACRES"] += 1
        else:
            self.call_count["ERS_proc_ACRES"] = 1

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
            Path(PROC_par).touch()
        valid_values = [0, 1]
        result = self._validate(type in valid_values, result)
        return result

    def pre_rc_RSAT(self, SAR_par: str, PROC_par: str, signal_data: str, rc_data: str, prefilt_dec = None, loff = None, nl = None, nr_samp = None, kaiser = None, filt_lm = None, nr_ext = None, fr_ext = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "pre_rc_RSAT", supplied_args))

        if "pre_rc_RSAT" in self.call_count:
            self.call_count["pre_rc_RSAT"] += 1
        else:
            self.call_count["pre_rc_RSAT"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if signal_data is not None:
            result = self._validate(Path(signal_data).exists(), result)
        if rc_data is not None:
            Path(rc_data).touch()
        return result

    def autof(self, SAR_par: str, PROC_par: str, rc_data: str, autofocus: str, SNR_min = None, prefilter = None, auto_az = None, az_offset = None, auto_bins = None, dop_ambig = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "autof", supplied_args))

        if "autof" in self.call_count:
            self.call_count["autof"] += 1
        else:
            self.call_count["autof"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if rc_data is not None:
            result = self._validate(Path(rc_data).exists(), result)
        if autofocus is not None:
            Path(autofocus).touch()
        valid_values = [0, 1]
        result = self._validate(dop_ambig in valid_values, result)
        return result

    def cat_raw(self, RAW_list: str, SAR_par: str, PROC_par: str, RAW_out: str, fill = None, loff = None, nl = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "cat_raw", supplied_args))

        if "cat_raw" in self.call_count:
            self.call_count["cat_raw"] += 1
        else:
            self.call_count["cat_raw"] = 1

        if RAW_list is not None:
            result = self._validate(Path(RAW_list).exists(), result)
        if SAR_par is not None:
            Path(SAR_par).touch()
        if PROC_par is not None:
            Path(PROC_par).touch()
        if RAW_out is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PROC_par is not None:
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

        if CEOS_SAR_ldr is not None:
            result = self._validate(Path(CEOS_SAR_ldr).exists(), result)
        if PROC_par is not None:
            Path(PROC_par).touch()
        return result

    def create_proc_par(self, SAR_par: str, PROC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "create_proc_par", supplied_args))

        if "create_proc_par" in self.call_count:
            self.call_count["create_proc_par"] += 1
        else:
            self.call_count["create_proc_par"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None and not Path(PROC_par).exists():
            Path(PROC_par).touch()
        return result

    def SIRC_proc(self, CEOS_SAR_leader: str, SAR_par: str, PROC_par: str, UTC_MET = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "SIRC_proc", supplied_args))

        if "SIRC_proc" in self.call_count:
            self.call_count["SIRC_proc"] += 1
        else:
            self.call_count["SIRC_proc"] = 1

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if SAR_par is not None:
            Path(SAR_par).touch()
        if PROC_par is not None:
            Path(PROC_par).touch()
        return result

    def ASAR_IM_proc(self, L0: str, INS: str, SAR_par: str, PROC_par: str, raw: str, ant_gain: str, loff = None, nl = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("MSP", "ASAR_IM_proc", supplied_args))

        if "ASAR_IM_proc" in self.call_count:
            self.call_count["ASAR_IM_proc"] += 1
        else:
            self.call_count["ASAR_IM_proc"] = 1

        if L0 is not None:
            result = self._validate(Path(L0).exists(), result)
        if INS is not None:
            result = self._validate(Path(INS).exists(), result)
        if SAR_par is not None:
            Path(SAR_par).touch()
        if PROC_par is not None:
            Path(PROC_par).touch()
        if raw is not None:
            Path(raw).touch()
        if ant_gain is not None:
            result = self._validate(Path(ant_gain).exists(), result)
        return result

    def dishgt(self, hgt: str, pwr: str, width, start_hgt = None, start_pwr = None, nlines = None, m_cycle = None, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dishgt", supplied_args))

        if "dishgt" in self.call_count:
            self.call_count["dishgt"] += 1
        else:
            self.call_count["dishgt"] = 1

        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        return result

    def rasdt_pwr(self, data: str, pwr: str, width, start_data = None, start_pwr = None, nlines = None, pixavr = None, pixavaz = None, cycle = None, scale = None, exp = None, LR = None, rasf: str = None, cc: str = None, start_cc = None, cc_min = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasdt_pwr", supplied_args))

        if "rasdt_pwr" in self.call_count:
            self.call_count["rasdt_pwr"] += 1
        else:
            self.call_count["rasdt_pwr"] = 1

        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        if cc is not None:
            result = self._validate(Path(cc).exists(), result)
        return result

    def rasdt_cmap(self, data: str, pwr: str, width, start_data = None, start_pwr = None, nlines = None, pixavr = None, pixavaz = None, min = None, max = None, mflg = None, cmap = None, scale = None, exp = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasdt_cmap", supplied_args))

        if "rasdt_cmap" in self.call_count:
            self.call_count["rasdt_cmap"] += 1
        else:
            self.call_count["rasdt_cmap"] = 1

        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        valid_values = [0, 1]
        result = self._validate(mflg in valid_values, result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def dis2hgt(self, hgt1: str, hgt2: str, width1, width2, start = None, nlines = None, roff = None, azoff = None, m_cycle = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2hgt", supplied_args))

        if "dis2hgt" in self.call_count:
            self.call_count["dis2hgt"] += 1
        else:
            self.call_count["dis2hgt"] = 1

        if hgt1 is not None:
            result = self._validate(Path(hgt1).exists(), result)
        if hgt2 is not None:
            result = self._validate(Path(hgt2).exists(), result)
        return result

    def gcp_ras(self, ras: str, GCP: str, mag = None, win_sz = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "gcp_ras", supplied_args))

        if "gcp_ras" in self.call_count:
            self.call_count["gcp_ras"] += 1
        else:
            self.call_count["gcp_ras"] = 1

        if ras is not None:
            result = self._validate(Path(ras).exists(), result)
        if GCP is not None:
            Path(GCP).touch()
        return result

    def cpx_math(self, d1: str, d2: str, d_out: str, width, mode, roff = None, loff = None, nr = None, nl = None, c_re = None, c_im = None, zflg = None, rflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "cpx_math", supplied_args))

        if "cpx_math" in self.call_count:
            self.call_count["cpx_math"] += 1
        else:
            self.call_count["cpx_math"] = 1

        if d1 is not None:
            result = self._validate(Path(d1).exists(), result)
        if d2 is not None:
            result = self._validate(Path(d2).exists(), result)
        if d_out is not None:
            Path(d_out).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def disflag(self, flag: str, width, start = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disflag", supplied_args))

        if "disflag" in self.call_count:
            self.call_count["disflag"] += 1
        else:
            self.call_count["disflag"] = 1

        if flag is not None:
            result = self._validate(Path(flag).exists(), result)
        return result

    def set_value(self, PAR_in: str, PAR_out: str, keyword, value, new_key = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "set_value", supplied_args))

        if "set_value" in self.call_count:
            self.call_count["set_value"] += 1
        else:
            self.call_count["set_value"] = 1

        if PAR_in is not None:
            result = self._validate(Path(PAR_in).exists(), result)
        if PAR_out is not None:
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

        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        if real is not None:
            Path(real).touch()
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(type in valid_values, result)
        return result

    def rascc(self, cc: str, pwr: str, width, start_cc = None, start_pwr = None, nlines = None, pixavr = None, pixavaz = None, cmin = None, cmax = None, scale = None, exp = None, LR = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rascc", supplied_args))

        if "rascc" in self.call_count:
            self.call_count["rascc"] += 1
        else:
            self.call_count["rascc"] = 1

        if cc is not None:
            result = self._validate(Path(cc).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def disgbyte(self, image: str, width, start = None, nlines = None, scale = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disgbyte", supplied_args))

        if "disgbyte" in self.call_count:
            self.call_count["disgbyte"] += 1
        else:
            self.call_count["disgbyte"] = 1

        if image is not None:
            result = self._validate(Path(image).exists(), result)
        return result

    def dis_dB(self, pwr: str, width, start = None, nlines = None, min_dB = None, max_dB = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis_dB", supplied_args))

        if "dis_dB" in self.call_count:
            self.call_count["dis_dB"] += 1
        else:
            self.call_count["dis_dB"] = 1

        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        return result

    def discc(self, cc: str, pwr: str, width, start_cc = None, start_pwr = None, nlines = None, min_corr = None, max_corr = None, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "discc", supplied_args))

        if "discc" in self.call_count:
            self.call_count["discc"] += 1
        else:
            self.call_count["discc"] = 1

        if cc is not None:
            result = self._validate(Path(cc).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        return result

    def svg_map(self, image: str, dem_par: str, svg: str, font = None, fsize = None, color = None, gcolor = None, majorx = None, majory = None, minorx = None, minory = None, thick = None, grid = None, gopac = None, gdash = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "svg_map", supplied_args))

        if "svg_map" in self.call_count:
            self.call_count["svg_map"] += 1
        else:
            self.call_count["svg_map"] = 1

        if image is not None:
            result = self._validate(Path(image).exists(), result)
        if dem_par is not None:
            result = self._validate(Path(dem_par).exists(), result)
        if svg is not None:
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

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
            Path(outfile).touch()
        return result

    def real_to_cpx(self, data1: str, data2: str = None, cpx: str = None, width = None, type = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "real_to_cpx", supplied_args))

        if "real_to_cpx" in self.call_count:
            self.call_count["real_to_cpx"] += 1
        else:
            self.call_count["real_to_cpx"] = 1

        if data1 is not None:
            result = self._validate(Path(data1).exists(), result)
        if data2 is not None:
            result = self._validate(Path(data2).exists(), result)
        if cpx is not None:
            Path(cpx).touch()
        return result

    def rasmph_pwr(self, cpx: str, pwr: str, width, start_cpx = None, start_pwr = None, nlines = None, pixavr = None, pixavaz = None, scale = None, exp = None, LR = None, rasf: str = None, cc: str = None, start_cc = None, cc_min = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasmph_pwr", supplied_args))

        if "rasmph_pwr" in self.call_count:
            self.call_count["rasmph_pwr"] += 1
        else:
            self.call_count["rasmph_pwr"] = 1

        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        if cc is not None:
            result = self._validate(Path(cc).exists(), result)
        return result

    def data2geotiff(self, DEM_par: str, data: str, type, GeoTIFF: str, nodata = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "data2geotiff", supplied_args))

        if "data2geotiff" in self.call_count:
            self.call_count["data2geotiff"] += 1
        else:
            self.call_count["data2geotiff"] = 1

        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if data is not None:
            result = self._validate(Path(data).exists(), result)
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(type in valid_values, result)
        if GeoTIFF is not None:
            Path(GeoTIFF).touch()
        return result

    def disras(self, ras: str, mag = None, win_sz = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disras", supplied_args))

        if "disras" in self.call_count:
            self.call_count["disras"] += 1
        else:
            self.call_count["disras"] = 1

        if ras is not None:
            result = self._validate(Path(ras).exists(), result)
        return result

    def disrmg(self, unw: str, pwr: str, width, start_unw = None, start_pwr = None, nlines = None, ph_scale = None, scale = None, exp = None, ph_offset = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disrmg", supplied_args))

        if "disrmg" in self.call_count:
            self.call_count["disrmg"] += 1
        else:
            self.call_count["disrmg"] = 1

        if unw is not None:
            result = self._validate(Path(unw).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        return result

    def dis2rmg(self, unw1: str, unw2: str, width1, width2, start = None, nlines = None, roff = None, azoff = None, ph_scale = None, ph_offset = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2rmg", supplied_args))

        if "dis2rmg" in self.call_count:
            self.call_count["dis2rmg"] += 1
        else:
            self.call_count["dis2rmg"] = 1

        if unw1 is not None:
            result = self._validate(Path(unw1).exists(), result)
        if unw2 is not None:
            result = self._validate(Path(unw2).exists(), result)
        return result

    def dis2ras(self, ras1: str, ras2: str, mag = None, win_sz = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2ras", supplied_args))

        if "dis2ras" in self.call_count:
            self.call_count["dis2ras"] += 1
        else:
            self.call_count["dis2ras"] = 1

        if ras1 is not None:
            result = self._validate(Path(ras1).exists(), result)
        if ras2 is not None:
            result = self._validate(Path(ras2).exists(), result)
        return result

    def dismph_pwr(self, cpx: str, pwr: str, width, start_cpx = None, start_pwr = None, nlines = None, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph_pwr", supplied_args))

        if "dismph_pwr" in self.call_count:
            self.call_count["dismph_pwr"] += 1
        else:
            self.call_count["dismph_pwr"] = 1

        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        return result

    def vec_math(self, d1: str, d2: str, d_out: str, width, mode, c1 = None, c2 = None, c3 = None, nflg = None, nl = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "vec_math", supplied_args))

        if "vec_math" in self.call_count:
            self.call_count["vec_math"] += 1
        else:
            self.call_count["vec_math"] = 1

        if d1 is not None:
            result = self._validate(Path(d1).exists(), result)
        if d2 is not None:
            result = self._validate(Path(d2).exists(), result)
        if d_out is not None:
            Path(d_out).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6]
        result = self._validate(mode in valid_values, result)
        return result

    def dismph_fft(self, cpx: str, width, start = None, nlines = None, scale = None, exp = None, nfft = None, mag = None, data_type = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph_fft", supplied_args))

        if "dismph_fft" in self.call_count:
            self.call_count["dismph_fft"] += 1
        else:
            self.call_count["dismph_fft"] = 1

        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def svg_arrow(self, dv_norm: str, dv_phi: str, width, svg: str, image = None, norm = None, gridx = None, gridy = None, color = None, thick = None, head = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "svg_arrow", supplied_args))

        if "svg_arrow" in self.call_count:
            self.call_count["svg_arrow"] += 1
        else:
            self.call_count["svg_arrow"] = 1

        if dv_norm is not None:
            result = self._validate(Path(dv_norm).exists(), result)
        if dv_phi is not None:
            result = self._validate(Path(dv_phi).exists(), result)
        if svg is not None:
            Path(svg).touch()
        return result

    def ras3pwr(self, d1: str, d2: str, d3: str, width, start = None, nlines = None, pixavr = None, pixavaz = None, scale1 = None, scale2 = None, scale3 = None, exp = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras3pwr", supplied_args))

        if "ras3pwr" in self.call_count:
            self.call_count["ras3pwr"] += 1
        else:
            self.call_count["ras3pwr"] = 1

        if d1 is not None:
            result = self._validate(Path(d1).exists(), result)
        if d2 is not None:
            result = self._validate(Path(d2).exists(), result)
        if d3 is not None:
            result = self._validate(Path(d3).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def disdt_pwr(self, data: str, pwr: str, width, start_data = None, start_pwr = None, nlines = None, cycle = None, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disdt_pwr", supplied_args))

        if "disdt_pwr" in self.call_count:
            self.call_count["disdt_pwr"] += 1
        else:
            self.call_count["disdt_pwr"] = 1

        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if pwr is not None:
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

        if cmp1 is not None:
            result = self._validate(Path(cmp1).exists(), result)
        if cmp2 is not None:
            result = self._validate(Path(cmp2).exists(), result)
        if cmp3 is not None:
            result = self._validate(Path(cmp3).exists(), result)
        if vec is not None:
            Path(vec).touch()
        return result

    def float2short(self, infile: str, outfile: str, scale = None, exp = None, neg = None, output = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "float2short", supplied_args))

        if "float2short" in self.call_count:
            self.call_count["float2short"] += 1
        else:
            self.call_count["float2short"] = 1

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
            Path(outfile).touch()
        return result

    def disdt_pwr24(self, data: str, pwr: str, width, start_data = None, start_pwr = None, nlines = None, cycle = None, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disdt_pwr24", supplied_args))

        if "disdt_pwr24" in self.call_count:
            self.call_count["disdt_pwr24"] += 1
        else:
            self.call_count["disdt_pwr24"] = 1

        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        return result

    def rashgt_shd(self, hgt: str, data: str, width, col_post, row_post = None, start = None, nlines = None, pixavr = None, pixavaz = None, theta0 = None, phi0 = None, color0 = None, cycle = None, LR = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rashgt_shd", supplied_args))

        if "rashgt_shd" in self.call_count:
            self.call_count["rashgt_shd"] += 1
        else:
            self.call_count["rashgt_shd"] = 1

        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def ras2ras(self, ras_in: str, ras_out: str, cmap = None, force24 = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras2ras", supplied_args))

        if "ras2ras" in self.call_count:
            self.call_count["ras2ras"] += 1
        else:
            self.call_count["ras2ras"] = 1

        if ras_in is not None:
            result = self._validate(Path(ras_in).exists(), result)
        if ras_out is not None:
            Path(ras_out).touch()
        valid_values = [0, 1]
        result = self._validate(force24 in valid_values, result)
        return result

    def ras24_float(self, f1: str, f2: str, f3: str, width, rasf: str, color_model = None, h0 = None, hrange = None, imin = None, imax = None, sat_min = None, sat_max = None, sc1 = None, A1 = None, B1 = None, cyclic1 = None, sc2 = None, A2 = None, B2 = None, start_f1 = None, start_f2 = None, nlines = None, pixavr = None, pixavaz = None, LR = None, General = None, start_f3 = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras24_float", supplied_args))

        if "ras24_float" in self.call_count:
            self.call_count["ras24_float"] += 1
        else:
            self.call_count["ras24_float"] = 1

        if f1 is not None:
            result = self._validate(Path(f1).exists(), result)
        if f2 is not None:
            result = self._validate(Path(f2).exists(), result)
        if f3 is not None:
            result = self._validate(Path(f3).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        valid_values = [0, 1, 2]
        result = self._validate(color_model in valid_values, result)
        return result

    def double2float(self, infile: str, outfile: str, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "double2float", supplied_args))

        if "double2float" in self.call_count:
            self.call_count["double2float"] += 1
        else:
            self.call_count["double2float"] = 1

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
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

    def float_math(self, d1: str, d2: str, d_out: str, width, mode, roff = None, loff = None, nr = None, nl = None, c0 = None, zflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "float_math", supplied_args))

        if "float_math" in self.call_count:
            self.call_count["float_math"] += 1
        else:
            self.call_count["float_math"] = 1

        if d1 is not None:
            result = self._validate(Path(d1).exists(), result)
        if d2 is not None:
            result = self._validate(Path(d2).exists(), result)
        if d_out is not None:
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

        if din is not None:
            result = self._validate(Path(din).exists(), result)
        if dout is not None:
            Path(dout).touch()
        valid_values = [0, 1, 2]
        result = self._validate(dtype in valid_values, result)
        return result

    def disras_dem_par(self, ras: str, DEM_par: str, mag = None, win_sz = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disras_dem_par", supplied_args))

        if "disras_dem_par" in self.call_count:
            self.call_count["disras_dem_par"] += 1
        else:
            self.call_count["disras_dem_par"] = 1

        if ras is not None:
            result = self._validate(Path(ras).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        return result

    def dismph(self, cpx: str, width, start = None, nlines = None, scale = None, exp = None, dtype = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph", supplied_args))

        if "dismph" in self.call_count:
            self.call_count["dismph"] += 1
        else:
            self.call_count["dismph"] = 1

        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def rasshd(self, DEM: str, width, col_post, row_post = None, start = None, nlines = None, pixavr = None, pixavaz = None, theta0 = None, phi0 = None, LR = None, rasf: str = None, dtype = None, zero_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasshd", supplied_args))

        if "rasshd" in self.call_count:
            self.call_count["rasshd"] += 1
        else:
            self.call_count["rasshd"] = 1

        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def dis2mph(self, cpx1: str, cpx2: str, width1, width2, start = None, nlines = None, roff = None, azoff = None, scale = None, exp = None, sc_abs1 = None, sc_abs2 = None, dtype = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2mph", supplied_args))

        if "dis2mph" in self.call_count:
            self.call_count["dis2mph"] += 1
        else:
            self.call_count["dis2mph"] = 1

        if cpx1 is not None:
            result = self._validate(Path(cpx1).exists(), result)
        if cpx2 is not None:
            result = self._validate(Path(cpx2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def dis2SLC(self, SLC1: str, SLC2: str, width1, width2, start = None, nlines = None, roff = None, azoff = None, scale = None, exp = None, data_type = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2SLC", supplied_args))

        if "dis2SLC" in self.call_count:
            self.call_count["dis2SLC"] += 1
        else:
            self.call_count["dis2SLC"] = 1

        if SLC1 is not None:
            result = self._validate(Path(SLC1).exists(), result)
        if SLC2 is not None:
            result = self._validate(Path(SLC2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def raspwr(self, pwr: str, width, start = None, nlines = None, pixavr = None, pixavaz = None, scale = None, exp = None, LR = None, rasf: str = None, data_type = None, hdrz = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "raspwr", supplied_args))

        if "raspwr" in self.call_count:
            self.call_count["raspwr"] += 1
        else:
            self.call_count["raspwr"] = 1

        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        valid_values = [0, 1, 2]
        result = self._validate(data_type in valid_values, result)
        return result

    def disSLC(self, SLC: str, width, start = None, nlines = None, scale = None, exp = None, data_type = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disSLC", supplied_args))

        if "disSLC" in self.call_count:
            self.call_count["disSLC"] += 1
        else:
            self.call_count["disSLC"] = 1

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def ras8_colormap(self, model, h0, hrange, ival, sat, cm: str, cm_ras: str = None, width = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras8_colormap", supplied_args))

        if "ras8_colormap" in self.call_count:
            self.call_count["ras8_colormap"] += 1
        else:
            self.call_count["ras8_colormap"] = 1

        valid_values = [0, 1, 2, 3]
        result = self._validate(model in valid_values, result)
        if cm is not None:
            Path(cm).touch()
        if cm_ras is not None:
            Path(cm_ras).touch()
        return result

    def cp_data(self, infile: str, outfile: str, lbytes, start = None, nlines = None, offset = None, file_ldr = None, offb = None, nbyte = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "cp_data", supplied_args))

        if "cp_data" in self.call_count:
            self.call_count["cp_data"] += 1
        else:
            self.call_count["cp_data"] = 1

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
            Path(outfile).touch()
        return result

    def ras_cpt(self, data: str, width, cpt: str, color_model = None, start = None, nlines = None, pixavr = None, pixavaz = None, LR = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras_cpt", supplied_args))

        if "ras_cpt" in self.call_count:
            self.call_count["ras_cpt"] += 1
        else:
            self.call_count["ras_cpt"] = 1

        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if cpt is not None:
            result = self._validate(Path(cpt).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(color_model in valid_values, result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def uchar2float(self, infile, outfile, scale, exp, offset = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "uchar2float", supplied_args))

        if "uchar2float" in self.call_count:
            self.call_count["uchar2float"] += 1
        else:
            self.call_count["uchar2float"] = 1

        return result

    def rasdt_pwr24(self, data: str, pwr: str, width, start_data = None, start_pwr = None, nlines = None, pixavr = None, pixavaz = None, cycle = None, scale = None, exp = None, LR = None, rasf: str = None, cc = None, start_cc = None, cc_min = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasdt_pwr24", supplied_args))

        if "rasdt_pwr24" in self.call_count:
            self.call_count["rasdt_pwr24"] += 1
        else:
            self.call_count["rasdt_pwr24"] = 1

        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
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

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
            Path(outfile).touch()
        return result

    def disbyte(self, image: str, width, start = None, nlines = None, scale = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disbyte", supplied_args))

        if "disbyte" in self.call_count:
            self.call_count["disbyte"] += 1
        else:
            self.call_count["disbyte"] = 1

        if image is not None:
            result = self._validate(Path(image).exists(), result)
        return result

    def rasrmg(self, unw: str, pwr: str, width, start_unw = None, start_pwr = None, nlines = None, pixavr = None, pixavaz = None, ph_scale = None, scale = None, exp = None, ph_offset = None, LR = None, rasf: str = None, cc: str = None, start_cc = None, cc_min = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasrmg", supplied_args))

        if "rasrmg" in self.call_count:
            self.call_count["rasrmg"] += 1
        else:
            self.call_count["rasrmg"] = 1

        if unw is not None:
            result = self._validate(Path(unw).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        if cc is not None:
            result = self._validate(Path(cc).exists(), result)
        return result

    def kml_plan(self, MLI_par: str, DEM_par: str, lookup_table: str, kml: str, geoid: str = None, geoid_par: str = None, extension = None, flight_path = None, t_event = None, pt_list: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "kml_plan", supplied_args))

        if "kml_plan" in self.call_count:
            self.call_count["kml_plan"] += 1
        else:
            self.call_count["kml_plan"] = 1

        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if lookup_table is not None:
            result = self._validate(Path(lookup_table).exists(), result)
        if kml is not None:
            Path(kml).touch()
        if geoid is not None:
            result = self._validate(Path(geoid).exists(), result)
        if geoid_par is not None:
            result = self._validate(Path(geoid_par).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(flight_path in valid_values, result)
        if pt_list is not None:
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

        if image is not None:
            result = self._validate(Path(image).exists(), result)
        if dem_par is not None:
            result = self._validate(Path(dem_par).exists(), result)
        if kml is not None:
            Path(kml).touch()
        return result

    def dismph_pk(self, cpx: str, width, start = None, nlines = None, scale = None, exp = None, data_type = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph_pk", supplied_args))

        if "dismph_pk" in self.call_count:
            self.call_count["dismph_pk"] += 1
        else:
            self.call_count["dismph_pk"] = 1

        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def ras_dB(self, pwr: str, width, start = None, nlines = None, pixavr = None, pixavaz = None, min_dB = None, max_dB = None, dB_offset = None, LR = None, rasf: str = None, abs_flag = None, inverse = None, channel = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras_dB", supplied_args))

        if "ras_dB" in self.call_count:
            self.call_count["ras_dB"] += 1
        else:
            self.call_count["ras_dB"] = 1

        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        valid_values = [0, 1]
        result = self._validate(abs_flag in valid_values, result)
        valid_values = [1, 2, 3]
        result = self._validate(channel in valid_values, result)
        return result

    def thres_data(self, data_in: str, width, data_out: str, t_data: str, t_min, t_max, dtype = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "thres_data", supplied_args))

        if "thres_data" in self.call_count:
            self.call_count["thres_data"] += 1
        else:
            self.call_count["thres_data"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        if t_data is not None:
            result = self._validate(Path(t_data).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(dtype in valid_values, result)
        return result

    def distree(self, flag: str, unw: str, cpx: str, width, start = None, nlines = None, ph_scale = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "distree", supplied_args))

        if "distree" in self.call_count:
            self.call_count["distree"] += 1
        else:
            self.call_count["distree"] = 1

        if flag is not None:
            result = self._validate(Path(flag).exists(), result)
        if unw is not None:
            result = self._validate(Path(unw).exists(), result)
        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        return result

    def dis_linear(self, pwr: str, width, start = None, nlines = None, min = None, max = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis_linear", supplied_args))

        if "dis_linear" in self.call_count:
            self.call_count["dis_linear"] += 1
        else:
            self.call_count["dis_linear"] = 1

        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        return result

    def disdem_par(self, DEM: str, DEM_par: str, start = None, nlines = None, exaggerate = None, theta0 = None, phi0 = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disdem_par", supplied_args))

        if "disdem_par" in self.call_count:
            self.call_count["disdem_par"] += 1
        else:
            self.call_count["disdem_par"] = 1

        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        return result

    def gcp_2ras(self, ras1: str, ras2: str, gcp: str, mag = None, win_sz = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "gcp_2ras", supplied_args))

        if "gcp_2ras" in self.call_count:
            self.call_count["gcp_2ras"] += 1
        else:
            self.call_count["gcp_2ras"] = 1

        if ras1 is not None:
            result = self._validate(Path(ras1).exists(), result)
        if ras2 is not None:
            result = self._validate(Path(ras2).exists(), result)
        if gcp is not None:
            Path(gcp).touch()
        return result

    def rashgt(self, hgt: str, pwr: str, width, start_hgt = None, start_pwr = None, nlines = None, pixavr = None, pixavaz = None, m_cycle = None, scale = None, exp = None, LR = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rashgt", supplied_args))

        if "rashgt" in self.call_count:
            self.call_count["rashgt"] += 1
        else:
            self.call_count["rashgt"] = 1

        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def ras_linear(self, pwr: str, width, start = None, nlines = None, pixavr = None, pixavaz = None, min = None, max = None, LR = None, rasf: str = None, inverse = None, channel = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras_linear", supplied_args))

        if "ras_linear" in self.call_count:
            self.call_count["ras_linear"] += 1
        else:
            self.call_count["ras_linear"] = 1

        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        valid_values = [1, 2, 3]
        result = self._validate(channel in valid_values, result)
        return result

    def tree_edit(self, flag: str, ras: str, mag = None, win_sz = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "tree_edit", supplied_args))

        if "tree_edit" in self.call_count:
            self.call_count["tree_edit"] += 1
        else:
            self.call_count["tree_edit"] = 1

        if flag is not None:
            result = self._validate(Path(flag).exists(), result)
        if ras is not None:
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

        if d1 is not None:
            result = self._validate(Path(d1).exists(), result)
        if d2 is not None:
            result = self._validate(Path(d2).exists(), result)
        if dout is not None:
            Path(dout).touch()
        return result

    def dis2pwr(self, pwr1: str, pwr2: str, width1, width2, start = None, nlines = None, roff = None, azoff = None, scale = None, exp = None, dtype = None, sc_abs1 = None, sc_abs2 = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2pwr", supplied_args))

        if "dis2pwr" in self.call_count:
            self.call_count["dis2pwr"] += 1
        else:
            self.call_count["dis2pwr"] = 1

        if pwr1 is not None:
            result = self._validate(Path(pwr1).exists(), result)
        if pwr2 is not None:
            result = self._validate(Path(pwr2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def dispwr(self, pwr: str, width, start = None, nlines = None, scale = None, exp = None, data_type = None, sc_abs = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dispwr", supplied_args))

        if "dispwr" in self.call_count:
            self.call_count["dispwr"] += 1
        else:
            self.call_count["dispwr"] = 1

        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def rasSLC(self, SLC: str, width, start = None, nlines = None, pixavr = None, pixavaz = None, scale = None, exp = None, LR = None, data_type = None, header = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasSLC", supplied_args))

        if "rasSLC" in self.call_count:
            self.call_count["rasSLC"] += 1
        else:
            self.call_count["rasSLC"] = 1

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def rasmph_pwr24(self, cpx: str, pwr: str, width, start_cpx = None, start_pwr = None, nlines = None, pixavr = None, pixavaz = None, scale = None, exp = None, LR = None, rasf: str = None, cc = None, start_cc = None, cc_min = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasmph_pwr24", supplied_args))

        if "rasmph_pwr24" in self.call_count:
            self.call_count["rasmph_pwr24"] += 1
        else:
            self.call_count["rasmph_pwr24"] = 1

        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def dismph_ub(self, cpx: str, width, start = None, nlines = None, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph_ub", supplied_args))

        if "dismph_ub" in self.call_count:
            self.call_count["dismph_ub"] += 1
        else:
            self.call_count["dismph_ub"] = 1

        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        return result

    def short2float(self, infile: str, outfile: str, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "short2float", supplied_args))

        if "short2float" in self.call_count:
            self.call_count["short2float"] += 1
        else:
            self.call_count["short2float"] = 1

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
            Path(outfile).touch()
        return result

    def ras8_color_scale(self, rasf: str, color_model = None, h0 = None, hrange = None, ival = None, sat = None, chip_width = None, gap = None, chip_height = None, nval = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras8_color_scale", supplied_args))

        if "ras8_color_scale" in self.call_count:
            self.call_count["ras8_color_scale"] += 1
        else:
            self.call_count["ras8_color_scale"] = 1

        if rasf is not None:
            Path(rasf).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(color_model in valid_values, result)
        return result

    def rasmph(self, cpx: str, width, start = None, nlines = None, pixavr = None, pixavaz = None, scale = None, exp = None, LR = None, rasf: str = None, data_type = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasmph", supplied_args))

        if "rasmph" in self.call_count:
            self.call_count["rasmph"] += 1
        else:
            self.call_count["rasmph"] = 1

        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        if rasf is not None:
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

    def polyras(self, ras: str, mag = None, win_sz = None, poly_file: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "polyras", supplied_args))

        if "polyras" in self.call_count:
            self.call_count["polyras"] += 1
        else:
            self.call_count["polyras"] = 1

        if ras is not None:
            result = self._validate(Path(ras).exists(), result)
        if poly_file is not None:
            Path(poly_file).touch()
        return result

    def rastree(self, flag: str, unw: str, cpx, width, start = None, nlines = None, ph_scale = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rastree", supplied_args))

        if "rastree" in self.call_count:
            self.call_count["rastree"] += 1
        else:
            self.call_count["rastree"] = 1

        if flag is not None:
            result = self._validate(Path(flag).exists(), result)
        if unw is not None:
            result = self._validate(Path(unw).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def disshd(self, DEM: str, width, col_post, row_post = None, start = None, nlines = None, theta0 = None, phi0 = None, data_type = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "disshd", supplied_args))

        if "disshd" in self.call_count:
            self.call_count["disshd"] += 1
        else:
            self.call_count["disshd"] = 1

        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        valid_values = [0, 1]
        result = self._validate(data_type in valid_values, result)
        return result

    def dis2gbyte(self, image1: str, image2: str, width1, width2, start = None, nlines = None, roff = None, azoff = None, scale = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2gbyte", supplied_args))

        if "dis2gbyte" in self.call_count:
            self.call_count["dis2gbyte"] += 1
        else:
            self.call_count["dis2gbyte"] = 1

        if image1 is not None:
            result = self._validate(Path(image1).exists(), result)
        if image2 is not None:
            result = self._validate(Path(image2).exists(), result)
        return result

    def ras8_float(self, f1: str, f2: str, width, rasf: str, color_model = None, h0 = None, hrange = None, imin = None, imax = None, sat = None, sc1 = None, A1 = None, B1 = None, cyclic1 = None, sc2 = None, A2 = None, B2 = None, start_f1 = None, start_f2 = None, nlines = None, pixavr = None, pixavaz = None, LR = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras8_float", supplied_args))

        if "ras8_float" in self.call_count:
            self.call_count["ras8_float"] += 1
        else:
            self.call_count["ras8_float"] = 1

        if f1 is not None:
            result = self._validate(Path(f1).exists(), result)
        if f2 is not None:
            result = self._validate(Path(f2).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(color_model in valid_values, result)
        return result

    def flip(self, infile: str, outfile: str, width, format = None, sense = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "flip", supplied_args))

        if "flip" in self.call_count:
            self.call_count["flip"] += 1
        else:
            self.call_count["flip"] = 1

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
            Path(outfile).touch()
        return result

    def create_array(self, output: str, width, nlines, dtype = None, val = None, val_im = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "create_array", supplied_args))

        if "create_array" in self.call_count:
            self.call_count["create_array"] += 1
        else:
            self.call_count["create_array"] = 1

        if output is not None:
            Path(output).touch()
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(dtype in valid_values, result)
        return result

    def ascii2float(self, data_in: str, width, data_out: str, loff = None, nl = None, coff = None, nv = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ascii2float", supplied_args))

        if "ascii2float" in self.call_count:
            self.call_count["ascii2float"] += 1
        else:
            self.call_count["ascii2float"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        return result

    def data2tiff(self, data: str, width, type, TIFF: str, nodata = None, xspacing = None, yspacing = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "data2tiff", supplied_args))

        if "data2tiff" in self.call_count:
            self.call_count["data2tiff"] += 1
        else:
            self.call_count["data2tiff"] = 1

        if data is not None:
            result = self._validate(Path(data).exists(), result)
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(type in valid_values, result)
        if TIFF is not None:
            Path(TIFF).touch()
        return result

    def mapshd(self, DEM: str, width, col_post, row_post, theta0, phi0, shade: str, dtype = None, zero_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "mapshd", supplied_args))

        if "mapshd" in self.call_count:
            self.call_count["mapshd"] += 1
        else:
            self.call_count["mapshd"] = 1

        if DEM is not None:
            result = self._validate(Path(DEM).exists(), result)
        if shade is not None:
            Path(shade).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def kml_pt(self, table: str, lat_col, lon_col, val1_col, val1_label, val2_col, val2_label, val3_col, val3_label, id_col, kml: str, icon_URL = None, logo_URL = None, legend_URL = None, color_model = None, h0 = None, hrange = None, imin = None, imax = None, sat_min = None, sat_max = None, sc1 = None, A1 = None, B1 = None, cyclic1 = None, sc2 = None, A2 = None, B2 = None, start_f1 = None, start_f2 = None, B3 = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "kml_pt", supplied_args))

        if "kml_pt" in self.call_count:
            self.call_count["kml_pt"] += 1
        else:
            self.call_count["kml_pt"] = 1

        if table is not None:
            result = self._validate(Path(table).exists(), result)
        if kml is not None:
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

        if vec is not None:
            result = self._validate(Path(vec).exists(), result)
        valid_values = [1, 2, 3]
        result = self._validate(index in valid_values, result)
        if cmp is not None:
            Path(cmp).touch()
        return result

    def dismph_pwr24(self, cpx: str, pwr: str, width, start_cpx = None, start_pwr = None, nlines = None, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dismph_pwr24", supplied_args))

        if "dismph_pwr24" in self.call_count:
            self.call_count["dismph_pwr24"] += 1
        else:
            self.call_count["dismph_pwr24"] = 1

        if cpx is not None:
            result = self._validate(Path(cpx).exists(), result)
        if pwr is not None:
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

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
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

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
            Path(outfile).touch()
        return result

    def dis2byte(self, image1: str, image2: str, width1, width2, start = None, nlines = None, roff = None, azoff = None, scale = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2byte", supplied_args))

        if "dis2byte" in self.call_count:
            self.call_count["dis2byte"] += 1
        else:
            self.call_count["dis2byte"] = 1

        if image1 is not None:
            result = self._validate(Path(image1).exists(), result)
        if image2 is not None:
            result = self._validate(Path(image2).exists(), result)
        return result

    def svg_poly(self, image: str, dem_par: str, poly: str, svg: str, width = None, nlines = None, thick = None, lcolor = None, lopac = None, pcolor = None, popac = None, tcolor = None, font = None, fsize = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "svg_poly", supplied_args))

        if "svg_poly" in self.call_count:
            self.call_count["svg_poly"] += 1
        else:
            self.call_count["svg_poly"] = 1

        if image is not None:
            result = self._validate(Path(image).exists(), result)
        if dem_par is not None:
            result = self._validate(Path(dem_par).exists(), result)
        if poly is not None:
            result = self._validate(Path(poly).exists(), result)
        if svg is not None:
            Path(svg).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        result = self._validate(thick in valid_values, result)
        return result

    def rasbyte(self, raw: str, width, start = None, nlines = None, pixavr = None, pixavaz = None, scale = None, LR = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "rasbyte", supplied_args))

        if "rasbyte" in self.call_count:
            self.call_count["rasbyte"] += 1
        else:
            self.call_count["rasbyte"] = 1

        if raw is not None:
            result = self._validate(Path(raw).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def ras_cpt_scale(self, rasf: str, cpt: str, color_model = None, width = None, nlines = None, start_value = None, end_value = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "ras_cpt_scale", supplied_args))

        if "ras_cpt_scale" in self.call_count:
            self.call_count["ras_cpt_scale"] += 1
        else:
            self.call_count["ras_cpt_scale"] = 1

        if rasf is not None:
            Path(rasf).touch()
        if cpt is not None:
            result = self._validate(Path(cpt).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(color_model in valid_values, result)
        return result

    def dis2cc(self, cc1: str, cc2: str, width1, width2, start = None, nlines = None, roff = None, azoff = None, cmin = None, cmax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "dis2cc", supplied_args))

        if "dis2cc" in self.call_count:
            self.call_count["dis2cc"] += 1
        else:
            self.call_count["dis2cc"] = 1

        if cc1 is not None:
            result = self._validate(Path(cc1).exists(), result)
        if cc2 is not None:
            result = self._validate(Path(cc2).exists(), result)
        return result

    def float2ascii(self, din: str, width, data_out: str, loff = None, nl = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "float2ascii", supplied_args))

        if "float2ascii" in self.call_count:
            self.call_count["float2ascii"] += 1
        else:
            self.call_count["float2ascii"] = 1

        if din is not None:
            result = self._validate(Path(din).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        return result

    def replace_values(self, f_in: str, value, new_value, f_out: str, width, rpl_flg = None, dtype = None, zflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("DISP", "replace_values", supplied_args))

        if "replace_values" in self.call_count:
            self.call_count["replace_values"] += 1
        else:
            self.call_count["replace_values"] = 1

        if f_in is not None:
            result = self._validate(Path(f_in).exists(), result)
        if f_out is not None:
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

        if sbi_unw is not None:
            result = self._validate(Path(sbi_unw).exists(), result)
        if SLCf_par is not None:
            result = self._validate(Path(SLCf_par).exists(), result)
        if SLCb_par is not None:
            result = self._validate(Path(SLCb_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if az_offset is not None:
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

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if GRD_par is not None:
            Path(GRD_par).touch()
        if GRD is not None:
            Path(GRD).touch()
        return result

    def par_ASF_SLC(self, CEOS_SAR_leader, SLC_par: str, CEOS_data: str = None, SLC: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ASF_SLC", supplied_args))

        if "par_ASF_SLC" in self.call_count:
            self.call_count["par_ASF_SLC"] += 1
        else:
            self.call_count["par_ASF_SLC"] = 1

        if SLC_par is not None:
            Path(SLC_par).touch()
        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if SLC is not None:
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

        if HDF5 is not None:
            result = self._validate(Path(HDF5).exists(), result)
        if trunk is not None:
            Path(trunk).touch()
        return result

    def ScanSAR_full_aperture_SLC(self, SLC1_tab: str, SLC2_tab: str, SLCR_tab: str = None, SLCR_dir = None, vmode = None, wflg = None, imode = None, order = None, n_ovr = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_full_aperture_SLC", supplied_args))

        if "ScanSAR_full_aperture_SLC" in self.call_count:
            self.call_count["ScanSAR_full_aperture_SLC"] += 1
        else:
            self.call_count["ScanSAR_full_aperture_SLC"] = 1

        if SLC1_tab is not None:
            result = self._validate(Path(SLC1_tab).exists(), result)
        if SLC2_tab is not None and not Path(SLC2_tab).exists():
            Path(SLC2_tab).touch()
        if SLCR_tab is not None:
            result = self._validate(Path(SLCR_tab).exists(), result)
        valid_values = [0, 1]
        result = self._validate(vmode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(wflg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(imode in valid_values, result)
        return result

    def init_offset(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, rlks = None, azlks = None, rpos = None, azpos = None, offr = None, offaz = None, thres = None, rwin = None, azwin = None, cflag = None, deramp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "init_offset", supplied_args))

        if "init_offset" in self.call_count:
            self.call_count["init_offset"] += 1
        else:
            self.call_count["init_offset"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        valid_values = [0, 1]
        result = self._validate(cflag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(deramp in valid_values, result)
        return result

    def bridge(self, int: str, flag: str, unw: str, bridge: str, width, xmin = None, xmax = None, ymin = None, ymax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "bridge", supplied_args))

        if "bridge" in self.call_count:
            self.call_count["bridge"] += 1
        else:
            self.call_count["bridge"] = 1

        if int is not None:
            result = self._validate(Path(int).exists(), result)
        if flag is not None:
            result = self._validate(Path(flag).exists(), result)
        if unw is not None and not Path(unw).exists():
            Path(unw).touch()
        if bridge is not None:
            result = self._validate(Path(bridge).exists(), result)
        return result

    def par_ERSDAC_PALSAR(self, VEXCEL_SLC_par, SLC_par: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ERSDAC_PALSAR", supplied_args))

        if "par_ERSDAC_PALSAR" in self.call_count:
            self.call_count["par_ERSDAC_PALSAR"] += 1
        else:
            self.call_count["par_ERSDAC_PALSAR"] = 1

        if SLC_par is not None:
            Path(SLC_par).touch()
        return result

    def SR_to_GRD(self, MLI_par: str, OFF_par: str, GRD_par: str, in_file: str, out_file: str, rlks = None, azlks = None, interp_mode = None, grd_rsp = None, grd_azsp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SR_to_GRD", supplied_args))

        if "SR_to_GRD" in self.call_count:
            self.call_count["SR_to_GRD"] += 1
        else:
            self.call_count["SR_to_GRD"] = 1

        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if GRD_par is not None and not Path(GRD_par).exists():
            Path(GRD_par).touch()
        if in_file is not None:
            result = self._validate(Path(in_file).exists(), result)
        if out_file is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        return result

    def offset_fit(self, offs: str, ccp: str, OFF_par: str, coffs: str = None, coffsets: str = None, thres = None, npoly = None, interact_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_fit", supplied_args))

        if "offset_fit" in self.call_count:
            self.call_count["offset_fit"] += 1
        else:
            self.call_count["offset_fit"] = 1

        if offs is not None:
            result = self._validate(Path(offs).exists(), result)
        if ccp is not None:
            result = self._validate(Path(ccp).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if coffs is not None:
            Path(coffs).touch()
        if coffsets is not None:
            Path(coffsets).touch()
        valid_values = [0, 1]
        result = self._validate(interact_flag in valid_values, result)
        return result

    def radcal_PRI(self, PRI: str, PRI_PAR: str, GRD: str, GRD_PAR: str, K_dB = None, inc_ref = None, roff = None, nr = None, loff = None, nl = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "radcal_PRI", supplied_args))

        if "radcal_PRI" in self.call_count:
            self.call_count["radcal_PRI"] += 1
        else:
            self.call_count["radcal_PRI"] = 1

        if PRI is not None:
            result = self._validate(Path(PRI).exists(), result)
        if PRI_PAR is not None:
            result = self._validate(Path(PRI_PAR).exists(), result)
        if GRD is not None:
            Path(GRD).touch()
        if GRD_PAR is not None:
            Path(GRD_PAR).touch()
        return result

    def gcp_phase(self, unw: str, OFF_par: str, gcp: str, gcp_ph: str, win_sz = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "gcp_phase", supplied_args))

        if "gcp_phase" in self.call_count:
            self.call_count["gcp_phase"] += 1
        else:
            self.call_count["gcp_phase"] = 1

        if unw is not None:
            result = self._validate(Path(unw).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if gcp is not None:
            result = self._validate(Path(gcp).exists(), result)
        if gcp_ph is not None:
            Path(gcp_ph).touch()
        return result

    def par_MSP(self, SAR_par: str, PROC_par: str, SLC_MLI_par: str = None, image_format = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_MSP", supplied_args))

        if "par_MSP" in self.call_count:
            self.call_count["par_MSP"] += 1
        else:
            self.call_count["par_MSP"] = 1

        if SAR_par is not None:
            result = self._validate(Path(SAR_par).exists(), result)
        if PROC_par is not None:
            result = self._validate(Path(PROC_par).exists(), result)
        if SLC_MLI_par is not None:
            Path(SLC_MLI_par).touch()
        valid_values = [0, 1, 2]
        result = self._validate(image_format in valid_values, result)
        return result

    def offset_pwr_tracking2(self, SLC1: str, SLC2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, offs: str, ccp: str, OFF_par2: str = None, offs2: str = None, rwin = None, azwin = None, offsets: str = None, n_ovr = None, thres = None, rstep = None, azstep = None, rstart = None, rstop = None, azstart = None, azstop = None, bw_frac = None, deramp = None, int_filt = None, pflag = None, pltflg = None, ccs: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_pwr_tracking2", supplied_args))

        if "offset_pwr_tracking2" in self.call_count:
            self.call_count["offset_pwr_tracking2"] += 1
        else:
            self.call_count["offset_pwr_tracking2"] = 1

        if SLC1 is not None:
            result = self._validate(Path(SLC1).exists(), result)
        if SLC2 is not None:
            result = self._validate(Path(SLC2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if offs is not None:
            Path(offs).touch()
        if ccp is not None:
            Path(ccp).touch()
        if OFF_par2 is not None:
            result = self._validate(Path(OFF_par2).exists(), result)
        if offs2 is not None:
            result = self._validate(Path(offs2).exists(), result)
        if offsets is not None:
            Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(deramp in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(int_filt in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        if ccs is not None:
            Path(ccs).touch()
        return result

    def offset_tracking(self, offs: str, ccp: str, SLC_par: str, OFF_par: str, disp_map: str, disp_val: str = None, mode = None, thres = None, poly_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_tracking", supplied_args))

        if "offset_tracking" in self.call_count:
            self.call_count["offset_tracking"] += 1
        else:
            self.call_count["offset_tracking"] = 1

        if offs is not None:
            result = self._validate(Path(offs).exists(), result)
        if ccp is not None:
            result = self._validate(Path(ccp).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if disp_map is not None:
            Path(disp_map).touch()
        if disp_val is not None:
            Path(disp_val).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(poly_flag in valid_values, result)
        return result

    def SLC_interp_ScanSAR(self, SLC2_tab: str, SLC2_par: str, SLC1_tab: str, SLC1_par: str, OFF_par: str, SLC2R_tab: str, SLC_2R: str = None, SLC2R_par: str = None, mode = None, order = None, SLC2R_dir = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_interp_ScanSAR", supplied_args))

        if "SLC_interp_ScanSAR" in self.call_count:
            self.call_count["SLC_interp_ScanSAR"] += 1
        else:
            self.call_count["SLC_interp_ScanSAR"] = 1

        if SLC2_tab is not None:
            result = self._validate(Path(SLC2_tab).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if SLC1_tab is not None:
            result = self._validate(Path(SLC1_tab).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if SLC2R_tab is not None and not Path(SLC2R_tab).exists():
            Path(SLC2R_tab).touch()
        if SLC_2R is not None:
            Path(SLC_2R).touch()
        if SLC2R_par is not None:
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

        if d1 is not None:
            result = self._validate(Path(d1).exists(), result)
        if d2 is not None:
            result = self._validate(Path(d2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def res_map(self, hgt: str, gr: str, data: str, SLC_par: str, OFF_par: str, res_hgt: str, res_data: str, nr = None, naz = None, azps_res = None, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "res_map", supplied_args))

        if "res_map" in self.call_count:
            self.call_count["res_map"] += 1
        else:
            self.call_count["res_map"] = 1

        if hgt is not None:
            result = self._validate(Path(hgt).exists(), result)
        if gr is not None:
            result = self._validate(Path(gr).exists(), result)
        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if res_hgt is not None:
            Path(res_hgt).touch()
        if res_data is not None:
            Path(res_data).touch()
        return result

    def RSAT2_vec(self, SLC_par: str, RSAT2_orb, nstate = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "RSAT2_vec", supplied_args))

        if "RSAT2_vec" in self.call_count:
            self.call_count["RSAT2_vec"] += 1
        else:
            self.call_count["RSAT2_vec"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        return result

    def par_S1_GRD(self, GeoTIFF: str, annotation_XML: str, calibration_XML: str, noise_XML: str, MLI_par: str, MLI: str, GRD_par: str = None, GRD: str = None, eflg = None, rps = None, noise_pwr = None, edge_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_S1_GRD", supplied_args))

        if "par_S1_GRD" in self.call_count:
            self.call_count["par_S1_GRD"] += 1
        else:
            self.call_count["par_S1_GRD"] = 1

        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if annotation_XML is not None:
            result = self._validate(Path(annotation_XML).exists(), result)
        if calibration_XML is not None:
            result = self._validate(Path(calibration_XML).exists(), result)
        if noise_XML is not None:
            result = self._validate(Path(noise_XML).exists(), result)
        if MLI_par is not None:
            Path(MLI_par).touch()
        if MLI is not None:
            Path(MLI).touch()
        if GRD_par is not None:
            Path(GRD_par).touch()
        if GRD is not None:
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

        if SLC1_tab is not None:
            result = self._validate(Path(SLC1_tab).exists(), result)
        if SLC2_tab is not None:
            result = self._validate(Path(SLC2_tab).exists(), result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(phflg in valid_values, result)
        return result

    def par_RSAT_SGF(self, CEOS_leader: str, CEOS_data: str, GRD_par: str, GRD: str, sc_dB = None, dt = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RSAT_SGF", supplied_args))

        if "par_RSAT_SGF" in self.call_count:
            self.call_count["par_RSAT_SGF"] += 1
        else:
            self.call_count["par_RSAT_SGF"] = 1

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if GRD_par is not None:
            Path(GRD_par).touch()
        if GRD is not None:
            Path(GRD).touch()
        return result

    def par_ICEYE_GRD(self, GeoTIFF: str, MLI_par: str, mli = None, GRD_par: str = None, GRD = None, rps = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ICEYE_GRD", supplied_args))

        if "par_ICEYE_GRD" in self.call_count:
            self.call_count["par_ICEYE_GRD"] += 1
        else:
            self.call_count["par_ICEYE_GRD"] = 1

        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if GRD_par is not None:
            Path(GRD_par).touch()
        return result

    def par_CS_SLC_TIF(self, GeoTIFF: str, XML: str, trunk: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_CS_SLC_TIF", supplied_args))

        if "par_CS_SLC_TIF" in self.call_count:
            self.call_count["par_CS_SLC_TIF"] += 1
        else:
            self.call_count["par_CS_SLC_TIF"] = 1

        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if XML is not None:
            result = self._validate(Path(XML).exists(), result)
        if trunk is not None:
            Path(trunk).touch()
        return result

    def DORIS_vec(self, SLC_PAR: str, DOR: str, nstate = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "DORIS_vec", supplied_args))

        if "DORIS_vec" in self.call_count:
            self.call_count["DORIS_vec"] += 1
        else:
            self.call_count["DORIS_vec"] = 1

        if SLC_PAR is not None and not Path(SLC_PAR).exists():
            Path(SLC_PAR).touch()
        if DOR is not None:
            result = self._validate(Path(DOR).exists(), result)
        return result

    def bpf(self, data_in: str, data_out: str, width, fc_x, bw_x, fc_y, bw_y, roff = None, azoff = None, nr = None, naz = None, data_type = None, f_mode = None, beta = None, fir_len = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "bpf", supplied_args))

        if "bpf" in self.call_count:
            self.call_count["bpf"] += 1
        else:
            self.call_count["bpf"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
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

        if image is not None:
            result = self._validate(Path(image).exists(), result)
        return result

    def clear_flag(self, flag: str, width, flag_bits, xmin, xmax, ymin = None, ymax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "clear_flag", supplied_args))

        if "clear_flag" in self.call_count:
            self.call_count["clear_flag"] += 1
        else:
            self.call_count["clear_flag"] = 1

        if flag is not None:
            result = self._validate(Path(flag).exists(), result)
        return result

    def par_NovaSAR_GRD(self, GeoTIFF: str, XML: str, polarization, MLI_par: str, MLI: str = None, GRD_par: str = None, GRD: str = None, rps = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_NovaSAR_GRD", supplied_args))

        if "par_NovaSAR_GRD" in self.call_count:
            self.call_count["par_NovaSAR_GRD"] += 1
        else:
            self.call_count["par_NovaSAR_GRD"] = 1

        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if XML is not None:
            result = self._validate(Path(XML).exists(), result)
        if MLI_par is not None:
            Path(MLI_par).touch()
        if MLI is not None:
            Path(MLI).touch()
        if GRD_par is not None:
            Path(GRD_par).touch()
        if GRD is not None:
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

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if slr is not None:
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

        if aux_data is not None:
            result = self._validate(Path(aux_data).exists(), result)
        if slc_Re is not None:
            result = self._validate(Path(slc_Re).exists(), result)
        if slc_Im is not None:
            result = self._validate(Path(slc_Im).exists(), result)
        if date is not None:
            result = self._validate(Path(date).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        return result

    def corr_flag(self, corr: str, flag: str, width, corr_thr, xmin = None, xmax = None, ymin = None, ymax = None, border = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "corr_flag", supplied_args))

        if "corr_flag" in self.call_count:
            self.call_count["corr_flag"] += 1
        else:
            self.call_count["corr_flag"] = 1

        if corr is not None:
            result = self._validate(Path(corr).exists(), result)
        if flag is not None and not Path(flag).exists():
            Path(flag).touch()
        return result

    def ptarg_cal_MLI(self, MLI_par: str, MLI: str, r_samp, az_samp, psigma, c_r_samp, c_az_samp, ptr_image: str, r_plot: str, az_plot: str, pcal: str, osf = None, win = None, pltflg = None, psz = None, csz = None, theta_inc = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ptarg_cal_MLI", supplied_args))

        if "ptarg_cal_MLI" in self.call_count:
            self.call_count["ptarg_cal_MLI"] += 1
        else:
            self.call_count["ptarg_cal_MLI"] = 1

        if MLI_par is not None:
            result = self._validate(Path(MLI_par).exists(), result)
        if MLI is not None:
            result = self._validate(Path(MLI).exists(), result)
        if ptr_image is not None:
            Path(ptr_image).touch()
        if r_plot is not None:
            Path(r_plot).touch()
        if az_plot is not None:
            Path(az_plot).touch()
        if pcal is not None:
            Path(pcal).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        return result

    def ORB_filt(self, SLC_par_in: str, SLC_par_out: str, interval = None, extra = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ORB_filt", supplied_args))

        if "ORB_filt" in self.call_count:
            self.call_count["ORB_filt"] += 1
        else:
            self.call_count["ORB_filt"] = 1

        if SLC_par_in is not None:
            result = self._validate(Path(SLC_par_in).exists(), result)
        if SLC_par_out is not None:
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

        if offs is not None:
            result = self._validate(Path(offs).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if offs_sub is not None:
            Path(offs_sub).touch()
        return result

    def radcal_MLI(self, MLI: str, MLI_PAR: str, OFF_par: str, CMLI: str, antenna: str = None, rloss_flag = None, ant_flag = None, refarea_flag = None, sc_dB = None, K_dB = None, pix_area: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "radcal_MLI", supplied_args))

        if "radcal_MLI" in self.call_count:
            self.call_count["radcal_MLI"] += 1
        else:
            self.call_count["radcal_MLI"] = 1

        if MLI is not None:
            result = self._validate(Path(MLI).exists(), result)
        if MLI_PAR is not None:
            result = self._validate(Path(MLI_PAR).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if CMLI is not None:
            Path(CMLI).touch()
        if antenna is not None:
            result = self._validate(Path(antenna).exists(), result)
        if pix_area is not None:
            Path(pix_area).touch()
        return result

    def S1_OPOD_vec(self, SLC_PAR: str, OPOD: str, nstate = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "S1_OPOD_vec", supplied_args))

        if "S1_OPOD_vec" in self.call_count:
            self.call_count["S1_OPOD_vec"] += 1
        else:
            self.call_count["S1_OPOD_vec"] = 1

        if SLC_PAR is not None and not Path(SLC_PAR).exists():
            Path(SLC_PAR).touch()
        if OPOD is not None:
            result = self._validate(Path(OPOD).exists(), result)
        return result

    def par_RCM_GRD(self, RCM_dir: str, polarization, radcal, noise, MLI_par: str = None, MLI: str = None, GRD_par: str = None, GRD: str = None, rps = None, noise_pwr: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RCM_GRD", supplied_args))

        if "par_RCM_GRD" in self.call_count:
            self.call_count["par_RCM_GRD"] += 1
        else:
            self.call_count["par_RCM_GRD"] = 1

        if RCM_dir is not None:
            result = self._validate(Path(RCM_dir).exists(), result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(radcal in valid_values, result)
        if MLI_par is not None:
            Path(MLI_par).touch()
        if MLI is not None:
            Path(MLI).touch()
        if GRD_par is not None:
            Path(GRD_par).touch()
        if GRD is not None:
            Path(GRD).touch()
        if noise_pwr is not None:
            Path(noise_pwr).touch()
        return result

    def offset_pwr_tracking(self, SLC1: str, SLC2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, offs: str, ccp: str, rwin = None, azwin = None, offsets: str = None, n_ovr = None, thres = None, rstep = None, azstep = None, rstart = None, rstop = None, azstart = None, azstop = None, lanczos = None, bw_frac = None, deramp = None, int_filt = None, pflag = None, pltflg = None, ccs: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_pwr_tracking", supplied_args))

        if "offset_pwr_tracking" in self.call_count:
            self.call_count["offset_pwr_tracking"] += 1
        else:
            self.call_count["offset_pwr_tracking"] = 1

        if SLC1 is not None:
            result = self._validate(Path(SLC1).exists(), result)
        if SLC2 is not None:
            result = self._validate(Path(SLC2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if offs is not None:
            Path(offs).touch()
        if ccp is not None:
            Path(ccp).touch()
        if offsets is not None:
            Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(deramp in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(int_filt in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        if ccs is not None:
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

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_par1 is not None:
            result = self._validate(Path(SLC_par1).exists(), result)
        if SLC_2 is not None:
            Path(SLC_2).touch()
        if SLC_par2 is not None:
            Path(SLC_par2).touch()
        return result

    def unw_model(self, interf: str, unw_model: str, unw: str, width, xinit = None, yinit = None, ref_ph = None, width_model = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "unw_model", supplied_args))

        if "unw_model" in self.call_count:
            self.call_count["unw_model"] += 1
        else:
            self.call_count["unw_model"] = 1

        if interf is not None:
            result = self._validate(Path(interf).exists(), result)
        if unw_model is not None:
            result = self._validate(Path(unw_model).exists(), result)
        if unw is not None:
            Path(unw).touch()
        return result

    def par_TX_SLC(self, annotation_XML: str, COSAR: str, SLC_par: str, SLC: str, pol = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_TX_SLC", supplied_args))

        if "par_TX_SLC" in self.call_count:
            self.call_count["par_TX_SLC"] += 1
        else:
            self.call_count["par_TX_SLC"] = 1

        if annotation_XML is not None:
            result = self._validate(Path(annotation_XML).exists(), result)
        if COSAR is not None:
            result = self._validate(Path(COSAR).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        return result

    def ScanSAR_burst_MLI(self, SLC_tab: str, MLI_tab: str, rlks, azlks, bflg = None, SLCR_tab: str = None, MLI_dir = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_burst_MLI", supplied_args))

        if "ScanSAR_burst_MLI" in self.call_count:
            self.call_count["ScanSAR_burst_MLI"] += 1
        else:
            self.call_count["ScanSAR_burst_MLI"] = 1

        if SLC_tab is not None:
            result = self._validate(Path(SLC_tab).exists(), result)
        if MLI_tab is not None:
            Path(MLI_tab).touch()
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        if SLCR_tab is not None:
            result = self._validate(Path(SLCR_tab).exists(), result)
        return result

    def base_copy(self, SLC1_par: str, baseline_1: str, SLC2_par: str, baseline_2: str, time_rev = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_copy", supplied_args))

        if "base_copy" in self.call_count:
            self.call_count["base_copy"] += 1
        else:
            self.call_count["base_copy"] = 1

        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if baseline_1 is not None:
            result = self._validate(Path(baseline_1).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if baseline_2 is not None:
            Path(baseline_2).touch()
        return result

    def SLC_interp_map(self, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, SLC_2R: str, SLC2R_par: str, OFF_par2: str, coffs_sm, loff = None, nlines = None, mode = None, order = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_interp_map", supplied_args))

        if "SLC_interp_map" in self.call_count:
            self.call_count["SLC_interp_map"] += 1
        else:
            self.call_count["SLC_interp_map"] = 1

        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if SLC_2R is not None:
            Path(SLC_2R).touch()
        if SLC2R_par is not None:
            Path(SLC2R_par).touch()
        if OFF_par2 is not None:
            result = self._validate(Path(OFF_par2).exists(), result)
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def base_perp(self, baseline: str, SLC1_par: str, OFF_par: str, time_rev = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_perp", supplied_args))

        if "base_perp" in self.call_count:
            self.call_count["base_perp"] += 1
        else:
            self.call_count["base_perp"] = 1

        if baseline is not None:
            result = self._validate(Path(baseline).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        return result

    def SLC_adf(self, SLC: str, ref_SLC: str, ref_SLC_par: str, SLC_filt: str, mode = None, alpha = None, nfft_r = None, nfft_az = None, r_step = None, az_step = None, mwin_r = None, mwin_az = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_adf", supplied_args))

        if "SLC_adf" in self.call_count:
            self.call_count["SLC_adf"] += 1
        else:
            self.call_count["SLC_adf"] = 1

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if ref_SLC is not None:
            result = self._validate(Path(ref_SLC).exists(), result)
        if ref_SLC_par is not None:
            result = self._validate(Path(ref_SLC_par).exists(), result)
        if SLC_filt is not None:
            Path(SLC_filt).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(mode in valid_values, result)
        return result

    def offset_SLC(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, offs: str, snr: str, rwin = None, azwin = None, offsets: str = None, n_ovr = None, nr = None, naz = None, thres = None, ISZ = None, pflag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_SLC", supplied_args))

        if "offset_SLC" in self.call_count:
            self.call_count["offset_SLC"] += 1
        else:
            self.call_count["offset_SLC"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if offs is not None:
            Path(offs).touch()
        if snr is not None:
            Path(snr).touch()
        if offsets is not None:
            Path(offsets).touch()
        return result

    def adf(self, interf: str, sm: str, cc: str, width, alpha = None, nfft = None, cc_win = None, step = None, loff = None, nlines = None, wfrac = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "adf", supplied_args))

        if "adf" in self.call_count:
            self.call_count["adf"] += 1
        else:
            self.call_count["adf"] = 1

        if interf is not None:
            result = self._validate(Path(interf).exists(), result)
        if sm is not None:
            Path(sm).touch()
        if cc is not None:
            Path(cc).touch()
        return result

    def PRC_vec(self, SLC_par: str, PRC: str, nstate = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "PRC_vec", supplied_args))

        if "PRC_vec" in self.call_count:
            self.call_count["PRC_vec"] += 1
        else:
            self.call_count["PRC_vec"] = 1

        if SLC_par is not None and not Path(SLC_par).exists():
            Path(SLC_par).touch()
        if PRC is not None:
            result = self._validate(Path(PRC).exists(), result)
        return result

    def adapt_filt(self, int: str, sm: str, width, low_SNR_thr = None, filt_width = None, xmin = None, xmax = None, ymin = None, ymax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "adapt_filt", supplied_args))

        if "adapt_filt" in self.call_count:
            self.call_count["adapt_filt"] += 1
        else:
            self.call_count["adapt_filt"] = 1

        if int is not None:
            result = self._validate(Path(int).exists(), result)
        if sm is not None:
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

        if product_XML is not None:
            result = self._validate(Path(product_XML).exists(), result)
        if lut_XML is not None:
            result = self._validate(Path(lut_XML).exists(), result)
        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if polarization is not None:
            result = self._validate(Path(polarization).exists(), result)
        if GRD_par is not None:
            Path(GRD_par).touch()
        if GRD is not None:
            Path(GRD).touch()
        return result

    def ScanSAR_burst_corners(self, SLC_par: str, TOPS_par: str, KML: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_burst_corners", supplied_args))

        if "ScanSAR_burst_corners" in self.call_count:
            self.call_count["ScanSAR_burst_corners"] += 1
        else:
            self.call_count["ScanSAR_burst_corners"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if TOPS_par is not None:
            result = self._validate(Path(TOPS_par).exists(), result)
        if KML is not None:
            Path(KML).touch()
        return result

    def par_RCM_SLC_ScanSAR(self, RCM_dir: str, polarization, radcal, noise_in, root_name: str, SLC_tab: str = None, beam = None, noise_out = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RCM_SLC_ScanSAR", supplied_args))

        if "par_RCM_SLC_ScanSAR" in self.call_count:
            self.call_count["par_RCM_SLC_ScanSAR"] += 1
        else:
            self.call_count["par_RCM_SLC_ScanSAR"] = 1

        if RCM_dir is not None:
            result = self._validate(Path(RCM_dir).exists(), result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(radcal in valid_values, result)
        if root_name is not None:
            Path(root_name).touch()
        if SLC_tab is not None:
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

        if OFF_par1 is not None:
            result = self._validate(Path(OFF_par1).exists(), result)
        if OFF_par2 is not None:
            result = self._validate(Path(OFF_par2).exists(), result)
        if OFF_par3 is not None:
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

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if MLI is not None:
            Path(MLI).touch()
        if MLI_par is not None:
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

        if SLC_tab1 is not None:
            result = self._validate(Path(SLC_tab1).exists(), result)
        if SLC_tab2 is not None:
            result = self._validate(Path(SLC_tab2).exists(), result)
        if SLC_tab3 is not None:
            result = self._validate(Path(SLC_tab3).exists(), result)
        return result

    def par_RSAT_SLC(self, CEOS_leader: str, SLC_par: str, CEOS_data: str, SLC: str = None, sc_dB = None, dt = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RSAT_SLC", supplied_args))

        if "par_RSAT_SLC" in self.call_count:
            self.call_count["par_RSAT_SLC"] += 1
        else:
            self.call_count["par_RSAT_SLC"] = 1

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if SLC is not None:
            Path(SLC).touch()
        return result

    def base_est_fft(self, interf: str, SLC1_par: str, OFF_par: str, baseline: str, nazfft = None, r_samp = None, az_line = None, nrfft = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_est_fft", supplied_args))

        if "base_est_fft" in self.call_count:
            self.call_count["base_est_fft"] += 1
        else:
            self.call_count["base_est_fft"] = 1

        if interf is not None:
            result = self._validate(Path(interf).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if baseline is not None:
            Path(baseline).touch()
        return result

    def par_RCM_GRC(self, RCM_dir: str, polarization, radcal, noise, SLC_par: str = None, SLC: str = None, GRC_par: str = None, GRC: str = None, rps = None, noise_pwr: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RCM_GRC", supplied_args))

        if "par_RCM_GRC" in self.call_count:
            self.call_count["par_RCM_GRC"] += 1
        else:
            self.call_count["par_RCM_GRC"] = 1

        if RCM_dir is not None:
            result = self._validate(Path(RCM_dir).exists(), result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(radcal in valid_values, result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        if GRC_par is not None:
            Path(GRC_par).touch()
        if GRC is not None:
            Path(GRC).touch()
        if noise_pwr is not None:
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

        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if baseline is not None:
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

        if HDF5 is not None:
            result = self._validate(Path(HDF5).exists(), result)
        if trunk is not None:
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

        if SLC_par is not None:
            Path(SLC_par).touch()
        return result

    def af_SLC(self, SLC_par: str, SLC: str, rwin = None, azwin = None, dr = None, daz = None, thres = None, a1_flg = None, b0_flg = None, offsets: str = None, n_ovr = None, roff = None, azoff = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "af_SLC", supplied_args))

        if "af_SLC" in self.call_count:
            self.call_count["af_SLC"] += 1
        else:
            self.call_count["af_SLC"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        valid_values = [0, 1]
        result = self._validate(a1_flg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(b0_flg in valid_values, result)
        if offsets is not None:
            Path(offsets).touch()
        return result

    def par_EORC_JERS_SLC(self, CEOS_SAR_leader: str, SLC_par: str, CEOS_data: str = None, slc: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_EORC_JERS_SLC", supplied_args))

        if "par_EORC_JERS_SLC" in self.call_count:
            self.call_count["par_EORC_JERS_SLC"] += 1
        else:
            self.call_count["par_EORC_JERS_SLC"] = 1

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if slc is not None:
            Path(slc).touch()
        return result

    def fill_gaps(self, data_in: str, width, data_out: str, dtype = None, method = None, max_dist = None, bp_flag = None, win = None, ds_method = None, ds_size = None, ds_data: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "fill_gaps", supplied_args))

        if "fill_gaps" in self.call_count:
            self.call_count["fill_gaps"] += 1
        else:
            self.call_count["fill_gaps"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        if ds_data is not None:
            Path(ds_data).touch()
        return result

    def rascc_mask_thinning(self, ras_in: str, in_file: str, width, ras_out: str, nmax = None, thresh_1 = None, thresh_nmax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "rascc_mask_thinning", supplied_args))

        if "rascc_mask_thinning" in self.call_count:
            self.call_count["rascc_mask_thinning"] += 1
        else:
            self.call_count["rascc_mask_thinning"] = 1

        if ras_in is not None:
            result = self._validate(Path(ras_in).exists(), result)
        if in_file is not None:
            result = self._validate(Path(in_file).exists(), result)
        if ras_out is not None:
            Path(ras_out).touch()
        return result

    def GRD_to_SR(self, GRD_par: str, SLC_par: str, OFF_par: str, in_file: str, out_file: str, rlks = None, azlks = None, interp_mode = None, sr_rsp = None, sr_azsp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "GRD_to_SR", supplied_args))

        if "GRD_to_SR" in self.call_count:
            self.call_count["GRD_to_SR"] += 1
        else:
            self.call_count["GRD_to_SR"] = 1

        if GRD_par is not None:
            result = self._validate(Path(GRD_par).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if in_file is not None:
            result = self._validate(Path(in_file).exists(), result)
        if out_file is not None:
            Path(out_file).touch()
        valid_values = [0, 1, 2]
        result = self._validate(interp_mode in valid_values, result)
        return result

    def multi_look2(self, SLC: str, SLC_par: str, MLI: str, MLI_par: str, r_dec, az_dec, rwin = None, azwin = None, wflg = None, lanczos = None, beta = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_look2", supplied_args))

        if "multi_look2" in self.call_count:
            self.call_count["multi_look2"] += 1
        else:
            self.call_count["multi_look2"] = 1

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if MLI is not None:
            Path(MLI).touch()
        if MLI_par is not None:
            Path(MLI_par).touch()
        valid_values = [0, 1]
        result = self._validate(wflg in valid_values, result)
        return result

    def dcomp_sirc(self, infile: str, outfile: str, samples, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "dcomp_sirc", supplied_args))

        if "dcomp_sirc" in self.call_count:
            self.call_count["dcomp_sirc"] += 1
        else:
            self.call_count["dcomp_sirc"] = 1

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
            Path(outfile).touch()
        return result

    def subtract_phase(self, interf_in: str, phase_file: str, interf_out: str, width, factor = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "subtract_phase", supplied_args))

        if "subtract_phase" in self.call_count:
            self.call_count["subtract_phase"] += 1
        else:
            self.call_count["subtract_phase"] = 1

        if interf_in is not None:
            result = self._validate(Path(interf_in).exists(), result)
        if phase_file is not None:
            result = self._validate(Path(phase_file).exists(), result)
        if interf_out is not None:
            Path(interf_out).touch()
        return result

    def dcomp_sirc_quad(self, infile: str, outfile: str, samples, parameter, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "dcomp_sirc_quad", supplied_args))

        if "dcomp_sirc_quad" in self.call_count:
            self.call_count["dcomp_sirc_quad"] += 1
        else:
            self.call_count["dcomp_sirc_quad"] = 1

        if infile is not None:
            result = self._validate(Path(infile).exists(), result)
        if outfile is not None:
            Path(outfile).touch()
        valid_values = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        result = self._validate(parameter in valid_values, result)
        return result

    def SLC_copy_ScanSAR(self, SLC1_tab: str, SLC2_tab: str, BURST_tab: str, dtype = None, SLC2_dir = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_copy_ScanSAR", supplied_args))

        if "SLC_copy_ScanSAR" in self.call_count:
            self.call_count["SLC_copy_ScanSAR"] += 1
        else:
            self.call_count["SLC_copy_ScanSAR"] = 1

        if SLC1_tab is not None:
            result = self._validate(Path(SLC1_tab).exists(), result)
        if SLC2_tab is not None and not Path(SLC2_tab).exists():
            Path(SLC2_tab).touch()
        if BURST_tab is not None:
            result = self._validate(Path(BURST_tab).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def par_RISAT_SLC(self, CEOS_leader: str, BAND_META: str, SLC_par: str, CEOS_image: str, SLC: str = None, line_dir = None, pix_dir = None, cal_flg = None, KdB = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RISAT_SLC", supplied_args))

        if "par_RISAT_SLC" in self.call_count:
            self.call_count["par_RISAT_SLC"] += 1
        else:
            self.call_count["par_RISAT_SLC"] = 1

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if BAND_META is not None:
            result = self._validate(Path(BAND_META).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if CEOS_image is not None:
            result = self._validate(Path(CEOS_image).exists(), result)
        if SLC is not None:
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

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if SLC_shift is not None:
            Path(SLC_shift).touch()
        if SLC_shift_par is not None:
            Path(SLC_shift_par).touch()
        return result

    def DELFT_vec2(self, SLC_par: str, DELFT_dir, nstate = None, interval = None, ODR = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "DELFT_vec2", supplied_args))

        if "DELFT_vec2" in self.call_count:
            self.call_count["DELFT_vec2"] += 1
        else:
            self.call_count["DELFT_vec2"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        return result

    def SLC_intf(self, SLC_1: str, SLC_2R: str, SLC1_par: str, SLC2R_par: str, OFF_par: str, interf: str, rlks, azlks, loff = None, nlines = None, sps_flg = None, azf_flg = None, rp1_flg = None, rp2_flg = None, SLC_1s = None, SLC_2Rs = None, SLC_1s_par = None, SLC_2Rs_par = None, az_beta = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_intf", supplied_args))

        if "SLC_intf" in self.call_count:
            self.call_count["SLC_intf"] += 1
        else:
            self.call_count["SLC_intf"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2R is not None:
            result = self._validate(Path(SLC_2R).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2R_par is not None:
            result = self._validate(Path(SLC2R_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if interf is not None:
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

    def ptarg_cal_SLC(self, SLC_par: str, SLC: str, r_samp, az_samp, psigma, c_r_samp, c_az_samp, ptr_image: str, r_plot: str, az_plot: str, pcal: str, osf = None, win = None, pltflg = None, psz = None, csz = None, c_image: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ptarg_cal_SLC", supplied_args))

        if "ptarg_cal_SLC" in self.call_count:
            self.call_count["ptarg_cal_SLC"] += 1
        else:
            self.call_count["ptarg_cal_SLC"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if ptr_image is not None:
            Path(ptr_image).touch()
        if r_plot is not None:
            Path(r_plot).touch()
        if az_plot is not None:
            Path(az_plot).touch()
        if pcal is not None:
            Path(pcal).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        if c_image is not None:
            Path(c_image).touch()
        return result

    def sbi_filt(self, SLC_1: str, SLC1_par: str, SLC2R_par: str, SLCf: str, SLCf_par: str, SLCb: str, SLCb_par: str, norm_sq, iwflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "sbi_filt", supplied_args))

        if "sbi_filt" in self.call_count:
            self.call_count["sbi_filt"] += 1
        else:
            self.call_count["sbi_filt"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2R_par is not None:
            result = self._validate(Path(SLC2R_par).exists(), result)
        if SLCf is not None:
            Path(SLCf).touch()
        if SLCf_par is not None:
            Path(SLCf_par).touch()
        if SLCb is not None:
            Path(SLCb).touch()
        if SLCb_par is not None:
            Path(SLCb_par).touch()
        valid_values = [0, 1]
        result = self._validate(iwflg in valid_values, result)
        return result

    def ph_slope_base(self, int_in: str, SLC_par: str, OFF_par: str, base: str, int_out: str, int_type = None, inverse = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ph_slope_base", supplied_args))

        if "ph_slope_base" in self.call_count:
            self.call_count["ph_slope_base"] += 1
        else:
            self.call_count["ph_slope_base"] = 1

        if int_in is not None:
            result = self._validate(Path(int_in).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if base is not None:
            result = self._validate(Path(base).exists(), result)
        if int_out is not None:
            Path(int_out).touch()
        return result

    def multi_look_ScanSAR(self, SLC_tab: str, MLI: str, MLI_par: str, rlks, azlks, bflg = None, SLCR_tab: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_look_ScanSAR", supplied_args))

        if "multi_look_ScanSAR" in self.call_count:
            self.call_count["multi_look_ScanSAR"] += 1
        else:
            self.call_count["multi_look_ScanSAR"] = 1

        if SLC_tab is not None:
            result = self._validate(Path(SLC_tab).exists(), result)
        if MLI is not None:
            Path(MLI).touch()
        if MLI_par is not None:
            Path(MLI_par).touch()
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        if SLCR_tab is not None:
            result = self._validate(Path(SLCR_tab).exists(), result)
        return result

    def par_ASNARO2(self, CEOS_data: str, CEOS_leader: str, SLC_par: str, SLC: str = None, reramp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ASNARO2", supplied_args))

        if "par_ASNARO2" in self.call_count:
            self.call_count["par_ASNARO2"] += 1
        else:
            self.call_count["par_ASNARO2"] = 1

        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        valid_values = [0, 1]
        result = self._validate(reramp in valid_values, result)
        return result

    def grasses(self, int: str, flag: str, unw: str, width, xmin = None, xmax = None, ymin = None, ymax = None, xinit = None, yinit = None, init_ph = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "grasses", supplied_args))

        if "grasses" in self.call_count:
            self.call_count["grasses"] += 1
        else:
            self.call_count["grasses"] = 1

        if int is not None:
            result = self._validate(Path(int).exists(), result)
        if flag is not None:
            result = self._validate(Path(flag).exists(), result)
        if unw is not None:
            Path(unw).touch()
        return result

    def mask_data(self, data_in: str, width, data_out: str, mask: str, dtype = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "mask_data", supplied_args))

        if "mask_data" in self.call_count:
            self.call_count["mask_data"] += 1
        else:
            self.call_count["mask_data"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        if mask is not None:
            result = self._validate(Path(mask).exists(), result)
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def par_ESA_JERS_SEASAT_SLC(self, CEOS_data: str, CEOS_leader: str, SLC_par: str, SLC: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ESA_JERS_SEASAT_SLC", supplied_args))

        if "par_ESA_JERS_SEASAT_SLC" in self.call_count:
            self.call_count["par_ESA_JERS_SEASAT_SLC"] += 1
        else:
            self.call_count["par_ESA_JERS_SEASAT_SLC"] = 1

        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        return result

    def ScanSAR_burst_overlap(self, SLC_tab: str, root_name: str, rlks, azlks, mode = None, bflg = None, SLCR_tab: str = None, dburst = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_burst_overlap", supplied_args))

        if "ScanSAR_burst_overlap" in self.call_count:
            self.call_count["ScanSAR_burst_overlap"] += 1
        else:
            self.call_count["ScanSAR_burst_overlap"] = 1

        if SLC_tab is not None:
            result = self._validate(Path(SLC_tab).exists(), result)
        if root_name is not None:
            Path(root_name).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        if SLCR_tab is not None:
            result = self._validate(Path(SLCR_tab).exists(), result)
        return result

    def ORB_prop_SLC(self, SLC_par: str, nstate = None, interval = None, extra = None, mode = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ORB_prop_SLC", supplied_args))

        if "ORB_prop_SLC" in self.call_count:
            self.call_count["ORB_prop_SLC"] += 1
        else:
            self.call_count["ORB_prop_SLC"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        return result

    def interp_ad(self, data_in: str, data_out: str, width, r_max = None, np_min = None, np_max = None, w_mode = None, dtype = None, cp_data = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "interp_ad", supplied_args))

        if "interp_ad" in self.call_count:
            self.call_count["interp_ad"] += 1
        else:
            self.call_count["interp_ad"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(w_mode in valid_values, result)
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(cp_data in valid_values, result)
        return result

    def par_RISAT_GRD(self, CEOS_leader: str, BAND_META: str, GRD_par: str, CEOS_image: str, GRD: str = None, line_dir = None, pix_dir = None, cal_flg = None, KdB = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RISAT_GRD", supplied_args))

        if "par_RISAT_GRD" in self.call_count:
            self.call_count["par_RISAT_GRD"] += 1
        else:
            self.call_count["par_RISAT_GRD"] = 1

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if BAND_META is not None:
            result = self._validate(Path(BAND_META).exists(), result)
        if GRD_par is not None:
            Path(GRD_par).touch()
        if CEOS_image is not None:
            result = self._validate(Path(CEOS_image).exists(), result)
        if GRD is not None:
            Path(GRD).touch()
        valid_values = [0, 1]
        result = self._validate(cal_flg in valid_values, result)
        return result

    def par_RSAT_SCW(self, CEOS_leader: str, CEOS_trailer: str, CEOS_data: str, GRD_par: str, GRD: str, sc_dB = None, dt = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RSAT_SCW", supplied_args))

        if "par_RSAT_SCW" in self.call_count:
            self.call_count["par_RSAT_SCW"] += 1
        else:
            self.call_count["par_RSAT_SCW"] = 1

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if CEOS_trailer is not None:
            result = self._validate(Path(CEOS_trailer).exists(), result)
        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if GRD_par is not None:
            Path(GRD_par).touch()
        if GRD is not None:
            Path(GRD).touch()
        return result

    def neutron(self, intensity: str, flag: str, width, n_thres, ymin = None, ymax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "neutron", supplied_args))

        if "neutron" in self.call_count:
            self.call_count["neutron"] += 1
        else:
            self.call_count["neutron"] = 1

        if intensity is not None:
            result = self._validate(Path(intensity).exists(), result)
        if flag is not None:
            result = self._validate(Path(flag).exists(), result)
        return result

    def SLC_mosaic_S1_TOPS(self, SLC_tab: str, SLC: str, SLC_par: str, rlks, azlks, bflg = None, SLCR_tab: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_mosaic_S1_TOPS", supplied_args))

        if "SLC_mosaic_S1_TOPS" in self.call_count:
            self.call_count["SLC_mosaic_S1_TOPS"] += 1
        else:
            self.call_count["SLC_mosaic_S1_TOPS"] = 1

        if SLC_tab is not None:
            result = self._validate(Path(SLC_tab).exists(), result)
        if SLC is not None:
            Path(SLC).touch()
        if SLC_par is not None:
            Path(SLC_par).touch()
        valid_values = [0, 1]
        result = self._validate(bflg in valid_values, result)
        if SLCR_tab is not None:
            result = self._validate(Path(SLCR_tab).exists(), result)
        return result

    def multi_look(self, SLC: str, SLC_par: str, MLI: str, MLI_par: str, rlks, azlks, loff = None, nlines = None, scale = None, exp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_look", supplied_args))

        if "multi_look" in self.call_count:
            self.call_count["multi_look"] += 1
        else:
            self.call_count["multi_look"] = 1

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if MLI is not None:
            Path(MLI).touch()
        if MLI_par is not None:
            Path(MLI_par).touch()
        return result

    def mosaic_WB(self, data_tab: str, dtype: str, data_out: str, data_par_out: str, sc_flg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "mosaic_WB", supplied_args))

        if "mosaic_WB" in self.call_count:
            self.call_count["mosaic_WB"] += 1
        else:
            self.call_count["mosaic_WB"] = 1

        if data_tab is not None:
            result = self._validate(Path(data_tab).exists(), result)
        if dtype is not None:
            result = self._validate(Path(dtype).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        if data_par_out is not None:
            Path(data_par_out).touch()
        valid_values = [0, 1]
        result = self._validate(sc_flg in valid_values, result)
        return result

    def ScanSAR_burst_to_mosaic(self, DATA_tab: str, mosaic: str, MLI_par: str, mflg = None, dtype = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_burst_to_mosaic", supplied_args))

        if "ScanSAR_burst_to_mosaic" in self.call_count:
            self.call_count["ScanSAR_burst_to_mosaic"] += 1
        else:
            self.call_count["ScanSAR_burst_to_mosaic"] = 1

        if DATA_tab is not None:
            result = self._validate(Path(DATA_tab).exists(), result)
        if mosaic is not None:
            Path(mosaic).touch()
        if MLI_par is not None:
            Path(MLI_par).touch()
        valid_values = [0, 1, 2]
        result = self._validate(mflg in valid_values, result)
        return result

    def par_UAVSAR_SLC(self, ann: str, SLC_MLC_in: str = None, SLC_MLI_par: str = None, SLC_MLI_out: str = None, image_type = None, image_format = None, DOP: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_UAVSAR_SLC", supplied_args))

        if "par_UAVSAR_SLC" in self.call_count:
            self.call_count["par_UAVSAR_SLC"] += 1
        else:
            self.call_count["par_UAVSAR_SLC"] = 1

        if ann is not None:
            result = self._validate(Path(ann).exists(), result)
        if SLC_MLC_in is not None:
            result = self._validate(Path(SLC_MLC_in).exists(), result)
        if SLC_MLI_par is not None:
            Path(SLC_MLI_par).touch()
        if SLC_MLI_out is not None:
            Path(SLC_MLI_out).touch()
        valid_values = [0, 2]
        result = self._validate(image_format in valid_values, result)
        if DOP is not None:
            result = self._validate(Path(DOP).exists(), result)
        return result

    def ScanSAR_burst_copy(self, SLC: str, SLC_par: str, TOPS_par: str, SLC_out: str, SLC_out_par: str, burst_num, drflg = None, SLC_par2: str = None, dtype = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ScanSAR_burst_copy", supplied_args))

        if "ScanSAR_burst_copy" in self.call_count:
            self.call_count["ScanSAR_burst_copy"] += 1
        else:
            self.call_count["ScanSAR_burst_copy"] = 1

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if TOPS_par is not None:
            result = self._validate(Path(TOPS_par).exists(), result)
        if SLC_out is not None:
            Path(SLC_out).touch()
        if SLC_out_par is not None:
            Path(SLC_out_par).touch()
        valid_values = [0, 1]
        result = self._validate(drflg in valid_values, result)
        if SLC_par2 is not None:
            Path(SLC_par2).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def hgt_map(self, unw: str, SLC_par: str, OFF_par: str, baseline: str, hgt: str, gr: str, ph_flag = None, loff = None, nlines = None, SLC2R_par = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "hgt_map", supplied_args))

        if "hgt_map" in self.call_count:
            self.call_count["hgt_map"] += 1
        else:
            self.call_count["hgt_map"] = 1

        if unw is not None:
            result = self._validate(Path(unw).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if baseline is not None:
            result = self._validate(Path(baseline).exists(), result)
        if hgt is not None:
            Path(hgt).touch()
        if gr is not None:
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

        if HDF5 is not None:
            result = self._validate(Path(HDF5).exists(), result)
        if trunk is not None:
            Path(trunk).touch()
        return result

    def par_TX_GRD(self, annotation_XML: str, COSAR, GRD_par, GRD: str, pol = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_TX_GRD", supplied_args))

        if "par_TX_GRD" in self.call_count:
            self.call_count["par_TX_GRD"] += 1
        else:
            self.call_count["par_TX_GRD"] = 1

        if annotation_XML is not None:
            result = self._validate(Path(annotation_XML).exists(), result)
        if GRD is not None:
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

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_par_in is not None:
            result = self._validate(Path(data_par_in).exists(), result)
        if data_tab is not None:
            result = self._validate(Path(data_tab).exists(), result)
        if dtype is not None:
            result = self._validate(Path(dtype).exists(), result)
        return result

    def base_init(self, SLC1_par: str, SLC2_par: str, OFF_par: str, interf: str, base, mflag = None, nrfft = None, nazfft = None, r_samp = None, az_line = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_init", supplied_args))

        if "base_init" in self.call_count:
            self.call_count["base_init"] += 1
        else:
            self.call_count["base_init"] = 1

        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if interf is not None:
            result = self._validate(Path(interf).exists(), result)
        valid_values = [0, 1, 2, 3, 4]
        result = self._validate(mflag in valid_values, result)
        return result

    def par_SIRC(self, CEOS_leader: str, SLC_par: str, UTC_MET = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_SIRC", supplied_args))

        if "par_SIRC" in self.call_count:
            self.call_count["par_SIRC"] += 1
        else:
            self.call_count["par_SIRC"] = 1

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        return result

    def rascc_mask(self, cc: str, pwr: str, width, start_cc = None, start_pwr = None, nlines = None, pixavr = None, pixavaz = None, cc_thres = None, pwr_thres = None, cc_min = None, cc_max = None, scale = None, exp = None, LR = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "rascc_mask", supplied_args))

        if "rascc_mask" in self.call_count:
            self.call_count["rascc_mask"] += 1
        else:
            self.call_count["rascc_mask"] = 1

        if cc is not None:
            result = self._validate(Path(cc).exists(), result)
        if pwr is not None:
            result = self._validate(Path(pwr).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def init_offset_orbit(self, SLC1_par: str, SLC2_par: str, OFF_par: str, rpos = None, azpos = None, cflag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "init_offset_orbit", supplied_args))

        if "init_offset_orbit" in self.call_count:
            self.call_count["init_offset_orbit"] += 1
        else:
            self.call_count["init_offset_orbit"] = 1

        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None and not Path(OFF_par).exists():
            Path(OFF_par).touch()
        valid_values = [0, 1]
        result = self._validate(cflag in valid_values, result)
        return result

    def interf_SLC(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, MLI_1: str, MLI_2: str, interf, nrlk = None, nazlk = None, loff = None, nltot = None, rfilt = None, azfilt = None, s_off = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "interf_SLC", supplied_args))

        if "interf_SLC" in self.call_count:
            self.call_count["interf_SLC"] += 1
        else:
            self.call_count["interf_SLC"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if MLI_1 is not None:
            Path(MLI_1).touch()
        if MLI_2 is not None:
            Path(MLI_2).touch()
        valid_values = [0, 1]
        result = self._validate(rfilt in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(azfilt in valid_values, result)
        return result

    def MLI_cat(self, MLI_1: str, MLI_2: str, MLI1_par: str, MLI2_par: str, MLI_3: str, MLI3_par: str, dtype = None, mflg = None, overlap = None, interp_mode = None, degree = None, extrapol = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "MLI_cat", supplied_args))

        if "MLI_cat" in self.call_count:
            self.call_count["MLI_cat"] += 1
        else:
            self.call_count["MLI_cat"] = 1

        if MLI_1 is not None:
            result = self._validate(Path(MLI_1).exists(), result)
        if MLI_2 is not None:
            result = self._validate(Path(MLI_2).exists(), result)
        if MLI1_par is not None:
            result = self._validate(Path(MLI1_par).exists(), result)
        if MLI2_par is not None:
            result = self._validate(Path(MLI2_par).exists(), result)
        if MLI_3 is not None:
            Path(MLI_3).touch()
        if MLI3_par is not None:
            Path(MLI3_par).touch()
        valid_values = [0, 1]
        result = self._validate(mflg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(extrapol in valid_values, result)
        return result

    def par_RCM_SLC(self, RCM_dir: str, polarization, radcal, noise, SLC_par: str, SLC: str, noise_pwr: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_RCM_SLC", supplied_args))

        if "par_RCM_SLC" in self.call_count:
            self.call_count["par_RCM_SLC"] += 1
        else:
            self.call_count["par_RCM_SLC"] = 1

        if RCM_dir is not None:
            result = self._validate(Path(RCM_dir).exists(), result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(radcal in valid_values, result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        if noise_pwr is not None:
            Path(noise_pwr).touch()
        return result

    def phase_slope(self, interf: str, slopes: str, width, win_sz = None, thres = None, xmin = None, xmax = None, ymin = None, ymax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "phase_slope", supplied_args))

        if "phase_slope" in self.call_count:
            self.call_count["phase_slope"] += 1
        else:
            self.call_count["phase_slope"] = 1

        if interf is not None:
            result = self._validate(Path(interf).exists(), result)
        if slopes is not None:
            Path(slopes).touch()
        return result

    def par_TX_ScanSAR(self, annot_XML: str, swath, SLC_par: str, SLC: str, TOPS_par: str, bwflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_TX_ScanSAR", supplied_args))

        if "par_TX_ScanSAR" in self.call_count:
            self.call_count["par_TX_ScanSAR"] += 1
        else:
            self.call_count["par_TX_ScanSAR"] = 1

        if annot_XML is not None:
            result = self._validate(Path(annot_XML).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        if TOPS_par is not None:
            Path(TOPS_par).touch()
        valid_values = [0, 1]
        result = self._validate(bwflg in valid_values, result)
        return result

    def ave_image(self, im_list: str, width, ave: str, start = None, nlines = None, pixav_x = None, pixav_y = None, zflag = None, nmin = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ave_image", supplied_args))

        if "ave_image" in self.call_count:
            self.call_count["ave_image"] += 1
        else:
            self.call_count["ave_image"] = 1

        if im_list is not None:
            result = self._validate(Path(im_list).exists(), result)
        if ave is not None:
            Path(ave).touch()
        valid_values = [0, 1]
        result = self._validate(zflag in valid_values, result)
        return result

    def multi_cpx(self, data_in: str, OFF_par_in: str, data_out: str, OFF_par_out: str, rlks = None, azlks = None, loff = None, nlines = None, roff = None, nsamp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_cpx", supplied_args))

        if "multi_cpx" in self.call_count:
            self.call_count["multi_cpx"] += 1
        else:
            self.call_count["multi_cpx"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if OFF_par_in is not None:
            result = self._validate(Path(OFF_par_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        if OFF_par_out is not None and not Path(OFF_par_out).exists():
            Path(OFF_par_out).touch()
        return result

    def ASAR_LO_phase_drift(self, SLC1_par: str, SLC2_par: str, OFF_par: str, ph_drift: str):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ASAR_LO_phase_drift", supplied_args))

        if "ASAR_LO_phase_drift" in self.call_count:
            self.call_count["ASAR_LO_phase_drift"] += 1
        else:
            self.call_count["ASAR_LO_phase_drift"] = 1

        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if ph_drift is not None:
            Path(ph_drift).touch()
        return result

    def radcal_pwr_stat(self, SLC_tab: str, SLC_tab_cal: str, plist: str, MSR_cal, PWR_cal, roff = None, loff = None, nr = None, nl = None, plist_out = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "radcal_pwr_stat", supplied_args))

        if "radcal_pwr_stat" in self.call_count:
            self.call_count["radcal_pwr_stat"] += 1
        else:
            self.call_count["radcal_pwr_stat"] = 1

        if SLC_tab is not None:
            result = self._validate(Path(SLC_tab).exists(), result)
        if SLC_tab_cal is not None:
            result = self._validate(Path(SLC_tab_cal).exists(), result)
        if plist is not None:
            result = self._validate(Path(plist).exists(), result)
        return result

    def par_ICEYE_SLC(self, HDF5: str, SLC_par: str, slc: str = None, dtype = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ICEYE_SLC", supplied_args))

        if "par_ICEYE_SLC" in self.call_count:
            self.call_count["par_ICEYE_SLC"] += 1
        else:
            self.call_count["par_ICEYE_SLC"] = 1

        if HDF5 is not None:
            result = self._validate(Path(HDF5).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if slc is not None:
            Path(slc).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def offset_SLC_tracking(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, offs: str, snr: str, rsw = None, azsw = None, offsets: str = None, n_ovr = None, thres = None, rstep = None, azstep = None, rstart = None, rstop = None, azstart = None, azstop = None, ISZ = None, pflag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_SLC_tracking", supplied_args))

        if "offset_SLC_tracking" in self.call_count:
            self.call_count["offset_SLC_tracking"] += 1
        else:
            self.call_count["offset_SLC_tracking"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if offs is not None:
            Path(offs).touch()
        if snr is not None:
            Path(snr).touch()
        if offsets is not None:
            Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        return result

    def tree_cc(self, flag: str, width, mbl = None, xmin = None, xmax = None, ymin = None, ymax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "tree_cc", supplied_args))

        if "tree_cc" in self.call_count:
            self.call_count["tree_cc"] += 1
        else:
            self.call_count["tree_cc"] = 1

        if flag is not None:
            result = self._validate(Path(flag).exists(), result)
        return result

    def MLI_copy(self, MLI_in: str, MLI_in_par: str, MLI_out: str, MLI_out_par: str, roff = None, nr = None, loff = None, nl = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "MLI_copy", supplied_args))

        if "MLI_copy" in self.call_count:
            self.call_count["MLI_copy"] += 1
        else:
            self.call_count["MLI_copy"] = 1

        if MLI_in is not None:
            result = self._validate(Path(MLI_in).exists(), result)
        if MLI_in_par is not None:
            result = self._validate(Path(MLI_in_par).exists(), result)
        if MLI_out is not None:
            Path(MLI_out).touch()
        if MLI_out_par is not None:
            Path(MLI_out_par).touch()
        return result

    def ORRM_vec(self, SLC_par: str, ORRM: str, nstate = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ORRM_vec", supplied_args))

        if "ORRM_vec" in self.call_count:
            self.call_count["ORRM_vec"] += 1
        else:
            self.call_count["ORRM_vec"] = 1

        if SLC_par is not None and not Path(SLC_par).exists():
            Path(SLC_par).touch()
        if ORRM is not None:
            result = self._validate(Path(ORRM).exists(), result)
        return result

    def SLC_ovr(self, SLC: str, SLC_par: str, SLC_ovr: str, SLC_ovr_par: str, r_ovr = None, az_ovr = None, mode = None, order = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_ovr", supplied_args))

        if "SLC_ovr" in self.call_count:
            self.call_count["SLC_ovr"] += 1
        else:
            self.call_count["SLC_ovr"] = 1

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if SLC_ovr is not None:
            Path(SLC_ovr).touch()
        if SLC_ovr_par is not None:
            Path(SLC_ovr_par).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def tree_gzw(self, flag: str, width, mbl = None, xmin = None, xmax = None, ymin = None, ymax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "tree_gzw", supplied_args))

        if "tree_gzw" in self.call_count:
            self.call_count["tree_gzw"] += 1
        else:
            self.call_count["tree_gzw"] = 1

        if flag is not None:
            result = self._validate(Path(flag).exists(), result)
        return result

    def mcf(self, interf: str, wgt: str, mask: str, unw: str, width, tri_mode = None, roff = None, loff = None, nr = None, nlines = None, npat_r = None, npat_az = None, ovrlap = None, r_init = None, az_init = None, init_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "mcf", supplied_args))

        if "mcf" in self.call_count:
            self.call_count["mcf"] += 1
        else:
            self.call_count["mcf"] = 1

        if interf is not None:
            result = self._validate(Path(interf).exists(), result)
        if wgt is not None:
            result = self._validate(Path(wgt).exists(), result)
        if mask is not None:
            result = self._validate(Path(mask).exists(), result)
        if unw is not None:
            Path(unw).touch()
        valid_values = [0, 1]
        result = self._validate(tri_mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(init_flag in valid_values, result)
        return result

    def par_ESA_ERS(self, CEOS_SAR_leader: str, SLC_par: str, CEOS_DAT: str = None, SLC: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ESA_ERS", supplied_args))

        if "par_ESA_ERS" in self.call_count:
            self.call_count["par_ESA_ERS"] += 1
        else:
            self.call_count["par_ESA_ERS"] = 1

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if CEOS_DAT is not None:
            result = self._validate(Path(CEOS_DAT).exists(), result)
        if SLC is not None:
            Path(SLC).touch()
        return result

    def SLC_interp(self, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, SLC_2R: str, SLC2R_par: str, loff = None, nlines = None, mode = None, order = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_interp", supplied_args))

        if "SLC_interp" in self.call_count:
            self.call_count["SLC_interp"] += 1
        else:
            self.call_count["SLC_interp"] = 1

        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if SLC_2R is not None:
            Path(SLC_2R).touch()
        if SLC2R_par is not None:
            Path(SLC2R_par).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        return result

    def par_S1_SLC(self, GeoTIFF: str, annotation_XML: str, calibration_XML: str, noise_XML: str, SLC_par: str, SLC: str, TOPS_par: str = None, dtype = None, sc_dB = None, noise_pwr = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_S1_SLC", supplied_args))

        if "par_S1_SLC" in self.call_count:
            self.call_count["par_S1_SLC"] += 1
        else:
            self.call_count["par_S1_SLC"] = 1

        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if annotation_XML is not None:
            result = self._validate(Path(annotation_XML).exists(), result)
        if calibration_XML is not None:
            result = self._validate(Path(calibration_XML).exists(), result)
        if noise_XML is not None:
            result = self._validate(Path(noise_XML).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        if TOPS_par is not None:
            Path(TOPS_par).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def par_ASAR(self, ASAR_ERS_file: str = None, output_name: str = None, K_dB = None, to = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_ASAR", supplied_args))

        if "par_ASAR" in self.call_count:
            self.call_count["par_ASAR"] += 1
        else:
            self.call_count["par_ASAR"] = 1

        if ASAR_ERS_file is not None:
            result = self._validate(Path(ASAR_ERS_file).exists(), result)
        if output_name is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if SLC_par is not None:
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

        if DATA is not None:
            result = self._validate(Path(DATA).exists(), result)
        return result

    def base_ls(self, SLC_par: str, OFF_par: str, gcp_ph: str, baseline: str, ph_flag = None, bc_flag = None, bn_flag = None, bcdot_flag = None, bndot_flag = None, bperp_min = None, SLC2R_par: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "base_ls", supplied_args))

        if "base_ls" in self.call_count:
            self.call_count["base_ls"] += 1
        else:
            self.call_count["base_ls"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if gcp_ph is not None:
            result = self._validate(Path(gcp_ph).exists(), result)
        if baseline is not None:
            result = self._validate(Path(baseline).exists(), result)
        if SLC2R_par is not None:
            result = self._validate(Path(SLC2R_par).exists(), result)
        return result

    def az_spec_SLC(self, SLC: str, SLC_par: str, spectrum: str, roff = None, namb = None, pltflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "az_spec_SLC", supplied_args))

        if "az_spec_SLC" in self.call_count:
            self.call_count["az_spec_SLC"] += 1
        else:
            self.call_count["az_spec_SLC"] = 1

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if spectrum is not None:
            Path(spectrum).touch()
        valid_values = [0, 1]
        result = self._validate(pltflg in valid_values, result)
        return result

    def SLC_copy(self, SLC_in: str, SLC_par_in: str, SLC_out: str, SLC_par_out: str, fcase = None, sc = None, roff = None, nr = None, loff = None, nl = None, swap = None, header_lines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_copy", supplied_args))

        if "SLC_copy" in self.call_count:
            self.call_count["SLC_copy"] += 1
        else:
            self.call_count["SLC_copy"] = 1

        if SLC_in is not None:
            result = self._validate(Path(SLC_in).exists(), result)
        if SLC_par_in is not None:
            result = self._validate(Path(SLC_par_in).exists(), result)
        if SLC_out is not None:
            Path(SLC_out).touch()
        if SLC_par_out is not None:
            Path(SLC_par_out).touch()
        valid_values = [1, 2, 3, 4]
        result = self._validate(fcase in valid_values, result)
        valid_values = [0, 1, 2]
        result = self._validate(swap in valid_values, result)
        return result

    def az_integrate(self, data: str, width: str, azi: str, cflg, scale = None, lz = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "az_integrate", supplied_args))

        if "az_integrate" in self.call_count:
            self.call_count["az_integrate"] += 1
        else:
            self.call_count["az_integrate"] = 1

        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if width is not None:
            result = self._validate(Path(width).exists(), result)
        if azi is not None:
            Path(azi).touch()
        valid_values = [0, 1]
        result = self._validate(cflg in valid_values, result)
        return result

    def SLC_cat(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, SLC_3: str, SLC3_par: str, dopflg = None, iflg = None, phflg = None, gainflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_cat", supplied_args))

        if "SLC_cat" in self.call_count:
            self.call_count["SLC_cat"] += 1
        else:
            self.call_count["SLC_cat"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if SLC_3 is not None:
            Path(SLC_3).touch()
        if SLC3_par is not None:
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

    def par_NovaSAR_SLC(self, GeoTIFF: str, XML: str, polarization, SLC_par: str, SLC: str = None, dtype = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_NovaSAR_SLC", supplied_args))

        if "par_NovaSAR_SLC" in self.call_count:
            self.call_count["par_NovaSAR_SLC"] += 1
        else:
            self.call_count["par_NovaSAR_SLC"] = 1

        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if XML is not None:
            result = self._validate(Path(XML).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        return result

    def SLC_corners(self, SLC_par: str, terra_alt: str = None, kml: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_corners", supplied_args))

        if "SLC_corners" in self.call_count:
            self.call_count["SLC_corners"] += 1
        else:
            self.call_count["SLC_corners"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if terra_alt is not None:
            result = self._validate(Path(terra_alt).exists(), result)
        if kml is not None:
            Path(kml).touch()
        return result

    def SLC_deramp(self, SLC_1: str, SLC_par1: str, SLC_2: str, SLC_par2: str, mode, dop_ph: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_deramp", supplied_args))

        if "SLC_deramp" in self.call_count:
            self.call_count["SLC_deramp"] += 1
        else:
            self.call_count["SLC_deramp"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_par1 is not None:
            result = self._validate(Path(SLC_par1).exists(), result)
        if SLC_2 is not None:
            Path(SLC_2).touch()
        if SLC_par2 is not None:
            Path(SLC_par2).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        if dop_ph is not None:
            Path(dop_ph).touch()
        return result

    def residue(self, int: str, flag: str, width, xmin = None, xmax = None, ymin = None, ymax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "residue", supplied_args))

        if "residue" in self.call_count:
            self.call_count["residue"] += 1
        else:
            self.call_count["residue"] = 1

        if int is not None:
            result = self._validate(Path(int).exists(), result)
        if flag is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PRI_par is not None:
            Path(PRI_par).touch()
        if CEOS_DAT is not None:
            result = self._validate(Path(CEOS_DAT).exists(), result)
        if PRI is not None:
            Path(PRI).touch()
        return result

    def create_offset(self, SLC1_par: str, SLC2_par: str, OFF_par: str, algorithm = None, rlks = None, azlks = None, iflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "create_offset", supplied_args))

        if "create_offset" in self.call_count:
            self.call_count["create_offset"] += 1
        else:
            self.call_count["create_offset"] = 1

        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None and not Path(OFF_par).exists():
            Path(OFF_par).touch()
        valid_values = [1, 2]
        result = self._validate(algorithm in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(iflg in valid_values, result)
        return result

    def multi_look_MLI(self, MLI_in: str, MLI_in_par: str, MLI_out: str, MLI_out_par: str, rlks, azlks, loff = None, nlines = None, scale = None, e_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_look_MLI", supplied_args))

        if "multi_look_MLI" in self.call_count:
            self.call_count["multi_look_MLI"] += 1
        else:
            self.call_count["multi_look_MLI"] = 1

        if MLI_in is not None:
            result = self._validate(Path(MLI_in).exists(), result)
        if MLI_in_par is not None:
            result = self._validate(Path(MLI_in_par).exists(), result)
        if MLI_out is not None:
            Path(MLI_out).touch()
        if MLI_out_par is not None:
            Path(MLI_out_par).touch()
        valid_values = [0, 1]
        result = self._validate(e_flag in valid_values, result)
        return result

    def multi_real(self, data_in: str, OFF_par_in: str, data_out: str, OFF_par_out: str, rlks = None, azlks = None, loff = None, nlines = None, roff = None, nsamp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "multi_real", supplied_args))

        if "multi_real" in self.call_count:
            self.call_count["multi_real"] += 1
        else:
            self.call_count["multi_real"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if OFF_par_in is not None:
            result = self._validate(Path(OFF_par_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        if OFF_par_out is not None and not Path(OFF_par_out).exists():
            Path(OFF_par_out).touch()
        return result

    def SLC_intf2(self, SLC_1: str, SLC_2R: str, SLC1_par: str, SLC2R_par: str, MLI_1: str, MLI_2R: str, MLI1_par: str, MLI2R_par: str, interf: str, cc: str, r_dec, az_dec, rwin = None, azwin = None, wflg = None, lanczos = None, beta = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "SLC_intf2", supplied_args))

        if "SLC_intf2" in self.call_count:
            self.call_count["SLC_intf2"] += 1
        else:
            self.call_count["SLC_intf2"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2R is not None:
            result = self._validate(Path(SLC_2R).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2R_par is not None:
            result = self._validate(Path(SLC2R_par).exists(), result)
        if MLI_1 is not None:
            Path(MLI_1).touch()
        if MLI_2R is not None:
            Path(MLI_2R).touch()
        if MLI1_par is not None:
            Path(MLI1_par).touch()
        if MLI2R_par is not None:
            Path(MLI2R_par).touch()
        if interf is not None:
            Path(interf).touch()
        if cc is not None:
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

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if GRD_par is not None:
            Path(GRD_par).touch()
        if GRD is not None:
            Path(GRD).touch()
        return result

    def offset_pwr(self, SLC1: str, SLC2: str, SLC1_par: str, SLC2_par: str, OFF_par: str, offs: str, ccp: str, rwin = None, azwin = None, offsets: str = None, n_ovr = None, nr = None, naz = None, thres = None, lanczos = None, bw_frac = None, deramp = None, int_filt = None, pflag = None, pltflg = None, ccs: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "offset_pwr", supplied_args))

        if "offset_pwr" in self.call_count:
            self.call_count["offset_pwr"] += 1
        else:
            self.call_count["offset_pwr"] = 1

        if SLC1 is not None:
            result = self._validate(Path(SLC1).exists(), result)
        if SLC2 is not None:
            result = self._validate(Path(SLC2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if OFF_par is not None:
            result = self._validate(Path(OFF_par).exists(), result)
        if offs is not None:
            Path(offs).touch()
        if ccp is not None:
            Path(ccp).touch()
        if offsets is not None:
            Path(offsets).touch()
        valid_values = [0, 1]
        result = self._validate(deramp in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(int_filt in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(pflag in valid_values, result)
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        if ccs is not None:
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

        if CEOS_Image is not None:
            result = self._validate(Path(CEOS_Image).exists(), result)
        if SLC_par is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if PRI_par is not None:
            Path(PRI_par).touch()
        if CEOS_DAT is not None:
            result = self._validate(Path(CEOS_DAT).exists(), result)
        if PRI is not None:
            Path(PRI).touch()
        return result

    def par_KC_PALSAR_slr(self, facter_m: str, CEOS_leader: str, SLC_par: str, pol, pls_mode, KC_data: str, pwr: str, fdtab: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_KC_PALSAR_slr", supplied_args))

        if "par_KC_PALSAR_slr" in self.call_count:
            self.call_count["par_KC_PALSAR_slr"] += 1
        else:
            self.call_count["par_KC_PALSAR_slr"] = 1

        if facter_m is not None:
            result = self._validate(Path(facter_m).exists(), result)
        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        valid_values = [1, 2, 3]
        result = self._validate(pls_mode in valid_values, result)
        if KC_data is not None:
            result = self._validate(Path(KC_data).exists(), result)
        if pwr is not None:
            Path(pwr).touch()
        if fdtab is not None:
            Path(fdtab).touch()
        return result

    def ptarg_SLC(self, SLC_par: str, SLC: str, r_samp, az_samp, ptr_image: str, r_plot: str, az_plot: str, ptr_par: str = None, osf = None, win = None, pltflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "ptarg_SLC", supplied_args))

        if "ptarg_SLC" in self.call_count:
            self.call_count["ptarg_SLC"] += 1
        else:
            self.call_count["ptarg_SLC"] = 1

        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if ptr_image is not None:
            Path(ptr_image).touch()
        if r_plot is not None:
            Path(r_plot).touch()
        if az_plot is not None:
            Path(az_plot).touch()
        if ptr_par is not None:
            Path(ptr_par).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(pltflg in valid_values, result)
        return result

    def par_EORC_PALSAR(self, CEOS_leader: str, SLC_par: str, CEOS_data: str, SLC: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_EORC_PALSAR", supplied_args))

        if "par_EORC_PALSAR" in self.call_count:
            self.call_count["par_EORC_PALSAR"] += 1
        else:
            self.call_count["par_EORC_PALSAR"] = 1

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if CEOS_data is not None:
            result = self._validate(Path(CEOS_data).exists(), result)
        if SLC is not None:
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

        if annotation_XML is not None:
            result = self._validate(Path(annotation_XML).exists(), result)
        return result

    def par_GF3_SLC(self, GeoTIFF: str, annotation_XML: str, SLC_par: str, SLC: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "par_GF3_SLC", supplied_args))

        if "par_GF3_SLC" in self.call_count:
            self.call_count["par_GF3_SLC"] += 1
        else:
            self.call_count["par_GF3_SLC"] = 1

        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if annotation_XML is not None:
            result = self._validate(Path(annotation_XML).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        return result

    def cc_wave(self, interf: str, MLI_1: str, MLI_2: str, cc: str, width, bx = None, by = None, wflg = None, xmin = None, xmax = None, ymin = None, ymax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "cc_wave", supplied_args))

        if "cc_wave" in self.call_count:
            self.call_count["cc_wave"] += 1
        else:
            self.call_count["cc_wave"] = 1

        if interf is not None:
            result = self._validate(Path(interf).exists(), result)
        if MLI_1 is not None:
            result = self._validate(Path(MLI_1).exists(), result)
        if MLI_2 is not None:
            result = self._validate(Path(MLI_2).exists(), result)
        if cc is not None:
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

        if product_XML is not None:
            result = self._validate(Path(product_XML).exists(), result)
        if lut_XML is not None:
            result = self._validate(Path(lut_XML).exists(), result)
        if GeoTIFF is not None:
            result = self._validate(Path(GeoTIFF).exists(), result)
        if polarization is not None:
            result = self._validate(Path(polarization).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        if SLC is not None:
            Path(SLC).touch()
        return result

    def residue_cc(self, int: str, flag: str, width, xmin = None, xmax = None, ymin = None, ymax = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "residue_cc", supplied_args))

        if "residue_cc" in self.call_count:
            self.call_count["residue_cc"] += 1
        else:
            self.call_count["residue_cc"] = 1

        if int is not None:
            result = self._validate(Path(int).exists(), result)
        if flag is not None:
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

        if CEOS_SAR_leader is not None:
            result = self._validate(Path(CEOS_SAR_leader).exists(), result)
        if SLC_par is not None:
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

        if CEOS_leader is not None:
            result = self._validate(Path(CEOS_leader).exists(), result)
        if CEOS_trailer is not None:
            result = self._validate(Path(CEOS_trailer).exists(), result)
        if SLC_par is not None:
            Path(SLC_par).touch()
        return result

    def fspf(self, data_in: str, data_out: str, width, dtype = None, r_max = None, spf_type = None, MLI_par = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "fspf", supplied_args))

        if "fspf" in self.call_count:
            self.call_count["fspf"] += 1
        else:
            self.call_count["fspf"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        valid_values = [0, 1, 2]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1, 2, 3, 4, 5]
        result = self._validate(spf_type in valid_values, result)
        return result

    def radcal_SLC(self, SLC: str, SLC_PAR: str, CSLC: str, CSLC_PAR: str, fcase = None, antenna = None, rloss_flag = None, ant_flag = None, refarea_flag = None, sc_dB = None, K_dB = None, pix_area: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("ISP", "radcal_SLC", supplied_args))

        if "radcal_SLC" in self.call_count:
            self.call_count["radcal_SLC"] += 1
        else:
            self.call_count["radcal_SLC"] = 1

        if SLC is not None:
            result = self._validate(Path(SLC).exists(), result)
        if SLC_PAR is not None:
            result = self._validate(Path(SLC_PAR).exists(), result)
        if CSLC is not None:
            Path(CSLC).touch()
        if CSLC_PAR is not None:
            Path(CSLC_PAR).touch()
        valid_values = [1, 2, 3, 4]
        result = self._validate(fcase in valid_values, result)
        if pix_area is not None:
            Path(pix_area).touch()
        return result

    def line_interp(self, input = None, output = None, width = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "line_interp", supplied_args))

        if "line_interp" in self.call_count:
            self.call_count["line_interp"] += 1
        else:
            self.call_count["line_interp"] = 1

        return result

    def product_cpx(self, f1: str, f2: str, f_out: str, width, start = None, nlines = None, conjg_flg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "product_cpx", supplied_args))

        if "product_cpx" in self.call_count:
            self.call_count["product_cpx"] += 1
        else:
            self.call_count["product_cpx"] = 1

        if f1 is not None:
            result = self._validate(Path(f1).exists(), result)
        if f2 is not None:
            result = self._validate(Path(f2).exists(), result)
        if f_out is not None:
            Path(f_out).touch()
        valid_values = [0, 1]
        result = self._validate(conjg_flg in valid_values, result)
        return result

    def ras_majority(self, ras_in: str, ras_out: str, filter_width = None, LR = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_majority", supplied_args))

        if "ras_majority" in self.call_count:
            self.call_count["ras_majority"] += 1
        else:
            self.call_count["ras_majority"] = 1

        if ras_in is not None:
            result = self._validate(Path(ras_in).exists(), result)
        if ras_out is not None:
            Path(ras_out).touch()
        return result

    def average_filter(self, din: str, dout: str, width, bx, by = None, wflg = None, min_pt = None, zflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "average_filter", supplied_args))

        if "average_filter" in self.call_count:
            self.call_count["average_filter"] += 1
        else:
            self.call_count["average_filter"] = 1

        if din is not None:
            result = self._validate(Path(din).exists(), result)
        if dout is not None:
            Path(dout).touch()
        valid_values = [0, 1, 2]
        result = self._validate(wflg in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def mask_class(self, class_map: str, file_in: str, file_out: str, format_flag, LR, selection_flag, n_class, class_1, class_n = None, null_value = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "mask_class", supplied_args))

        if "mask_class" in self.call_count:
            self.call_count["mask_class"] += 1
        else:
            self.call_count["mask_class"] = 1

        if class_map is not None:
            result = self._validate(Path(class_map).exists(), result)
        if file_in is not None:
            result = self._validate(Path(file_in).exists(), result)
        if file_out is not None:
            Path(file_out).touch()
        valid_values = [0, 1, 2, 3]
        result = self._validate(format_flag in valid_values, result)
        return result

    def ras_ratio_dB(self, pwr1: str, pwr2: str, width, start_pwr1 = None, start_pwr2 = None, nlines = None, pixavr = None, pixavaz = None, min_cc = None, max_cc = None, scale = None, exp = None, LR = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_ratio_dB", supplied_args))

        if "ras_ratio_dB" in self.call_count:
            self.call_count["ras_ratio_dB"] += 1
        else:
            self.call_count["ras_ratio_dB"] = 1

        if pwr1 is not None:
            result = self._validate(Path(pwr1).exists(), result)
        if pwr2 is not None:
            result = self._validate(Path(pwr2).exists(), result)
        if rasf is not None:
            Path(rasf).touch()
        return result

    def linear_to_dB(self, data_in: str, data_out: str, width, inverse_flag = None, null_value = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "linear_to_dB", supplied_args))

        if "linear_to_dB" in self.call_count:
            self.call_count["linear_to_dB"] += 1
        else:
            self.call_count["linear_to_dB"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        valid_values = [0, 1]
        result = self._validate(inverse_flag in valid_values, result)
        return result

    def histogram(self, data_in: str, width, polygon: str, hist: str, stat: str, min, max, nbins = None, mode = None, lin_log = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "histogram", supplied_args))

        if "histogram" in self.call_count:
            self.call_count["histogram"] += 1
        else:
            self.call_count["histogram"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if polygon is not None:
            result = self._validate(Path(polygon).exists(), result)
        if hist is not None:
            Path(hist).touch()
        if stat is not None:
            Path(stat).touch()
        valid_values = [0, 1]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(lin_log in valid_values, result)
        return result

    def gamma_map(self, input_data: str, output_data: str, width, nlooks, bx, by = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "gamma_map", supplied_args))

        if "gamma_map" in self.call_count:
            self.call_count["gamma_map"] += 1
        else:
            self.call_count["gamma_map"] = 1

        if input_data is not None:
            result = self._validate(Path(input_data).exists(), result)
        if output_data is not None:
            Path(output_data).touch()
        return result

    def m_chi(self, s0: str, m: str, s2chi: str, S_par: str, c1: str, c2: str = None, c3: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "m-chi", supplied_args))

        if "m-chi" in self.call_count:
            self.call_count["m-chi"] += 1
        else:
            self.call_count["m-chi"] = 1

        if s0 is not None:
            result = self._validate(Path(s0).exists(), result)
        if m is not None:
            result = self._validate(Path(m).exists(), result)
        if s2chi is not None:
            result = self._validate(Path(s2chi).exists(), result)
        if S_par is not None:
            result = self._validate(Path(S_par).exists(), result)
        if c1 is not None:
            Path(c1).touch()
        if c2 is not None:
            Path(c2).touch()
        if c3 is not None:
            Path(c3).touch()
        return result

    def reallks(self, image: str, ML_image: str, width, rlks = None, azlks = None, start = None, nlines = None, r_start = None, nsamp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "reallks", supplied_args))

        if "reallks" in self.call_count:
            self.call_count["reallks"] += 1
        else:
            self.call_count["reallks"] = 1

        if image is not None:
            result = self._validate(Path(image).exists(), result)
        if ML_image is not None:
            Path(ML_image).touch()
        return result

    def frost(self, pwr1: str, pwr1_frost: str, width, fx = None, sx = None, power = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "frost", supplied_args))

        if "frost" in self.call_count:
            self.call_count["frost"] += 1
        else:
            self.call_count["frost"] = 1

        if pwr1 is not None:
            result = self._validate(Path(pwr1).exists(), result)
        if pwr1_frost is not None:
            Path(pwr1_frost).touch()
        return result

    def mt_lee_filt(self, im_list: str, ref_image: str, width, winsz, L_ref, L, cthres, out_list: str, ref_out: str = None, b_coeff: str = None, filt_num: str = None, msr: str = None, ctr: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "mt_lee_filt", supplied_args))

        if "mt_lee_filt" in self.call_count:
            self.call_count["mt_lee_filt"] += 1
        else:
            self.call_count["mt_lee_filt"] = 1

        if im_list is not None:
            result = self._validate(Path(im_list).exists(), result)
        if ref_image is not None:
            result = self._validate(Path(ref_image).exists(), result)
        if out_list is not None:
            result = self._validate(Path(out_list).exists(), result)
        if ref_out is not None:
            Path(ref_out).touch()
        if b_coeff is not None:
            Path(b_coeff).touch()
        if filt_num is not None:
            Path(filt_num).touch()
        if msr is not None:
            Path(msr).touch()
        if ctr is not None:
            Path(ctr).touch()
        return result

    def multi_stat(self, im_list: str, width, im_out: str, mode, rank, nmin = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "multi_stat", supplied_args))

        if "multi_stat" in self.call_count:
            self.call_count["multi_stat"] += 1
        else:
            self.call_count["multi_stat"] = 1

        if im_list is not None:
            result = self._validate(Path(im_list).exists(), result)
        if im_out is not None:
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

        if data1 is not None:
            result = self._validate(Path(data1).exists(), result)
        valid_values = [1]
        result = self._validate(func in valid_values, result)
        if data2 is not None:
            Path(data2).touch()
        return result

    def temp_filt(self, data_tab: str, width, waz = None, wr = None, wt_flag = None, zero_flag = None, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "temp_filt", supplied_args))

        if "temp_filt" in self.call_count:
            self.call_count["temp_filt"] += 1
        else:
            self.call_count["temp_filt"] = 1

        if data_tab is not None:
            result = self._validate(Path(data_tab).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(wt_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def cpxlks(self, CMPLX: str, ML_CMPLX: str, width, rlks = None, azlks = None, start = None, nlines = None, r_start = None, nsamp = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "cpxlks", supplied_args))

        if "cpxlks" in self.call_count:
            self.call_count["cpxlks"] += 1
        else:
            self.call_count["cpxlks"] = 1

        if CMPLX is not None:
            result = self._validate(Path(CMPLX).exists(), result)
        if ML_CMPLX is not None:
            Path(ML_CMPLX).touch()
        return result

    def ras_to_rgb(self, red_channel: str, green_channel: str, blue_channel: str, ras_out: str, LR = None, null_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_to_rgb", supplied_args))

        if "ras_to_rgb" in self.call_count:
            self.call_count["ras_to_rgb"] += 1
        else:
            self.call_count["ras_to_rgb"] = 1

        if red_channel is not None:
            result = self._validate(Path(red_channel).exists(), result)
        if green_channel is not None:
            result = self._validate(Path(green_channel).exists(), result)
        if blue_channel is not None:
            result = self._validate(Path(blue_channel).exists(), result)
        if ras_out is not None:
            Path(ras_out).touch()
        valid_values = [0, 1]
        result = self._validate(null_flag in valid_values, result)
        return result

    def ave2pwr(self, pwr1: str, pwr2: str, pwr_out: str, width, scale_factor = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ave2pwr", supplied_args))

        if "ave2pwr" in self.call_count:
            self.call_count["ave2pwr"] += 1
        else:
            self.call_count["ave2pwr"] = 1

        if pwr1 is not None:
            result = self._validate(Path(pwr1).exists(), result)
        if pwr2 is not None:
            result = self._validate(Path(pwr2).exists(), result)
        if pwr_out is not None:
            Path(pwr_out).touch()
        return result

    def ras_m_chi(self, s1: str, c1: str, c2: str, c3: str, width, start = None, nlines = None, pixavr = None, pixavaz = None, scale = None, exp = None, rasf: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_m-chi", supplied_args))

        if "ras_m-chi" in self.call_count:
            self.call_count["ras_m-chi"] += 1
        else:
            self.call_count["ras_m-chi"] = 1

        if s1 is not None:
            result = self._validate(Path(s1).exists(), result)
        if c1 is not None:
            result = self._validate(Path(c1).exists(), result)
        if c2 is not None:
            result = self._validate(Path(c2).exists(), result)
        if c3 is not None:
            result = self._validate(Path(c3).exists(), result)
        if rasf is not None:
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

        if pwr1 is not None:
            result = self._validate(Path(pwr1).exists(), result)
        if inc is not None:
            result = self._validate(Path(inc).exists(), result)
        if gamma is not None:
            Path(gamma).touch()
        return result

    def hsi_color_scale(self, file_out: str, nval = None, chip_width = None, gap = None, height = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "hsi_color_scale", supplied_args))

        if "hsi_color_scale" in self.call_count:
            self.call_count["hsi_color_scale"] += 1
        else:
            self.call_count["hsi_color_scale"] = 1

        if file_out is not None:
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

        if unw is not None:
            result = self._validate(Path(unw).exists(), result)
        if cpx is not None:
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

        if SLC_HH is not None:
            result = self._validate(Path(SLC_HH).exists(), result)
        if SLC_VV is not None:
            result = self._validate(Path(SLC_VV).exists(), result)
        if SLC_HV is not None:
            result = self._validate(Path(SLC_HV).exists(), result)
        if SLC_HH_par is not None:
            result = self._validate(Path(SLC_HH_par).exists(), result)
        if SLC_VV_par is not None:
            result = self._validate(Path(SLC_VV_par).exists(), result)
        if SLC_HV_par is not None:
            result = self._validate(Path(SLC_HV_par).exists(), result)
        if P is not None:
            Path(P).touch()
        return result

    def m_alpha(self, s0: str, m: str, alpha: str, S_par: str, c1: str, c2: str = None, c3: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "m-alpha", supplied_args))

        if "m-alpha" in self.call_count:
            self.call_count["m-alpha"] += 1
        else:
            self.call_count["m-alpha"] = 1

        if s0 is not None:
            result = self._validate(Path(s0).exists(), result)
        if m is not None:
            result = self._validate(Path(m).exists(), result)
        if alpha is not None:
            result = self._validate(Path(alpha).exists(), result)
        if S_par is not None:
            result = self._validate(Path(S_par).exists(), result)
        if c1 is not None:
            Path(c1).touch()
        if c2 is not None:
            Path(c2).touch()
        if c3 is not None:
            Path(c3).touch()
        return result

    def stokes(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, S: str, S_par: str, rlks, azlks, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "stokes", supplied_args))

        if "stokes" in self.call_count:
            self.call_count["stokes"] += 1
        else:
            self.call_count["stokes"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if S is not None:
            Path(S).touch()
        if S_par is not None:
            Path(S_par).touch()
        return result

    def histogram_ras(self, ras_in: str, polygon, histograms, mean_stdev = None, percent = None, lr_flag = None, start = None, stop = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "histogram_ras", supplied_args))

        if "histogram_ras" in self.call_count:
            self.call_count["histogram_ras"] += 1
        else:
            self.call_count["histogram_ras"] = 1

        if ras_in is not None:
            result = self._validate(Path(ras_in).exists(), result)
        return result

    def product(self, data_1: str, data_2: str, product: str, width, bx = None, by = None, wgt_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "product", supplied_args))

        if "product" in self.call_count:
            self.call_count["product"] += 1
        else:
            self.call_count["product"] = 1

        if data_1 is not None:
            result = self._validate(Path(data_1).exists(), result)
        if data_2 is not None:
            result = self._validate(Path(data_2).exists(), result)
        if product is not None:
            Path(product).touch()
        valid_values = [0, 1, 2]
        result = self._validate(wgt_flag in valid_values, result)
        return result

    def mt_lee_filt_cpx(self, cpx_list: str, ref_image: str, width, winsz, L_ref, cthres, out_list: str, ref_out: str = None, b_coeff: str = None, filt_num: str = None, msr: str = None, ctr: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "mt_lee_filt_cpx", supplied_args))

        if "mt_lee_filt_cpx" in self.call_count:
            self.call_count["mt_lee_filt_cpx"] += 1
        else:
            self.call_count["mt_lee_filt_cpx"] = 1

        if cpx_list is not None:
            result = self._validate(Path(cpx_list).exists(), result)
        if ref_image is not None:
            result = self._validate(Path(ref_image).exists(), result)
        if out_list is not None:
            result = self._validate(Path(out_list).exists(), result)
        if ref_out is not None:
            Path(ref_out).touch()
        if b_coeff is not None:
            Path(b_coeff).touch()
        if filt_num is not None:
            Path(filt_num).touch()
        if msr is not None:
            Path(msr).touch()
        if ctr is not None:
            Path(ctr).touch()
        return result

    def polcovar(self, SLC_1: str, SLC_2: str, SLC_3: str, SLC1_par: str, SLC2_par: str, SLC3_par: str, C: str, C_par: str, rlks, azlks, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "polcovar", supplied_args))

        if "polcovar" in self.call_count:
            self.call_count["polcovar"] += 1
        else:
            self.call_count["polcovar"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC_3 is not None:
            result = self._validate(Path(SLC_3).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if SLC3_par is not None:
            result = self._validate(Path(SLC3_par).exists(), result)
        if C is not None:
            Path(C).touch()
        if C_par is not None:
            Path(C_par).touch()
        return result

    def stokes_qm(self, S: str, S_par: str, m: str = None, s2chi: str = None, s2psi: str = None, m_l: str = None, m_c: str = None, lp_ratio: str = None, cp_ratio: str = None, mu: str = None, delta: str = None, alpha: str = None, phi: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "stokes_qm", supplied_args))

        if "stokes_qm" in self.call_count:
            self.call_count["stokes_qm"] += 1
        else:
            self.call_count["stokes_qm"] = 1

        if S is not None:
            result = self._validate(Path(S).exists(), result)
        if S_par is not None:
            result = self._validate(Path(S_par).exists(), result)
        if m is not None:
            Path(m).touch()
        if s2chi is not None:
            Path(s2chi).touch()
        if s2psi is not None:
            Path(s2psi).touch()
        if m_l is not None:
            Path(m_l).touch()
        if m_c is not None:
            Path(m_c).touch()
        if lp_ratio is not None:
            Path(lp_ratio).touch()
        if cp_ratio is not None:
            Path(cp_ratio).touch()
        if mu is not None:
            Path(mu).touch()
        if delta is not None:
            Path(delta).touch()
        if alpha is not None:
            Path(alpha).touch()
        if phi is not None:
            Path(phi).touch()
        return result

    def frame(self, data_in: str, data_out: str, width, dtype, dx1, dx2, dy1, dy2, null_flag = None, all_flag = None, null_value = None, frame_value = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "frame", supplied_args))

        if "frame" in self.call_count:
            self.call_count["frame"] += 1
        else:
            self.call_count["frame"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
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

    def polcoh(self, SLC_1: str, SLC_2: str, SLC_3: str, SLC1_par: str, SLC2_par: str, SLC3_par: str, T: str, T_par: str, rlks, azlks, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "polcoh", supplied_args))

        if "polcoh" in self.call_count:
            self.call_count["polcoh"] += 1
        else:
            self.call_count["polcoh"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC_3 is not None:
            result = self._validate(Path(SLC_3).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if SLC3_par is not None:
            result = self._validate(Path(SLC3_par).exists(), result)
        if T is not None:
            Path(T).touch()
        if T_par is not None:
            Path(T_par).touch()
        return result

    def lin_comb(self, nfiles, f1: str, f2: str, constant = None, factor1 = None, factor2 = None, f_out: str = None, width = None, start = None, nlines = None, pixav_x = None, pixav_y = None, zero_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "lin_comb", supplied_args))

        if "lin_comb" in self.call_count:
            self.call_count["lin_comb"] += 1
        else:
            self.call_count["lin_comb"] = 1

        if f1 is not None:
            result = self._validate(Path(f1).exists(), result)
        if f2 is not None:
            result = self._validate(Path(f2).exists(), result)
        if f_out is not None:
            Path(f_out).touch()
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def multi_class_mapping(self, nfiles, f1: str, f2: str, fn: str = None, classf: str = None, ras_out: str = None, width = None, start = None, nlines = None, pixav_x = None, pixav_y = None, LR = None, color_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "multi_class_mapping", supplied_args))

        if "multi_class_mapping" in self.call_count:
            self.call_count["multi_class_mapping"] += 1
        else:
            self.call_count["multi_class_mapping"] = 1

        if f1 is not None:
            result = self._validate(Path(f1).exists(), result)
        if f2 is not None:
            result = self._validate(Path(f2).exists(), result)
        if fn is not None:
            result = self._validate(Path(fn).exists(), result)
        if classf is not None:
            result = self._validate(Path(classf).exists(), result)
        if ras_out is not None:
            Path(ras_out).touch()
        valid_values = [0, 1]
        result = self._validate(color_flag in valid_values, result)
        return result

    def takecut(self, data_in: str, width, report: str, mode, pos, pr_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "takecut", supplied_args))

        if "takecut" in self.call_count:
            self.call_count["takecut"] += 1
        else:
            self.call_count["takecut"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if report is not None:
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

        if data is not None:
            result = self._validate(Path(data).exists(), result)
        if polygon is not None:
            result = self._validate(Path(polygon).exists(), result)
        if report is not None:
            Path(report).touch()
        return result

    def temp_lin_var(self, data_tab: str, mean: str, stdev: str, width, waz = None, wr = None, wt_flag = None, zero_flag = None, loff = None, nlines = None, norm_pow = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "temp_lin_var", supplied_args))

        if "temp_lin_var" in self.call_count:
            self.call_count["temp_lin_var"] += 1
        else:
            self.call_count["temp_lin_var"] = 1

        if data_tab is not None:
            result = self._validate(Path(data_tab).exists(), result)
        if mean is not None:
            Path(mean).touch()
        if stdev is not None:
            Path(stdev).touch()
        valid_values = [0, 1, 2]
        result = self._validate(wt_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def ave_cpx(self, cpx_list: str, width, ave: str, start = None, nlines = None, zflag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ave_cpx", supplied_args))

        if "ave_cpx" in self.call_count:
            self.call_count["ave_cpx"] += 1
        else:
            self.call_count["ave_cpx"] = 1

        if cpx_list is not None:
            result = self._validate(Path(cpx_list).exists(), result)
        if ave is not None:
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

        if mask_1 is not None:
            result = self._validate(Path(mask_1).exists(), result)
        if mask_2 is not None:
            result = self._validate(Path(mask_2).exists(), result)
        if mask_out is not None:
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

        if SLC_HH is not None:
            result = self._validate(Path(SLC_HH).exists(), result)
        if SLC_HV is not None:
            result = self._validate(Path(SLC_HV).exists(), result)
        if SLC_VH is not None:
            result = self._validate(Path(SLC_VH).exists(), result)
        if SLC_VV is not None:
            result = self._validate(Path(SLC_VV).exists(), result)
        if SLC_HH_par is not None:
            result = self._validate(Path(SLC_HH_par).exists(), result)
        if SLC_HV_par is not None:
            result = self._validate(Path(SLC_HV_par).exists(), result)
        if SLC_VH_par is not None:
            result = self._validate(Path(SLC_VH_par).exists(), result)
        if SLC_VV_par is not None:
            result = self._validate(Path(SLC_VV_par).exists(), result)
        if CP is not None:
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

    def single_class_mapping(self, nfiles, f1: str, lt1, ut1, fn: str = None, ltn = None, utn = None, ras_out: str = None, width = None, start = None, nlines = None, pixav_x = None, pixav_y = None, LR = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "single_class_mapping", supplied_args))

        if "single_class_mapping" in self.call_count:
            self.call_count["single_class_mapping"] += 1
        else:
            self.call_count["single_class_mapping"] = 1

        if f1 is not None:
            result = self._validate(Path(f1).exists(), result)
        if fn is not None:
            result = self._validate(Path(fn).exists(), result)
        if ras_out is not None:
            Path(ras_out).touch()
        return result

    def drawthat(self, ras_in: str, ras_out: str, pt_list: str, mode = None, r = None, g = None, b = None, xs = None, zflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "drawthat", supplied_args))

        if "drawthat" in self.call_count:
            self.call_count["drawthat"] += 1
        else:
            self.call_count["drawthat"] = 1

        if ras_in is not None:
            result = self._validate(Path(ras_in).exists(), result)
        if ras_out is not None:
            Path(ras_out).touch()
        if pt_list is not None:
            result = self._validate(Path(pt_list).exists(), result)
        valid_values = [0, 1, 2]
        result = self._validate(mode in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def m_delta(self, s0: str, m: str, delta: str, S_par: str, c1: str, c2: str = None, c3: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "m-delta", supplied_args))

        if "m-delta" in self.call_count:
            self.call_count["m-delta"] += 1
        else:
            self.call_count["m-delta"] = 1

        if s0 is not None:
            result = self._validate(Path(s0).exists(), result)
        if m is not None:
            result = self._validate(Path(m).exists(), result)
        if delta is not None:
            result = self._validate(Path(delta).exists(), result)
        if S_par is not None:
            result = self._validate(Path(S_par).exists(), result)
        if c1 is not None:
            Path(c1).touch()
        if c2 is not None:
            Path(c2).touch()
        if c3 is not None:
            Path(c3).touch()
        return result

    def temp_filt_ad(self, data_tab: str, width, zero_flag = None, loffset = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "temp_filt_ad", supplied_args))

        if "temp_filt_ad" in self.call_count:
            self.call_count["temp_filt_ad"] += 1
        else:
            self.call_count["temp_filt_ad"] = 1

        if data_tab is not None:
            result = self._validate(Path(data_tab).exists(), result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def bm3d(self, data_in: str, width, data_out: str, dtype = None, profile = None, looks = None, sigma = None, block_size = None, s_dist = None, step = None, d_max = None, t1d = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "bm3d", supplied_args))

        if "bm3d" in self.call_count:
            self.call_count["bm3d"] += 1
        else:
            self.call_count["bm3d"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
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

    def restore_float(self, input = None, output = None, width = None, interp_limit = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "restore_float", supplied_args))

        if "restore_float" in self.call_count:
            self.call_count["restore_float"] += 1
        else:
            self.call_count["restore_float"] = 1

        return result

    def haalpha(self, alpha: str, beta: str, gamma: str, SLC_par: str, anisotropy: str, entropy: str, lambda1: str, lambda2: str, lambda3: str, MLI_par: str, rlks, azlks, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "haalpha", supplied_args))

        if "haalpha" in self.call_count:
            self.call_count["haalpha"] += 1
        else:
            self.call_count["haalpha"] = 1

        if alpha is not None:
            Path(alpha).touch()
        if beta is not None:
            result = self._validate(Path(beta).exists(), result)
        if gamma is not None:
            result = self._validate(Path(gamma).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if anisotropy is not None:
            Path(anisotropy).touch()
        if entropy is not None:
            Path(entropy).touch()
        if lambda1 is not None:
            Path(lambda1).touch()
        if lambda2 is not None:
            Path(lambda2).touch()
        if lambda3 is not None:
            Path(lambda3).touch()
        if MLI_par is not None:
            Path(MLI_par).touch()
        return result

    def wolf(self, SLC_1: str, SLC_2: str, SLC1_par: str, SLC2_par: str, J: str, J_par: str, rlks, azlks, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "wolf", supplied_args))

        if "wolf" in self.call_count:
            self.call_count["wolf"] += 1
        else:
            self.call_count["wolf"] = 1

        if SLC_1 is not None:
            result = self._validate(Path(SLC_1).exists(), result)
        if SLC_2 is not None:
            result = self._validate(Path(SLC_2).exists(), result)
        if SLC1_par is not None:
            result = self._validate(Path(SLC1_par).exists(), result)
        if SLC2_par is not None:
            result = self._validate(Path(SLC2_par).exists(), result)
        if J is not None:
            Path(J).touch()
        if J_par is not None:
            Path(J_par).touch()
        return result

    def takethat_dem_par(self, data_in: str, width, positions: str, DEM_par: str, report: str, mode = None, zero_flag = None, nn_flag = None, print_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "takethat_dem_par", supplied_args))

        if "takethat_dem_par" in self.call_count:
            self.call_count["takethat_dem_par"] += 1
        else:
            self.call_count["takethat_dem_par"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if positions is not None:
            result = self._validate(Path(positions).exists(), result)
        if DEM_par is not None:
            result = self._validate(Path(DEM_par).exists(), result)
        if report is not None:
            Path(report).touch()
        return result

    def lee(self, input_data: str, output_data: str, width, nlooks, bx, by = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "lee", supplied_args))

        if "lee" in self.call_count:
            self.call_count["lee"] += 1
        else:
            self.call_count["lee"] = 1

        if input_data is not None:
            result = self._validate(Path(input_data).exists(), result)
        if output_data is not None:
            Path(output_data).touch()
        return result

    def enh_lee(self, input_data: str, output_data: str, width, nlooks, damp, bx, by = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "enh_lee", supplied_args))

        if "enh_lee" in self.call_count:
            self.call_count["enh_lee"] += 1
        else:
            self.call_count["enh_lee"] = 1

        if input_data is not None:
            result = self._validate(Path(input_data).exists(), result)
        if output_data is not None:
            Path(output_data).touch()
        return result

    def ratio(self, d1: str, d2: str, ratio: str, width, bx = None, by = None, wgt_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ratio", supplied_args))

        if "ratio" in self.call_count:
            self.call_count["ratio"] += 1
        else:
            self.call_count["ratio"] = 1

        if d1 is not None:
            result = self._validate(Path(d1).exists(), result)
        if d2 is not None:
            result = self._validate(Path(d2).exists(), result)
        if ratio is not None:
            Path(ratio).touch()
        valid_values = [0, 1, 2]
        result = self._validate(wgt_flag in valid_values, result)
        return result

    def cc_ad(self, interf: str, pwr1: str, pwr2: str, slope: str, texture: str, cc_ad: str, width, box_min = None, box_max = None, wgt_flag = None, loff = None, nl = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "cc_ad", supplied_args))

        if "cc_ad" in self.call_count:
            self.call_count["cc_ad"] += 1
        else:
            self.call_count["cc_ad"] = 1

        if interf is not None:
            result = self._validate(Path(interf).exists(), result)
        if pwr1 is not None:
            result = self._validate(Path(pwr1).exists(), result)
        if pwr2 is not None:
            result = self._validate(Path(pwr2).exists(), result)
        if slope is not None:
            result = self._validate(Path(slope).exists(), result)
        if texture is not None:
            result = self._validate(Path(texture).exists(), result)
        if cc_ad is not None:
            Path(cc_ad).touch()
        valid_values = [0, 1]
        result = self._validate(wgt_flag in valid_values, result)
        return result

    def ras_ras(self, ras_in: str, ras_out: str, col_looks = None, row_looks = None, LR = None, r_lin_log = None, g_lin_log = None, b_lin_log = None, force24 = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_ras", supplied_args))

        if "ras_ras" in self.call_count:
            self.call_count["ras_ras"] += 1
        else:
            self.call_count["ras_ras"] = 1

        if ras_in is not None:
            result = self._validate(Path(ras_in).exists(), result)
        if ras_out is not None:
            Path(ras_out).touch()
        valid_values = [0, 1]
        result = self._validate(force24 in valid_values, result)
        return result

    def lin_comb_cpx(self, nfiles, f1: str, f2: str, constant_r = None, constant_i = None, factor1_r = None, factor1_i = None, factor2_r = None, factor2_i = None, f_out: str = None, width = None, start = None, nlines = None, pixav_x = None, pixav_y = None, zero_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "lin_comb_cpx", supplied_args))

        if "lin_comb_cpx" in self.call_count:
            self.call_count["lin_comb_cpx"] += 1
        else:
            self.call_count["lin_comb_cpx"] = 1

        if f1 is not None:
            result = self._validate(Path(f1).exists(), result)
        if f2 is not None:
            result = self._validate(Path(f2).exists(), result)
        if f_out is not None:
            Path(f_out).touch()
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def median_filter(self, din: str, dout: str, width, bx, by = None, min_pt = None, zflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "median_filter", supplied_args))

        if "median_filter" in self.call_count:
            self.call_count["median_filter"] += 1
        else:
            self.call_count["median_filter"] = 1

        if din is not None:
            result = self._validate(Path(din).exists(), result)
        if dout is not None:
            Path(dout).touch()
        valid_values = [0, 1]
        result = self._validate(zflg in valid_values, result)
        return result

    def takethat(self, data_in: str, width, positions: str, report: str, mode = None, zero_flag = None, nn_flag = None, print_flag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "takethat", supplied_args))

        if "takethat" in self.call_count:
            self.call_count["takethat"] += 1
        else:
            self.call_count["takethat"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if positions is not None:
            result = self._validate(Path(positions).exists(), result)
        if report is not None:
            Path(report).touch()
        return result

    def cc_monitoring(self, nfiles, f1: str, f2: str, ras_out: str = None, width = None, cc_thresh = None, start = None, nlines = None, pixav_x = None, pixav_y = None, LR = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "cc_monitoring", supplied_args))

        if "cc_monitoring" in self.call_count:
            self.call_count["cc_monitoring"] += 1
        else:
            self.call_count["cc_monitoring"] = 1

        if f1 is not None:
            result = self._validate(Path(f1).exists(), result)
        if f2 is not None:
            result = self._validate(Path(f2).exists(), result)
        if ras_out is not None:
            Path(ras_out).touch()
        return result

    def temp_log_var(self, data_tab: str, mean: str, stdev: str, width, waz = None, wr = None, wt_flag = None, zero_flag = None, loff = None, nlines = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "temp_log_var", supplied_args))

        if "temp_log_var" in self.call_count:
            self.call_count["temp_log_var"] += 1
        else:
            self.call_count["temp_log_var"] = 1

        if data_tab is not None:
            result = self._validate(Path(data_tab).exists(), result)
        if mean is not None:
            Path(mean).touch()
        if stdev is not None:
            Path(stdev).touch()
        valid_values = [0, 1, 2]
        result = self._validate(wt_flag in valid_values, result)
        valid_values = [0, 1]
        result = self._validate(zero_flag in valid_values, result)
        return result

    def ras_to_hsi(self, HUE: str, SATURATION: str, INTENSITY: str, ras_out: str, LR = None, cflg = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "ras_to_hsi", supplied_args))

        if "ras_to_hsi" in self.call_count:
            self.call_count["ras_to_hsi"] += 1
        else:
            self.call_count["ras_to_hsi"] = 1

        if HUE is not None:
            result = self._validate(Path(HUE).exists(), result)
        if SATURATION is not None:
            result = self._validate(Path(SATURATION).exists(), result)
        if INTENSITY is not None:
            result = self._validate(Path(INTENSITY).exists(), result)
        if ras_out is not None:
            Path(ras_out).touch()
        valid_values = [0, 1]
        result = self._validate(cflg in valid_values, result)
        return result

    def edge_detection(self, data_in: str, width, data_out: str, dtype = None, op_flg = None, sigma_x = None, sigma_y = None, T1 = None, T2 = None, min_seg_size = None, max_reg_len = None, max_reg_std = None, max_reg_dist = None, seg_out: str = None, line_filt = None, max_line_std = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "edge_detection", supplied_args))

        if "edge_detection" in self.call_count:
            self.call_count["edge_detection"] += 1
        else:
            self.call_count["edge_detection"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if data_out is not None:
            Path(data_out).touch()
        valid_values = [0, 1]
        result = self._validate(dtype in valid_values, result)
        valid_values = [0, 1, 2]
        result = self._validate(op_flg in valid_values, result)
        if seg_out is not None:
            Path(seg_out).touch()
        valid_values = [0, 1]
        result = self._validate(line_filt in valid_values, result)
        return result

    def texture(self, data_in: str, format_flag, texture: str, width, type = None, bx = None, by = None, r_looks = None, az_looks = None, weights_flag = None, data_in_mean: str = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "texture", supplied_args))

        if "texture" in self.call_count:
            self.call_count["texture"] += 1
        else:
            self.call_count["texture"] = 1

        if data_in is not None:
            result = self._validate(Path(data_in).exists(), result)
        if texture is not None:
            Path(texture).touch()
        valid_values = [0, 1]
        result = self._validate(type in valid_values, result)
        valid_values = [0, 1, 2]
        result = self._validate(weights_flag in valid_values, result)
        if data_in_mean is not None:
            result = self._validate(Path(data_in_mean).exists(), result)
        return result

    def diplane_helix(self, LL: str, RR: str, SLC_par: str, diplane: str, helix: str, MLI_par: str, rlks, azlks, loff = None, nlines = None, scale = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "diplane_helix", supplied_args))

        if "diplane_helix" in self.call_count:
            self.call_count["diplane_helix"] += 1
        else:
            self.call_count["diplane_helix"] = 1

        if LL is not None:
            result = self._validate(Path(LL).exists(), result)
        if RR is not None:
            result = self._validate(Path(RR).exists(), result)
        if SLC_par is not None:
            result = self._validate(Path(SLC_par).exists(), result)
        if diplane is not None:
            Path(diplane).touch()
        if helix is not None:
            Path(helix).touch()
        if MLI_par is not None:
            Path(MLI_par).touch()
        return result

    def lin_comb_ref(self, f1: str, f2: str, constant, factor1, factor2, f_out: str, width, roff = None, loff = None, nr = None, nl = None, zflag = None):
        supplied_args = locals()
        result = (0, "", "")

        self.call_sequence.append(PyGammaCall("LAT", "lin_comb_ref", supplied_args))

        if "lin_comb_ref" in self.call_count:
            self.call_count["lin_comb_ref"] += 1
        else:
            self.call_count["lin_comb_ref"] = 1

        if f1 is not None:
            result = self._validate(Path(f1).exists(), result)
        if f2 is not None:
            result = self._validate(Path(f2).exists(), result)
        if f_out is not None:
            Path(f_out).touch()
        valid_values = [0, 1]
        result = self._validate(zflag in valid_values, result)
        return result
