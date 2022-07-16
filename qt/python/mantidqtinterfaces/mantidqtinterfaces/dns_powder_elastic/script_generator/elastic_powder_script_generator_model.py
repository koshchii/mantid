# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +

"""
DNS script generator for elastic powder data.
"""

from mantidqtinterfaces.dns_powder_elastic.data_structures.dns_dataset import \
    DNSDataset
from mantidqtinterfaces.dns_powder_tof.helpers.list_range_converters import \
    get_normation
from mantidqtinterfaces.dns_powder_tof.script_generator. \
    common_script_generator_model import DNSScriptGeneratorModel


class DNSElasticPowderScriptGeneratorModel(DNSScriptGeneratorModel):
    # pylint: disable=too-many-instance-attributes
    # having the options as instance attributes, is much better readable
    # none of them are public
    def __init__(self, parent):
        super().__init__(parent)
        self._script = []
        self._plot_list = []
        self._sample_data = None
        self._standard_data = None
        self._loop = None
        self._spac = None
        self._vanac = None
        self._nicrc = None
        self._sampb = None
        self._backfac = None
        self._ign_vana = None
        self._sum_sfnsf = None
        self._nonmag = None
        self._xyz = None
        self._corrections = None
        self._export_path = None
        self._ascii = None
        self._nexus = None
        self._norm = None

    def _check_all_fields_there(self, fields):
        ts = all(field in fields
                 for field in ['x_sf', 'y_sf', 'z_sf', 'z_nsf'])
        if not ts:
            self.raise_error("Not all fields necessary found for XYZ "
                             "analysis.")
        return ts

    @staticmethod
    def _get_xyz_field_string(sample):
        st = "# do xyz polarization analysis for magnetic powders \n"
        st += f"{sample}_nuclear_coh, {sample}_magnetic, {sample}_spin_incoh"
        st += "= xyz_separation(\n"
        st += f"    x_sf='{sample}_x_sf',"
        st += f" y_sf='{sample}_y_sf',"
        st += f" z_sf='{sample}_z_sf',"
        st += f" z_nsf='{sample}_z_nsf')"
        return [st]

    @staticmethod
    def _get_to_matrix_string(sample, sample_type):
        return (f"ConvertMDHistoToMatrixWorkspace('{sample}_{sample_type}', "
                f"Outputworkspace='mat_{sample}_{sample_type}',"
                " Normalization='NoNormalization')"
                "")

    def _get_sep_save_string(self, sample, sample_type):
        lines = ""
        if self._ascii:
            lines += (f"SaveAscii('mat_{sample}_{sample_type}', "
                      f"'{self._export_path}/{sample}_{sample_type}.csv', "
                      "WriteSpectrumID=False)"
                      "")
        if self._nexus:
            lines += (f"\nSaveNexus('mat_{sample}_{sample_type}',"
                      f" '{self._export_path}/{sample}_{sample_type}.nxs')")
        return lines

    def _get_xyz_lines(self):
        lines = []
        if not self._xyz:
            return lines
        for sample, fields in self._sample_data.data_dic.items():
            if not self._check_all_fields_there(fields):
                break
            lines += self._get_xyz_field_string(sample)
            for sample_type in ['nuclear_coh', 'magnetic', 'spin_incoh']:
                lines += [self._get_to_matrix_string(sample, sample_type)]
                self._plot_list.append(f'mat_{sample}_{sample_type}')
                lines += [self._get_sep_save_string(sample, sample_type)]
        return lines

    def _get_nsf_sf_pairs(self):
        nsf_sf_pairs = []
        for sample, fields in self._sample_data.data_dic.items():
            for field in fields.keys():
                if (field.endswith('_sf')
                        and (f'{field[:-2]}nsf' in fields.keys())):
                    nsf_sf_pairs.append(f"{sample}_{field[:-3]}")
            if not nsf_sf_pairs:
                self.raise_error("No pair of SF and NSF measurements"
                                 f" found for {sample}")
        return nsf_sf_pairs

    def _get_nonmag_lines(self):
        lines = []
        if not self._nonmag:
            return lines

        nsf_sf_pairs = self._get_nsf_sf_pairs()
        lines += [
            "# separation of coherent and incoherent scattering of non"
            " magnetic sample"
        ]
        for sample_field in sorted(nsf_sf_pairs):
            # normally only z is measured but just in case
            lines += [
                f"{sample_field}_nuclear_coh, {sample_field}_spin_incoh = "
                "non_mag_sep"
                f"('{sample_field}_sf', '{sample_field}_nsf')"
            ]
            for sample_type in ['nuclear_coh', 'spin_incoh']:
                lines += [self._get_to_matrix_string(sample_field, sample_type)]
                self._plot_list.append(f'mat_{sample_field}_{sample_type}')
                lines += [self._get_sep_save_string(sample_field, sample_type)]

        return lines

    def _setup_sample_data(self, paths, f_selector):
        self._sample_data = DNSDataset(data=f_selector['full_data'],
                                       path=paths['data_dir'],
                                       is_sample=True)
        self._plot_list = [
            'mat_' + x for x in self._sample_data.create_subtract()
        ]

    def _setup_standard_data(self, paths, f_selector):
        if self._corrections:
            self._standard_data = DNSDataset(data=f_selector['standard_data'],
                                             path=paths['standards_dir'],
                                             is_sample=False,
                                             fields=self._sample_data.fields)

    def _interpolate_standard(self):
        self._standard_data.interpolate_standard(
            banks=self._sample_data.banks,
            script_name=self._sample_data.script_name,
            parent=self)

    @staticmethod
    def _get_header_lines():
        lines = [
            "from mantidqtinterfaces.dns_powder_elastic.scripts."
            "md_powder_elastic import "
            "load_all, background_subtraction",
            "from mantidqtinterfaces.dns_powder_elastic.scripts.md_powder"
            "_elastic import "
            "vanadium_correction, flipping_ratio_correction",
            "from mantidqtinterfaces.dns_powder_elastic.scripts.md_powder"
            "_elastic import "
            "non_mag_sep, xyz_separation",
            "from mantid.simpleapi import ConvertMDHistoToMatrixWorkspace,"
            " mtd",
            "from mantid.simpleapi import SaveAscii, SaveNexus",
            "",
        ]
        return lines

    def _get_sample_data_lines(self):
        return [f'sample_data = {self._sample_data.format_dataset()}']

    def _get_standard_data_lines(self):
        if self._corrections:
            return [
                f'standard_data = {self._standard_data.format_dataset()}'
            ]
        return ['']

    def _get_binning_lines(self):
        binning = [
            self._sample_data.ttheta.bin_edge_min,
            self._sample_data.ttheta.bin_edge_max,
            self._sample_data.ttheta.nbins
        ]
        return ["", f"binning = {binning}"]

    def _get_load_data_lines(self):
        lines = [
            "wss_sample = load_all(sample_data, binning, "
            f"normalize_to='{self._norm}')"
        ]
        if self._corrections:
            lines += ["wss_standard = load_all(standard_data, "
                      f"binning, normalize_to='{self._norm}')"]
        return lines

    def _set_loop(self):
        if len(self._sample_data.keys()) == 1:
            self._loop = "for workspace in wss_sample['" \
                         f"{list(self._sample_data.keys())[0]}']:"
            self._spac = "\n" + " " * 4
        else:
            self._loop = "for sample, workspacelist in wss_sample.items(): " \
                         "\n    for workspace in workspacelist:"
            self._spac = "\n" + " " * 8

    def _get_subtract_bg_from_standard_lines(self):

        if self._vanac or self._nicrc:
            return [
                "",
                "# subtract background from vanadium and nicr",
                "for sample, workspacelist in wss_standard.items(): "
                "\n    for workspace in workspacelist:"
                "\n        background_subtraction(workspace)",
            ]
        return []

    def _get_background_string(self):
        if self._sampb:
            return f"{self._spac}background_subtraction(workspace, " \
                   f"factor={self._backfac})"
        return ""

    def _get_vana_cstring(self):
        if self._vanac:
            return f"{self._spac}vanadium_correction(workspace," \
                   " vana_set=standard_data['vana'], " \
                   f"ignore_vana_fields={self._ign_vana}, " \
                   f"sum_vana_sf_nsf={self._sum_sfnsf})"
        return ""

    def _get_sample_corrections_lines(self):
        background_string = self._get_background_string()
        vana_cstring = self._get_vana_cstring()
        if self._sampb or self._vanac:
            return [
                "", "# correct sample data",
                f"{self._loop}{background_string}{vana_cstring}"
            ]
        return []

    def _get_nicr_cor_lines(self):
        if self._nicrc:
            return [
                f"{self._loop}{self._spac}flipping_ratio_correction(workspace)"
            ]
        return []

    def _get_ascii_save_string(self):
        if self._ascii:
            return (f"{self._spac}"
                    "SaveAscii('mat_{}'.format(workspace), "
                    f"'{self._export_path}/"
                    "{}.csv'.format(workspace)"
                    ", WriteSpectrumID=False)")
        return ''

    def _get_nexus_save_string(self):
        if self._nexus:
            return (f"{self._spac}SaveNexus"
                    "('mat_{}'.format(workspace), "
                    f"'{self._export_path}/"
                    "{}.nxs'.format("
                    "workspace))")
        return ""

    def _get_convert_to_matrix_lines(self, savestring):
        lines = [
            f"{self._loop}{self._spac}ConvertMDHistoToMatrixWorkspace"
            "(workspace, Outputworkspace='mat_{}'.format(workspace), "
            f"Normalization='NoNormalization'){savestring}", ""
        ]
        return lines

    def _get_save_string(self):
        ascii_string = self._get_ascii_save_string()
        nexus_string = self._get_nexus_save_string()
        savestrings = [x for x in [ascii_string, nexus_string] if x]
        return "".join(savestrings)

    def script_maker(self, options, paths, file_selector=None):
        self._script = []

        # options shortcut
        self._vanac = options['corrections'] and options['det_efficiency']
        self._nicrc = options['corrections'] and options['flipping_ratio']
        self._sampb = (options['corrections']
                       and options['subtract_background_from_sample'])
        self._backfac = options['background_factor']
        self._ign_vana = str(options['ignore_vana_fields'])
        self._sum_sfnsf = str(options['sum_vana_sf_nsf'])
        self._nonmag = options["separation"] and options['separation_coh_inc']
        self._xyz = options["separation"] and options['separation_xyz']
        self._corrections = (self._sampb or self._vanac or self._nicrc)
        self._export_path = paths["export_dir"]
        self._ascii = (paths["ascii"] and paths["export"]
                       and bool(self._export_path))
        self._nexus = (paths["nexus"] and paths["export"]
                       and bool(self._export_path))
        self._norm = get_normation(options)

        self._setup_sample_data(paths, file_selector)
        self._setup_standard_data(paths, file_selector)
        self._interpolate_standard()
        self._set_loop()

        # starting writing script
        self._add_lines_to_script(self._get_header_lines())
        self._add_lines_to_script(self._get_sample_data_lines())
        self._add_lines_to_script(self._get_standard_data_lines())
        self._add_lines_to_script(self._get_binning_lines())
        self._add_lines_to_script(self._get_load_data_lines())
        self._add_lines_to_script(self._get_subtract_bg_from_standard_lines())
        self._add_lines_to_script(self._get_sample_corrections_lines())
        self._add_lines_to_script(self._get_nicr_cor_lines())
        save_string = self._get_save_string()
        self._add_lines_to_script(
            self._get_convert_to_matrix_lines(save_string))
        self._add_lines_to_script(self._get_xyz_lines())
        self._add_lines_to_script(self._get_nonmag_lines())

        return self._script, ''

    def get_plot_list(self):
        return self._plot_list
