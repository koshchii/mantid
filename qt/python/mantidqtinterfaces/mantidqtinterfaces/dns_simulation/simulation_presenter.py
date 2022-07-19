# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2021 ISIS Rutherford Appleton Laboratory UKRI,
#     NScD Oak Ridge National Laboratory, European Spallation Source
#     & Institut Laue - Langevin
# SPDX - License - Identifier: GPL - 3.0 +
"""
Presenter for DNS simulation
"""

from mantidqtinterfaces.dns_powder_tof.data_structures.dns_observer import \
    DNSObserver


class DNSSimulationPresenter(DNSObserver):
    TTHLIMIT = 5  # limit for 2theta difference under which entries marked

    # as matching, 5deg is detector distance at DNS

    def __init__(self, name=None, parent=None, view=None, model=None,
                 subwidgets=None):
        # pylint: disable=too-many-arguments
        super().__init__(parent=parent, name=name, view=view, model=model)
        self._refls = None
        for widget in subwidgets:
            self.view.add_tab(widget.view, widget.name, -1)

        self.own_dict['cifset'] = False
        self._set_ki()
        self.sub_presenters = [wid.presenter for wid in subwidgets]
        # connect Signals
        self.view.sig_cif_set.connect(self._cif_set)
        self.view.sig_unitcell_changed.connect(self._unitcell_changed)
        self.view.sig_wavelength_changed.connect(self._set_ki)
        self.view.sig_calculate_clicked.connect(self._calculate)
        self.view.sig_sample_rot_changed.connect(self._sample_rot_changed)

        self.view.sig_inplane_unique_switched.connect(
            self.request_to_subwidget)
        # self.view._content.lE_cif_filename.setText(
        #    'C:/Users/thomasm/Documents/YFe2O4_icsd.cif')
        # self.view.sig_cif_set.emit(
        #    'C:/Users/thomasm/Documents/YFe2O4_icsd.cif')

    def request_to_subwidget(self):
        self.get_option_dict()
        sub_dict = self._get_sub_dict()
        for presenter in self.sub_presenters:
            presenter.process_request(sub_dict)

    def back_call_from_tableitem_clicked(self, det_rot, sample_rot):
        """this function is called from subwidget table_presenter
           if the sample_rot field is clicked to set off of identified
           reflection"""
        self.get_option_dict()
        if not self.own_dict.get('fix_omega', False):
            omegaoffset = self.model.get_oof_from_ident(
                det_rot, sample_rot, self.own_dict['sample_rot'],
                self.own_dict['det_rot'])
            self.view.set_omega_offset(omegaoffset)

    def _sample_rot_changed(self):
        self.get_option_dict()
        if not self.own_dict['fix_omega']:
            self.view.set_omega_offset(0)

    def _get_and_validate(self):
        """validates hkl input, needs to be 3 numbers and hkl1 and hlk2
           must not be parallel"""
        self.get_option_dict()
        hkl1 = self.own_dict['hkl1_v']
        hkl2 = self.own_dict['hkl2_v']
        validate = self.model.validate_hkl(hkl1, hkl2)
        if not validate[0]:
            self.raise_error(validate[1])
        return validate[0]

    def _calculate(self):
        """calculates reflections list and corresponding plots """
        # self._cif_set('C:/Users/thomasm/Documents/YFe2O4_icsd.cif')
        if not self._get_and_validate():
            return
        try:
            self._refls = self.model.get_refls_and_set_orientation(
                self.own_dict)
        except ValueError:
            self.raise_error('Spacegroup not valid, use HM'
                             'Symbol or IT Number.')
            return
        self._perp_inplane()
        self._d_tooltip()
        self.request_to_subwidget()

    def _get_sub_dict(self):
        inplane = self.parent.presenter.own_dict['inplane']
        unique = self.parent.presenter.own_dict['unique']
        filtered_refls = self.model.filter_refls(self._refls, inplane, unique)
        sub_dict = {'filtered_refls': filtered_refls,
                    'tth_limit': self.TTHLIMIT,
                    'refls': self._refls,
                    'hkl1_v': self.own_dict['hkl1_v'],
                    'hkl2_p_v': self.own_dict['hkl2_p_v'],
                    'inplane_refls': self.model.return_reflections_in_map(
                        self.own_dict['hkl1_v'], self.own_dict['hkl2_p_v'],
                        self._refls),
                    'wavelength': self.own_dict['wavelength'],
                    'orilat': self.model.get_orilat(),
                    'oof_set': self.own_dict['fix_omega'],
                    }
        return sub_dict

    def _cif_set(self, filename):
        load_dict = self.model.load_cif(filename)
        self.own_dict.update(load_dict)
        self.view.set_state(load_dict)
        self.own_dict['cifset'] = True

    def _d_tooltip(self):
        self.get_option_dict()
        hkl1_v = self.own_dict['hkl1_v']
        hkl2_v = self.own_dict['hkl2_v']
        hkl2_p_v = self.own_dict['hkl2_p_v']
        d_hkl1, d_hkl2, d_hkl2_p = self.model.get_ds(hkl1_v, hkl2_v, hkl2_p_v)
        self.view.set_d_tooltip(d_hkl1, d_hkl2, d_hkl2_p)

    def _perp_inplane(self):
        """returns vector perpendicular to hkl1 in the scatteringt plane """
        q2_p = self.model.get_hkl2_p()
        self.own_dict['hkl2_p_v'] = q2_p
        self.view.set_hkl2_p(q2_p)

    def _set_ki(self):
        self.get_option_dict()
        ki = self.model.get_ki(self.own_dict['wavelength'])
        self.view.set_ki(ki)

    def _set_spacegroup(self, spacegroup):
        self.view.set_spacegroup(spacegroup)

    def _unitcell_changed(self):
        self.own_dict['cifset'] = False

    def get_option_dict(self):
        """Return own options from view"""
        if self.view is not None:
            self.own_dict.update(self.view.get_state())
        hklv_dict = self.model.get_hkl_vector_dict(self.own_dict['hkl1'],
                                                   self.own_dict['hkl2'])
        self.own_dict.update(hklv_dict)
        return self.own_dict
