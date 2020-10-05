density_xlim = {
    'fwhm_g' : [0.8, 1.5],
    'fwhm_r' : [0.8, 1.2],
    'fwhm_i' : [0.75, 1.2],
    'fwhm_z' : [0.75, 1.2],
    'airmass_g' : [1., 1.4],
    'airmass_r' : [1., 1.4],
    'airmass_i' : [1., 1.4],
    'airmass_z' : [1., 1.4],
    'depth_g' : [22, 24],
    'depth_r' : [22, 24],
    'depth_i' : [22, 24],
    'depth_z' : [22, 24],
    'skyvar_i' : [5, 15],
    'skybrite_i' : [2000, 5000],
    'det_frac_i' : [0.75, 1.],
    'sig_zp_g' : [.005, .015],
    'sig_zp_r' : [.005, .015],
    'sig_zp_i' : [.005, .015],
    'sig_zp_z' : [.005, .015],
    'exp_time_i' : [0, 800]
}

density_dx = {
    'fwhm_g' : 0.025,
    'fwhm_r' : 0.0125,
    'fwhm_i' : 0.0125,
    'fwhm_z' : 0.0125,
    'airmass_g' : .02,
    'airmass_r' : .02,
    'airmass_i' : .02,
    'airmass_z' : .02,
    'depth_g': 0.5,
    'depth_r': 0.5,
    'depth_i': 0.5,
    'depth_z': 0.5,
    'skyvar_i' : 0.25,
    'skybrite_i' : 100,
    'det_frac_i' : .05,
    'sig_zp_g' : .0005,
    'sig_zp_r' : .0005,
    'sig_zp_i' : .0005,
    'sig_zp_z' : .0005,
    'exp_time_i' : 25
}

trend_xlim = {
    'fwhm_g' : [0.8, 1.5],
    'fwhm_r' : [0.8, 1.2],
    'fwhm_i' : [0.8, 1.2],
    'fwhm_z' : [0.8, 1.2],
    'airmass_g' : [1., 1.4],
    'airmass_r' : [1., 1.4],
    'airmass_i' : [1., 1.4],
    'airmass_z' : [1., 1.4],
    'skyvar_i' : [5, 15],
    'depth_g' : [22, 24],
    'depth_r' : [22, 24],
    'depth_i' : [22, 24],
    'depth_z' : [22, 24],
    'skybrite_i' : [2000, 5000],
    'det_frac_i' : [0.75, 1.],
    'sig_zp_g' : [.005, .015],
    'sig_zp_r' : [.005, .015],
    'sig_zp_i' : [.005, .015],
    'sig_zp_z' : [.005, .015],
    'exp_time_i' : [0, 800]
}

trend_dx = {
    'fwhm_g' : 0.05,
    'fwhm_r' : 0.05,
    'fwhm_i' : 0.05,
    'fwhm_z' : 0.05,
    'airmass_g' : .05,
    'airmass_r' : .05,
    'airmass_i' : .05,
    'airmass_z' : .05,
    'depth_g': 0.5,
    'depth_r': 0.5,
    'depth_i': 0.5,
    'depth_z': 0.5,
    'skyvar_i' : 1,
    'skybrite_i' : 250,
    'det_frac_i' : .05,
    'sig_zp_g' : .00125,
    'sig_zp_r' : .00125,
    'sig_zp_i' : .00125,
    'sig_zp_z' : .00125,
    'exp_time_i' : 100
}

cols = ['meas_cm_g_1', 'meas_cm_g_2', 'meas_cm_T', 'meas_cm_fracdev', 'meas_cm_TdByTe', 'meas_cm_flux_s2n_i']

mag_colname = 'meas_cm_mag_deredden'
# mag_colname = 'meas_cm_mag'
