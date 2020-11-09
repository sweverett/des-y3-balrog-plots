# Original:
# xlim = {
#     'meas_cm_mag_deredden' : [16, 26],
#     'meas_cm_mag' : [17, 26],
#     'g-r' : [-2, 4],
#     'r-i' : [-2, 4],
#     'i-z' : [-2, 4],
#     'meas_cm_g_1' : [-1, 1],
#     'meas_cm_g_2' : [-1, 1],
#     'meas_cm_T' : [-5, 100],
#     'meas_cm_fracdev': [0, 1],
#     'meas_cm_flux_s2n_i' : [0, 100],
#     'meas_cm_TdByTe' : [0, 100]
# }

# dx = {
#     'meas_cm_mag_deredden' : 0.25,
#     'meas_cm_mag' : 0.25,
#     'g-r' : .25,
#     'r-i' : .25,
#     'i-z' : .25,
#     'meas_cm_g_1' : 0.1,
#     'meas_cm_g_2' : 0.1,
#     'meas_cm_T' : 5,
#     'meas_cm_fracdev' : 0.05,
#     'meas_cm_flux_s2n_i' : 5,
#     'meas_cm_TdByTe' : 5
# }

xlim = {
    'meas_cm_mag_deredden' : [18, 25],
    'meas_cm_mag' : [18, 25],
    'g-r' : [-0.5, 3],
    'r-i' : [-0.5, 1.5],
    'i-z' : [-2, 4],
    'meas_cm_g_1' : [-1, 1],
    'meas_cm_g_2' : [-1, 1],
    'meas_cm_T' : [0, 20],
    'meas_cm_fracdev': [0, 1],
    'meas_cm_flux_s2n_i' : [0, 150],
    'meas_cm_TdByTe' : [0, 20]
}

dx = {
    'meas_cm_mag_deredden' : 0.25,
    'meas_cm_mag' : 0.25,
    'g-r' : .2,
    'r-i' : .1,
    'i-z' : .2,
    'meas_cm_g_1' : 0.1,
    'meas_cm_g_2' : 0.1,
    'meas_cm_T' : 1,
    'meas_cm_fracdev' : 0.05,
    'meas_cm_flux_s2n_i' : 7.5,
    'meas_cm_TdByTe' : 1
}

cols = ['meas_cm_g_1', 'meas_cm_g_2', 'meas_cm_T', 'meas_cm_fracdev', 'meas_cm_TdByTe', 'meas_cm_flux_s2n_i']

mag_colname = 'meas_cm_mag_deredden'
# mag_colname = 'meas_cm_mag'

# Used for a constant number of bins in trend plots if not using dx_trends
Nxbins = 30

# percentile info
percentile_min = 1
percentile_max = 99