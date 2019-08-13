###This file generates two .txt files. One for plume, one for tail. Each with flux, escape rate, and covered plane area
import numpy as np
import h5py
import matplotlib.pyplot as plt
from general_functions import *
from spherical_flux import *
from sliceplot import *

def flux_escape_calculation(ds_key_lr, ds_key_hr, lr_plume_txt_name, lr_tail_txt_name, hr_plume_txt_name, hr_tail_txt_name):
    ds_names_all, ds_types_all = get_datasets(load_key='R2349')
    ds_keys = [ds_key_lr] #these lines establish the path to the dataset
    ds_keys_hr = [ds_key_hr] # path to high res data
    ds_names = {dsk:dsn for dsk, dsn in ds_names_all.items() if dsk in ds_keys}
    ds_types = {dsk:[dsk] for dsk in ds_keys}
    ds_names_hr = {dsk:dsn for dsk, dsn in ds_names_all.items() if dsk in ds_keys_hr}
    ds_types_hr = {dsk:[dsk] for dsk in ds_keys_hr}

    fields_num_dens = ['O2_p1_number_density']
    fields_vx = ['O2_p1_velocity_x']
    fields_vz = ['O2_p1_velocity_z']

    z_planes = np.arange(1.5,3.1,0.1)
    x_planes = np.arange(-1.5,-3.1,-0.1)
    print 'Setting up plume txt files'
    setup_vals = []
    f_lr_setup = open(lr_plume_txt_name,'ab')
    f_lr_setup_vals = (setup_vals)
    np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',', header='flux, escape rate, covered area',fmt='%10.2e', comments='#')
    f_lr_setup.close()

    f_hr_setup = open(hr_plume_txt_name,'ab')
    f_hr_setup_vals = (setup_vals)
    np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='flux, escape rate, covered area',fmt='%10.2e', comments='#')
    f_hr_setup.close()

###Plume calculation
    print 'calculating plume planes'
    for i in range(len(z_planes)):
        print 'calculating plume plane {}'.format(z_planes[i])
        x = np.linspace(-1*z_planes[i], z_planes[i], 40)
        y = np.linspace(-1*z_planes[i], z_planes[i], 40)
        z = np.repeat(z_planes[i],1600)
        meshgrid_x, meshgrid_y = np.meshgrid(x,y)
        flat_meshgrid_x = meshgrid_x.flatten()
        flat_meshgrid_y = meshgrid_y.flatten()
        coords = np.array([flat_meshgrid_x, flat_meshgrid_y, z])

        indxs = get_path_idxs(coords, ds_names, ds_types)
        indxs_hr = get_path_idxs(coords, ds_names_hr, ds_types_hr)
        vx_data_i = get_all_data(ds_names, ds_types, indxs, fields_vx)
        vx_data_i_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields_vx)
        vz_data_i = get_all_data(ds_names, ds_types, indxs, fields_vz)
        vz_data_i_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields_vz)
        num_dens_data_i = get_all_data(ds_names, ds_types, indxs, fields_num_dens)
        num_dens_data_i_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields_num_dens)

        vx_y_axis_values = vx_data_i['O2_p1_velocity_x'][ds_key_lr]
        vx_y_axis_values_hr = vx_data_i_hr['O2_p1_velocity_x'][ds_key_hr]
        vz_y_axis_values = vz_data_i['O2_p1_velocity_z'][ds_key_lr]
        vz_y_axis_values_hr = vz_data_i_hr['O2_p1_velocity_z'][ds_key_hr]
        num_dens_y_axis_values = num_dens_data_i['O2_p1_number_density'][ds_key_lr]
        num_dens_y_axis_values_hr = num_dens_data_i_hr['O2_p1_number_density'][ds_key_hr]

        flux_z_vals_lr = vz_y_axis_values * num_dens_y_axis_values * 1e5
        flux_z_vals_hr = vz_y_axis_values_hr * num_dens_y_axis_values_hr * 1e5

        cv_a = (np.abs(x[0]) + x[-1]) * (np.abs(y[0]) + y[-1])
        covered_area = cv_a * (3390.0**2) * 1e10 #covered area of slice in cm^2
        escape_rate_lr = np.mean(flux_z_vals_lr) * covered_area
        escape_rate_hr = np.mean(flux_z_vals_hr) * covered_area

        f_lr = open(lr_plume_txt_name,'ab')
        dataout_lr = np.column_stack((np.mean(flux_z_vals_lr), escape_rate_lr, cv_a))
        np.savetxt(f_lr, dataout_lr, delimiter=',',fmt='%10.2e', comments='#')
        f_lr.close()

        f_hr = open(hr_plume_txt_name,'ab')
        dataout_hr = np.column_stack((np.mean(flux_z_vals_hr), escape_rate_hr, cv_a))
        np.savetxt(f_hr, dataout_hr, delimiter=',',fmt='%10.2e', comments='#')
        f_hr.close()

###Tail Calculation
    print 'Setting up tail txt files'
    setup_vals = []
    f_lr_setup = open(lr_tail_txt_name,'ab')
    f_lr_setup_vals = (setup_vals)
    np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',', header='flux, escape rate, covered area',fmt='%10.2e', comments='#')
    f_lr_setup.close()

    f_hr_setup = open(hr_tail_txt_name,'ab')
    f_hr_setup_vals = (setup_vals)
    np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='flux, escape rate, covered area',fmt='%10.2e', comments='#')
    f_hr_setup.close()

    print 'calculating tail planes'
    for i in range(len(x_planes)):
        print 'calculating tail plane {}'.format(x_planes[i])
        x = np.repeat(x_planes[i],1600)
        y = np.linspace(x_planes[i], -1*x_planes[i], 40)
        z = np.linspace(x_planes[i], -1*x_planes[i], 40)
        meshgrid_y, meshgrid_z = np.meshgrid(y,z)
        flat_meshgrid_y = meshgrid_y.flatten()
        flat_meshgrid_z = meshgrid_z.flatten()
        coords = np.array([x, flat_meshgrid_y, flat_meshgrid_z])

        indxs = get_path_idxs(coords, ds_names, ds_types)
        indxs_hr = get_path_idxs(coords, ds_names_hr, ds_types_hr)
        vx_data_i = get_all_data(ds_names, ds_types, indxs, fields_vx)
        vx_data_i_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields_vx)
        vz_data_i = get_all_data(ds_names, ds_types, indxs, fields_vz)
        vz_data_i_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields_vz)
        num_dens_data_i = get_all_data(ds_names, ds_types, indxs, fields_num_dens)
        num_dens_data_i_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields_num_dens)

        vx_y_axis_values = vx_data_i['O2_p1_velocity_x'][ds_key_lr]
        vx_y_axis_values_hr = vx_data_i_hr['O2_p1_velocity_x'][ds_key_hr]
        vz_y_axis_values = vz_data_i['O2_p1_velocity_z'][ds_key_lr]
        vz_y_axis_values_hr = vz_data_i_hr['O2_p1_velocity_z'][ds_key_hr]
        num_dens_y_axis_values = num_dens_data_i['O2_p1_number_density'][ds_key_lr]
        num_dens_y_axis_values_hr = num_dens_data_i_hr['O2_p1_number_density'][ds_key_hr]

        flux_x_vals_lr = vx_y_axis_values * num_dens_y_axis_values * 1e5
        flux_x_vals_hr = vx_y_axis_values_hr * num_dens_y_axis_values_hr * 1e5

        cv_a = (np.abs(y[0]) + y[-1]) * (np.abs(z[0]) + z[-1])
        covered_area = cv_a * (3390.0**2) * 1e10 #covered area of slice in cm^2
        escape_rate_lr = np.mean(flux_x_vals_lr) * covered_area
        escape_rate_hr = np.mean(flux_x_vals_hr) * covered_area

        f_lr = open(lr_tail_txt_name,'ab')
        dataout_lr = np.column_stack((np.mean(flux_x_vals_lr), escape_rate_lr, cv_a))
        np.savetxt(f_lr, dataout_lr, delimiter=',',fmt='%10.2e', comments='#')
        f_lr.close()

        f_hr = open(hr_tail_txt_name,'ab')
        dataout_hr = np.column_stack((np.mean(flux_x_vals_hr), escape_rate_hr, cv_a))
        np.savetxt(f_hr, dataout_hr, delimiter=',',fmt='%10.2e', comments='#')
        f_hr.close()
