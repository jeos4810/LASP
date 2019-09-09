###This file generates two .txt files. One for plume, one for tail. Each with flux, escape rate, and covered plane area
import numpy as np
import h5py
import matplotlib.pyplot as plt
from general_functions import *
from spherical_flux import *
from sliceplot import *

def flux_escape_calculation_all_ts(lr_tail_txt_name, hr_tail_txt_name):
    ds_names_all, ds_types_all = get_datasets(load_key='R2349')
    ds_keys_lr = ['batsrus_mf_10km_3deg_50','batsrus_mf_10km_3deg_55','batsrus_mf_10km_3deg_60','batsrus_mf_10km_3deg_65',
                 'batsrus_mf_10km_3deg_70','batsrus_mf_10km_3deg_75','batsrus_mf_10km_3deg_80','batsrus_mf_10km_3deg_85',
                 'batsrus_mf_10km_3deg_90','batsrus_mf_10km_3deg_95']
    ds_keys_hr = ['batsrus_mf_5km_15deg_50','batsrus_mf_5km_15deg_55','batsrus_mf_5km_15deg_60','batsrus_mf_5km_15deg_65',
                 'batsrus_mf_5km_15deg_70','batsrus_mf_5km_15deg_75','batsrus_mf_5km_15deg_80','batsrus_mf_5km_15deg_85',
                 'batsrus_mf_5km_15deg_90','batsrus_mf_5km_15deg_95']

    fields_num_dens = ['O2_p1_number_density']
    fields_vx = ['O2_p1_velocity_x']
    #fields_vz = ['O2_p1_velocity_z']
    print 'Setting up tail txt files'
    setup_vals = []
    f_lr_setup = open(lr_tail_txt_name,'ab')
    f_lr_setup_vals = (setup_vals)
    np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',', header='plane, escape rate, covered area',fmt='%10.2e', comments='#')
    f_lr_setup.close()

    f_hr_setup = open(hr_tail_txt_name,'ab')
    f_hr_setup_vals = (setup_vals)
    np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='plane, escape rate, covered area',fmt='%10.2e', comments='#')
    f_hr_setup.close()

    x_planes = np.arange(-1.8,-3.0,-0.1)
    y = np.linspace(-2.0, 2.0, 40)
    z = np.linspace(-2.0, 2.0, 40)
    meshgrid_y, meshgrid_z = np.meshgrid(y,z)
    flat_meshgrid_y = meshgrid_y.flatten()
    flat_meshgrid_z = meshgrid_z.flatten()
    for i in range(len(ds_keys_lr)):
        f_lr_setup = open(lr_tail_txt_name,'ab')
        f_lr_setup_vals = (setup_vals)
        np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',',header='---timestep {}'.format(ds_keys_lr[i]) ,fmt='%10.2e', comments='#')
        f_lr_setup.close()

        f_hr_setup = open(hr_tail_txt_name,'ab')
        f_hr_setup_vals = (setup_vals)
        np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='---timestep {}'.format(ds_keys_hr[i]),fmt='%10.2e', comments='#')
        f_hr_setup.close()
        ds_key_lr = [ds_keys_lr[i]] #these lines establish the path to the dataset
        ds_key_hr = [ds_keys_hr[i]] # path to high res data
        ds_names = {dsk:dsn for dsk, dsn in ds_names_all.items() if dsk in ds_key_lr}
        ds_types = {dsk:[dsk] for dsk in ds_key_lr}
        ds_names_hr = {dsk:dsn for dsk, dsn in ds_names_all.items() if dsk in ds_key_hr}
        ds_types_hr = {dsk:[dsk] for dsk in ds_key_hr}
        print ds_key_lr
###Tail Calculation
        print 'calculating tail planes for {} and {}'.format(ds_keys_lr[i], ds_keys_hr[i])
        for n in range(len(x_planes)):
            print 'calculating tail plane {}'.format(x_planes[n])
            x = np.repeat(x_planes[n],1600)
            coords = np.array([x, flat_meshgrid_y, flat_meshgrid_z])

            indxs = get_path_idxs(coords, ds_names, ds_types)
            indxs_hr = get_path_idxs(coords, ds_names_hr, ds_types_hr)
            vx_data_i = get_all_data(ds_names, ds_types, indxs, fields_vx)
            vx_data_i_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields_vx)
            num_dens_data_i = get_all_data(ds_names, ds_types, indxs, fields_num_dens)
            num_dens_data_i_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields_num_dens)

            vx_y_axis_values = vx_data_i['O2_p1_velocity_x'][ds_keys_lr[i]]
            vx_y_axis_values_hr = vx_data_i_hr['O2_p1_velocity_x'][ds_keys_hr[i]]
            num_dens_y_axis_values = num_dens_data_i['O2_p1_number_density'][ds_keys_lr[i]]
            num_dens_y_axis_values_hr = num_dens_data_i_hr['O2_p1_number_density'][ds_keys_hr[i]]

            flux_x_vals_lr = vx_y_axis_values * num_dens_y_axis_values * 1e5
            flux_x_vals_hr = vx_y_axis_values_hr * num_dens_y_axis_values_hr * 1e5

            cv_a = (np.abs(y[0]) + y[-1]) * (np.abs(z[0]) + z[-1])
            covered_area = cv_a * (3390.0**2) * 1e10 #covered area of slice in cm^2
            escape_rate_lr = np.mean(flux_x_vals_lr) * covered_area
            escape_rate_hr = np.mean(flux_x_vals_hr) * covered_area

            f_lr = open(lr_tail_txt_name,'ab')
            dataout_lr = np.column_stack((x_planes[n], escape_rate_lr, cv_a))
            np.savetxt(f_lr, dataout_lr, delimiter=',',fmt='%10.2e', comments='#')
            f_lr.close()

            f_hr = open(hr_tail_txt_name,'ab')
            dataout_hr = np.column_stack((x_planes[n], escape_rate_hr, cv_a))
            np.savetxt(f_hr, dataout_hr, delimiter=',',fmt='%10.2e', comments='#')
            f_hr.close()

def flux_escape_calculation_all_ts_plume(lr_plume_txt_name, hr_plume_txt_name):
    ds_names_all, ds_types_all = get_datasets(load_key='R2349')
    ds_keys_lr = ['batsrus_mf_10km_3deg_50','batsrus_mf_10km_3deg_55','batsrus_mf_10km_3deg_60','batsrus_mf_10km_3deg_65',
                 'batsrus_mf_10km_3deg_70','batsrus_mf_10km_3deg_75','batsrus_mf_10km_3deg_80','batsrus_mf_10km_3deg_85',
                 'batsrus_mf_10km_3deg_90','batsrus_mf_10km_3deg_95']
    ds_keys_hr = ['batsrus_mf_5km_15deg_50','batsrus_mf_5km_15deg_55','batsrus_mf_5km_15deg_60','batsrus_mf_5km_15deg_65',
                 'batsrus_mf_5km_15deg_70','batsrus_mf_5km_15deg_75','batsrus_mf_5km_15deg_80','batsrus_mf_5km_15deg_85',
                 'batsrus_mf_5km_15deg_90','batsrus_mf_5km_15deg_95']

    fields_num_dens = ['O2_p1_number_density']
    #fields_vx = ['O2_p1_velocity_x']
    fields_vz = ['O2_p1_velocity_z']
    print 'Setting up plume txt files'
    setup_vals = []
    f_lr_setup = open(lr_plume_txt_name,'ab')
    f_lr_setup_vals = (setup_vals)
    np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',', header='plane, escape rate, covered area',fmt='%10.2e', comments='#')
    f_lr_setup.close()

    f_hr_setup = open(hr_plume_txt_name,'ab')
    f_hr_setup_vals = (setup_vals)
    np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='plane, escape rate, covered area',fmt='%10.2e', comments='#')
    f_hr_setup.close()

    z_planes = np.arange(1.8, 3.0, 0.1)
    x = np.linspace(-2.0, 2.0, 40)
    y = np.linspace(-2.0, 2.0, 40)
    meshgrid_x, meshgrid_y = np.meshgrid(x,y)
    flat_meshgrid_x = meshgrid_x.flatten()
    flat_meshgrid_y = meshgrid_y.flatten()
    for i in range(len(ds_keys_lr)):
        f_lr_setup = open(lr_plume_txt_name,'ab')
        f_lr_setup_vals = (setup_vals)
        np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',',header='---timestep {}'.format(ds_keys_lr[i]) ,fmt='%10.2e', comments='#')
        f_lr_setup.close()

        f_hr_setup = open(hr_plume_txt_name,'ab')
        f_hr_setup_vals = (setup_vals)
        np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='---timestep {}'.format(ds_keys_hr[i]),fmt='%10.2e', comments='#')
        f_hr_setup.close()
        ds_key_lr = [ds_keys_lr[i]] #these lines establish the path to the dataset
        ds_key_hr = [ds_keys_hr[i]] # path to high res data
        ds_names = {dsk:dsn for dsk, dsn in ds_names_all.items() if dsk in ds_key_lr}
        ds_types = {dsk:[dsk] for dsk in ds_key_lr}
        ds_names_hr = {dsk:dsn for dsk, dsn in ds_names_all.items() if dsk in ds_key_hr}
        ds_types_hr = {dsk:[dsk] for dsk in ds_key_hr}
        print ds_key_lr
###Tail Calculation
        print 'calculating plume planes for {} and {}'.format(ds_keys_lr[i], ds_keys_hr[i])
        for n in range(len(z_planes)):
            print 'calculating plume plane {}'.format(z_planes[n])
            z = np.repeat(z_planes[n],1600)
            coords = np.array([flat_meshgrid_x, flat_meshgrid_y, z])

            indxs = get_path_idxs(coords, ds_names, ds_types)
            indxs_hr = get_path_idxs(coords, ds_names_hr, ds_types_hr)
            vz_data_i = get_all_data(ds_names, ds_types, indxs, fields_vz)
            vz_data_i_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields_vz)
            num_dens_data_i = get_all_data(ds_names, ds_types, indxs, fields_num_dens)
            num_dens_data_i_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields_num_dens)

            vz_y_axis_values = vz_data_i['O2_p1_velocity_z'][ds_keys_lr[i]]
            vz_y_axis_values_hr = vz_data_i_hr['O2_p1_velocity_z'][ds_keys_hr[i]]
            num_dens_y_axis_values = num_dens_data_i['O2_p1_number_density'][ds_keys_lr[i]]
            num_dens_y_axis_values_hr = num_dens_data_i_hr['O2_p1_number_density'][ds_keys_hr[i]]

            flux_z_vals_lr = vz_y_axis_values * num_dens_y_axis_values * 1e5
            flux_z_vals_hr = vz_y_axis_values_hr * num_dens_y_axis_values_hr * 1e5

            cv_a = (np.abs(x[0]) + x[-1]) * (np.abs(y[0]) + y[-1])
            covered_area = cv_a * (3390.0**2) * 1e10 #covered area of slice in cm^2
            escape_rate_lr = np.mean(flux_z_vals_lr) * covered_area
            escape_rate_hr = np.mean(flux_z_vals_hr) * covered_area

            f_lr = open(lr_plume_txt_name,'ab')
            dataout_lr = np.column_stack((z_planes[n], escape_rate_lr, cv_a))
            np.savetxt(f_lr, dataout_lr, delimiter=',',fmt='%10.2e', comments='#')
            f_lr.close()

            f_hr = open(hr_plume_txt_name,'ab')
            dataout_hr = np.column_stack((z_planes[n], escape_rate_hr, cv_a))
            np.savetxt(f_hr, dataout_hr, delimiter=',',fmt='%10.2e', comments='#')
            f_hr.close()

def flux_escape_calculation_all_ts_spheres(lr_plume_txt_name, hr_plume_txt_name, lr_tail_txt_name, hr_tail_txt_name, lr_total_esc_txt_name, hr_total_esc_txt_name):
    def conversion_cart_to_sphere(x,y,z):
        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arctan2(y,x)
        phi = np.arccos(z/r)
        return r,theta,phi

    ds_names_all, ds_types_all = get_datasets(load_key='R2349')
    ds_keys_lr = ['batsrus_mf_10km_3deg_50','batsrus_mf_10km_3deg_55','batsrus_mf_10km_3deg_60','batsrus_mf_10km_3deg_65',
                 'batsrus_mf_10km_3deg_70','batsrus_mf_10km_3deg_75','batsrus_mf_10km_3deg_80','batsrus_mf_10km_3deg_85',
                 'batsrus_mf_10km_3deg_90','batsrus_mf_10km_3deg_95']
    ds_keys_hr = ['batsrus_mf_5km_15deg_50','batsrus_mf_5km_15deg_55','batsrus_mf_5km_15deg_60','batsrus_mf_5km_15deg_65',
                 'batsrus_mf_5km_15deg_70','batsrus_mf_5km_15deg_75','batsrus_mf_5km_15deg_80','batsrus_mf_5km_15deg_85',
                 'batsrus_mf_5km_15deg_90','batsrus_mf_5km_15deg_95']

    fields_flux = ['O2_p1_flux']

    print 'Setting up plume txt files'
    setup_vals = []
    f_lr_setup = open(lr_plume_txt_name,'ab')
    f_lr_setup_vals = (setup_vals)
    np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',', header='r_dist(Rm), escape rate, covered area',fmt='%10.2e',  comments='#')
    f_lr_setup.close()

    f_hr_setup = open(hr_plume_txt_name,'ab')
    f_hr_setup_vals = (setup_vals)
    np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='r_dist(Rm), escape rate, covered area',fmt='%10.2e', comments='#')
    f_hr_setup.close()

    print 'Setting up tail txt files'
    setup_vals = []
    f_lr_setup = open(lr_tail_txt_name,'ab')
    f_lr_setup_vals = (setup_vals)
    np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',', header='r_dist(Rm), escape rate, covered area',fmt='%10.2e', comments='#')
    f_lr_setup.close()

    f_hr_setup = open(hr_tail_txt_name,'ab')
    f_hr_setup_vals = (setup_vals)
    np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='r_dist(Rm), escape rate, covered area',fmt='%10.2e', comments='#')
    f_hr_setup.close()

    print 'Setting up total esc txt files'
    setup_vals = []
    f_lr_setup = open(lr_total_esc_txt_name,'ab')
    f_lr_setup_vals = (setup_vals)
    np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',', header='r_dist(Rm), escape rate, covered area',fmt='%10.2e', comments='#')
    f_lr_setup.close()

    f_hr_setup = open(hr_total_esc_txt_name,'ab')
    f_hr_setup_vals = (setup_vals)
    np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='r_dist(Rm), escape rate, covered area',fmt='%10.2e', comments='#')
    f_hr_setup.close()

    r_planes = np.arange(2.0,3.1,0.1)
    for i in range(len(ds_keys_lr)):
        f_lr_setup = open(lr_plume_txt_name,'ab')
        f_lr_setup_vals = (setup_vals)
        np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',',header='---timestep {}'.format(ds_keys_lr[i])     ,fmt='%10.2e', comments='#')
        f_lr_setup.close()

        f_hr_setup = open(hr_plume_txt_name,'ab')
        f_hr_setup_vals = (setup_vals)
        np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='---timestep  {}'.format(ds_keys_hr[i]),fmt='%10.2e', comments='#')
        f_hr_setup.close()

        f_lr_setup = open(lr_tail_txt_name,'ab')
        f_lr_setup_vals = (setup_vals)
        np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',',header='---timestep {}'.format(ds_keys_lr[i])     ,fmt='%10.2e', comments='#')
        f_lr_setup.close()

        f_hr_setup = open(hr_tail_txt_name,'ab')
        f_hr_setup_vals = (setup_vals)
        np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='---timestep  {}'.format(ds_keys_hr[i]),fmt='%10.2e', comments='#')
        f_hr_setup.close()

        f_lr_setup = open(lr_total_esc_txt_name,'ab')
        f_lr_setup_vals = (setup_vals)
        np.savetxt(f_lr_setup, f_lr_setup_vals, delimiter=',',header='---timestep {}'.format(ds_keys_lr[i])     ,fmt='%10.2e', comments='#')
        f_lr_setup.close()

        f_hr_setup = open(hr_total_esc_txt_name,'ab')
        f_hr_setup_vals = (setup_vals)
        np.savetxt(f_hr_setup, f_hr_setup_vals, delimiter=',', header='---timestep  {}'.format(ds_keys_hr[i]),fmt='%10.2e', comments='#')
        f_hr_setup.close()

        ds_key_lr = [ds_keys_lr[i]] #these lines establish the path to the dataset
        ds_key_hr = [ds_keys_hr[i]] # path to high res data
        ds_names = {dsk:dsn for dsk, dsn in ds_names_all.items() if dsk in ds_key_lr}
        ds_types = {dsk:[dsk] for dsk in ds_key_lr}
        ds_names_hr = {dsk:dsn for dsk, dsn in ds_names_all.items() if dsk in ds_key_hr}
        ds_types_hr = {dsk:[dsk] for dsk in ds_key_hr}
        #print ds_key_lr

        print 'calculating escape for {} and {}'.format(ds_keys_lr[i], ds_keys_hr[i])

        for n in range(len(r_planes)):
            xy, coords, rhat, area = create_sphere_mesh(r=r_planes[n],d_angle=5)
            r_vals, theta, phi = conversion_cart_to_sphere(coords[0],coords[1],coords[2])
            lon = (180.0/np.pi) * theta
            lat = 90.0 - ((180.0/np.pi) * phi)
            plume_lat_req = lat >= 45.0
            tail_lat_req = lat <= 44.0
            indxs_lr = get_path_idxs(coords, ds_names_lr, ds_types_lr)
            indxs_hr = get_path_idxs(coords, ds_names_hr, ds_types_hr)
            data_lr = get_all_data(ds_names_lr, ds_types_lr, indxs_lr, fields=fields_flux, normal=rhat)
            data_hr = get_all_data(ds_names_hr, ds_types_hr, indxs_hr, fields=fields_flux, normal=rhat)
            fdat_flux_lr = data_lr['O2_p1_flux'][ds_key_lr]
            fdat_flux_hr = data_hr['O2_p1_flux'][ds_key_hr]
            fdat_flux_lr = fdat_flux_lr.reshape(xy[0].shape)
            fdat_flux_hr = fdat_flux_hr.reshape(xy[0].shape)
            flat_fdat_flux_lr = fdat_flux_lr.flatten()
            flat_fdat_flux_hr = fdat_flux_hr.flatten()

            flat_fdat_plume_flux_lr = flat_fdat_flux_lr[plume_lat_req]
            flat_fdat_plume_flux_hr = flat_fdat_flux_hr[plume_lat_req]
            flat_fdat_tail_flux_lr = flat_fdat_flux_lr[tail_lat_req]
            flat_fdat_tail_flux_hr = flat_fdat_flux_hr[tail_lat_req]

            escape_plume_fdat_lr = flat_fdat_plume_flux_lr * area[plume_lat_req]
            escape_plume_fdat_hr = flat_fdat_plume_flux_hr * area[plume_lat_req]
            escape_tail_fdat_lr = flat_fdat_tail_flux_lr * area[tail_lat_req]
            escape_tail_fdat_hr = flat_fdat_tail_flux_hr * area[tail_lat_req]
            final_plume_escape_lr = np.sum(escape_plume_fdat_lr)
            final_plume_escape_hr = np.sum(escape_plume_fdat_hr)
            final_tail_escape_lr = np.sum(escape_tail_fdat_lr)
            final_tail_escape_hr = np.sum(escape_tail_fdat_hr)

            total_escape_lr = flat_fdat_flux_lr * area
            total_escape_hr = flat_fdat_flux_hr * area
            final_total_escape_lr = np.sum(total_escape_lr)
            final_total_escape_hr = np.sum(total_escape_hr)

            f_lr = open(lr_plume_txt_name,'ab')
            dataout_lr = np.column_stack((r_planes[n], final_plume_escape_lr, np.sum(area[plume_lat_req])))
            np.savetxt(f_lr, dataout_lr, delimiter=',',fmt='%10.2e', comments='#')
            f_lr.close()

            f_hr = open(hr_plume_txt_name,'ab')
            dataout_hr = np.column_stack((r_planes[n], final_plume_escape_hr, np.sum(area[plume_lat_req])))
            np.savetxt(f_hr, dataout_hr, delimiter=',',fmt='%10.2e', comments='#')
            f_hr.close()

            f_lr = open(lr_tail_txt_name,'ab')
            dataout_lr = np.column_stack((r_planes[n], final_tail_escape_lr, np.sum(area[tail_lat_req])))
            np.savetxt(f_lr, dataout_lr, delimiter=',',fmt='%10.2e', comments='#')
            f_lr.close()

            f_hr = open(hr_tail_txt_name,'ab')
            dataout_hr = np.column_stack((r_planes[n], final_tail_escape_hr, np.sum(area[tail_lat_req])))
            np.savetxt(f_hr, dataout_hr, delimiter=',',fmt='%10.2e', comments='#')
            f_hr.close()

            f_lr = open(lr_total_esc_txt_name,'ab')
            dataout_lr = np.column_stack((r_planes[n], final_total_escape_lr, np.sum(area)))
            np.savetxt(f_lr, dataout_lr, delimiter=',',fmt='%10.2e', comments='#')
            f_lr.close()

            f_hr = open(hr_total_esc_txt_name,'ab')
            dataout_hr = np.column_stack((r_planes[n], final_total_escape_hr, np.sum(area)))
            np.savetxt(f_hr, dataout_hr, delimiter=',',fmt='%10.2e', comments='#')
            f_hr.close()
