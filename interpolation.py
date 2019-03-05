import numpy as np
import h5py
import matplotlib.pyplot as plt
import getopt
from general_functions import *
from sliceplot import *
from patch_sampling import *
from sphereslice import *
from spherical_flux import *
from patch_sampling import *
import sys
plt.style.use(['seaborn-poster'])

ds_names_all, ds_types_all = get_datasets('R2349', False)
ds_keys = ['batsrus_mf_10km_3deg_95'] #these lines establish the path to the dataset
ds_names = {dsk:dsn for dsk, dsn in ds_names_all.items() if dsk in ds_keys}
ds_types = {dsk:[dsk] for dsk in ds_keys}

fields =['O2_p1_velocity_x','O2_p1_velocity_y','O2_p1_velocity_z','O2_p1_velocity_total','O2_p1_number_density',
         'O_p1_velocity_x', 'O_p1_velocity_y', 'O_p1_velocity_z', 'O_p1_velocity_total', 'O_p1_number_density']

f = h5py.File('batsrus_3d_multi_fluid_95000_10k_3deg.h5','r')

dat_x_values = f['x'][:]
dat_y_values = f['y'][:]
dat_z_values = f['z'][:]

def conversion_sphere_to_cart(r,theta,phi):
    '''converts spherical coords to x,y,z
    INPUTS:
    r: in Rm
    theta: (rad)
    phi: (rad)'''

    x = r * np.cos(theta) * np.sin(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(phi)
    return x, y, z

def conversion_cart_to_sphere(x,y,z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(y,x)
    phi = np.arccos(z/r)
    return r,theta,phi
    #return '{0:05f} , {1:05f} , {2:05f}'.format(r,theta,phi)
dat_r_values, dat_theta_values, dat_phi_values = conversion_cart_to_sphere(dat_x_values, dat_y_values, dat_z_values) #create spherical dat dat_r_values
###CALCULATING NUMBER OF SPHERES AND LOCATIONS###
r_sphere_100km = (3390.0+100.0) / 3390.0 #100km above surface
one_km = (3390.0 + 1.0) / 3390.0 #5km in Rm
two_km = (3390.0 + 2.0) / 3390.0 #5km in Rm
three_km = (3390.0 + 3.0) / 3390.0 #5km in Rm
four_km = (3390.0 + 4.0) / 3390.0 #5km in Rm
five_km = (3390.0 + 5.0) / 3390.0 #5km in Rm

circumfrence = 2*np.pi*r_sphere_100km #circumfrence in martian radii
#print circumfrence
a = 3490 #km #100km above surface
b = 3490 #km
c = 2 #km #diameter
h = np.sqrt((a**2) - 0.25*(c**2))
corner_angle = np.arctan(2*h / c) #corner angle of triangle in rad
corner_angle_deg = corner_angle * 180/np.pi #corner angle of triangle in deg
central_angle = 180 - 2*corner_angle_deg #central angle in degrees
arc_length = 2*np.pi*r_sphere_100km* (central_angle / 360) #arc length between spheres (Rm)
#print arc_length
number_of_spheres = circumfrence / arc_length #number of 5km sphere along 100km line
print "computing {} spheres".format(np.floor(number_of_spheres))

r_sphere_100km_array = np.repeat(r_sphere_100km, np.floor(number_of_spheres)) #create array of r values to match theta and phi
theta_values = np.linspace(0.0, 2*np.pi, np.floor(number_of_spheres)) #theta values in rad, 2192 is for 100km
phi_values = np.linspace(0.0,np.pi, np.floor(number_of_spheres))#phi values in rad,90N to -90S

reshape_number_of_spheres_factor = np.floor(number_of_spheres)
reshape_number_of_spheres_factor = reshape_number_of_spheres_factor.astype(int)
interpolation_spheres_coords_cartesian = np.array([]) #create an array of sphere locations to interpolate within, in cartesian coords
for i in range(len(theta_values)):
    interpolation_spheres_coords_cartesian = np.append(interpolation_spheres_coords_cartesian,conversion_sphere_to_cart(r_sphere_100km_array[i], theta_values[i], phi_values[i]))
interpolation_spheres_coords_cartesian = interpolation_spheres_coords_cartesian.reshape(reshape_number_of_spheres_factor,3)
print 'there are {} sets of coords for interpolation spheres'.format(len(interpolation_spheres_coords_cartesian))

###FINDING POINTS WITHIN SPHERES###
def check_sphere_for_points(x,y,z,cx,cy,cz,r):
    '''x,y,z: point coords
       cx,cy,cz: center of sphere
       r: radius of sphere (Rm)'''

    goodpoints = np.array([]) #x,y,z, distance to center coords for points within the sphere
    distance = (x - cx)**2 + (y - cy)**2 + (z - cz)**2 #distance from center of sphere
    if distance <= r**2: #condition for point to be inside of sphere
        goodpoints = np.append(goodpoints, np.array([x,y,z, distance]))#append points within sphere to the 'goodpoints' array

    return goodpoints

print 'computing points in all spheres...could take a while'
Points_in_spheres = np.array([])
#len_dat_x_values = len(dat_x_values)
#for i in len_dat_x_values:
for i in range(len(dat_x_values)): #loop through all data and return coords within the sphere
    x = dat_x_values[i]
    y = dat_y_values[i]
    z = dat_z_values[i]
    for j in range(reshape_number_of_spheres_factor):
        Points_in_spheres = np.append(Points_in_spheres,
                                  check_sphere_for_points(x, y, z,
                                                          interpolation_spheres_coords_cartesian[j,0], #x coord of each sphere
                                                          interpolation_spheres_coords_cartesian[j,1], #y coord of each sphere
                                                          interpolation_spheres_coords_cartesian[j,2], #z coord of each sphere
                                                          r=five_km))
Points_in_spheres = Points_in_spheres.reshape(len(Points_in_spheres)/4, 4)
print 'found {} points in all spheres'.format(len(Points_in_spheres))

###PULLING OUT FIELD VALUES AT POINTS WITHIN SPHERE
print 'finding indices of all points in spheres'
all_points_in_spheres_indices = np.array([]) # array for indices
for i in range(len(Points_in_sphere)):
    if (dat_x_values[i] == dat_x_values[i-1]) & (dat_y_values[i] == dat_y_values[i-1]) & (dat_z_values[i] == dat_z_values[i-1]):
            pass #skip all duplicate points at end of .dat file
    else:
        all_points_in_spheres_indices = np.append(all_points_in_spheres_indices,
                                                 np.argwhere((dat_x_values==Points_in_spheres[i][0]) &
                                                             (dat_y_values==Points_in_spheres[i][1]) &
                                                             (dat_z_values==Points_in_spheres[i][2])))
all_points_in_spheres_indices = all_points_in_spheres_indices.astype(int)
        #appending indices of .dat file values that match the x,y,z coords found in the sphere
print 'finding and removing duplicates'
Duplicate_indices = np.array([])
for i in range(len(all_points_in_spheres_indices)):
    if (dat_x_values[all_points_in_spheres_indices[i]] == dat_x_values[all_points_in_spheres_indices[i-1]]) & (dat_y_values[all_points_in_spheres_indices[i]] == dat_y_values[all_points_in_spheres_indices[i-1]]) & (dat_z_values[all_points_in_spheres_indices[i]] == dat_z_values[all_points_in_spheres_indices[i-1]]):
        Duplicate_indices = np.append(Duplicate_indices, all_points_in_spheres_indices[i])
        #Duplicate_indices = np.append(Duplicate_indices, Points_in_one_sphere_indices[i])
Duplicate_indices = Duplicate_indices.astype(int)

List_points_in_spheres_indices = all_points_in_spheres_indices.tolist() #convert all indices to list
List_Duplicate_indices = Duplicate_indices.tolist() #convert duplicate indices to list
for l in range(len(List_Duplicate_indices)): #loop through indices, removing duplicates
    List_points_in_spheres_indices.remove(List_Duplicate_indices[l])
Final_points_in_spheres_indices = np.array(List_points_in_spheres_indices) #convert good indices back to array
if len(Points_in_spheres) != len(Final_points_in_spheres_indices):
    sys.exit("Points arrays dont match in length")

###INTERPOLATING###
print 'interpolating and finding new values'
dat_O2_p1_num_dens_values = f['O2_p1_number_density'][:][Final_points_in_spheres_indices] #Pulling out values at new indices within interp sphere
dat_O2_p1_v_x = f['O2_p1_velocity_x'][:][Final_points_in_spheres_indices]
dat_O2_p1_v_y = f['O2_p1_velocity_y'][:][Final_points_in_spheres_indices]
dat_O2_p1_v_z = f['O2_p1_velocity_z'][:][Final_points_in_spheres_indices]
dat_O2_p1_vtot = np.sqrt(dat_O2_p1_vx**2 + dat_O2_p1_vy**2 + dat_O2_p1_vz**2)

dat_O_p1_num_dens_values = f['O_p1_number_density'][:][Final_points_in_spheres_indices]
dat_O_p1_v_x = f['O_p1_velocity_x'][:][Final_points_in_spheres_indices]
dat_O_p1_v_y = f['O_p1_velocity_y'][:][Final_points_in_spheres_indices]
dat_O_p1_v_z = f['O_p1_velocity_z'][:][Final_points_in_spheres_indices]
dat_O_p1_vtot = np.sqrt(dat_O_p1_vx**2 + dat_O_p1_vy**2 + dat_O_p1_vz**2)

interpolation_fields = np.array([dat_O2_p1_num_dens_values, dat_O2_p1_vx, dat_O2_p1_vy, dat_O2_p1_vz, dat_O2_p1_vtot,
                             dat_O_p1_num_dens_values, dat_O_p1_vx, dat_O_p1_vy, dat_O_p1_vz, dat_O_p1_vtot])

distances = Points_in_spheres[:,3] #array of distances from center of sphere to points
weights = 1.0 / distances #1/di factor

interpolated_values_array = np.array([])
for n in range(len(interpolation_fields)):
    numerator = np.sum(weights*interpolation_fields[n])
    denominator = np.sum(weights)
    interpolated_values_array = np.append(interpolated_values_array, numerator / denominator)
combine_coords_value = np.concatenate((interpolation_spheres_coords_cartesian,interpolated_values_array),axis=0)
reshaped_combined_values = combine_coords_value.reshape(1,len(combine_coords_value))
###WRITING VALUES TO A FILE###
#with open("test_interpolation_file.txt") as file:
#    file.write(interpolated_value)
#    file.close()
file=open('test_interpolation_file.txt','ab') #opens file to append without overwriting
dataout = interpolation_spheres_coords_cartesian#np.column_stack((interpolation_spheres_coords_cartesian, interpolated_value))
np.savetxt(file, dataout, delimiter=',',header='#', fmt='%-10.10e',comments='#') #saves columnated values as sphere center coords, interpolated value
file.close()
