# Example function to evaluate the AEP of WindIO reference wind farms using
# the software FLORIS
# 
# Assumptions:
# - The Jensen wake deficit model is chosen, therefore TI has no influence
# - The wind rose is provided @hub-height and only 1 turbine grid point is
#   considered in the solver, therefore shear and veer are negligible
#
# Note: Differences in AEP predicition compared to PyWake result from different
#   approaches for the interpolation of wind speed bins (from 30° to finer resultion)
# 
# Created: 2023-08-04 by Samuel Kainz
# Last modified: 2023-09-12 by Samuel Kainz

# --------------------
# INPUT
# Specify wind speed step width in [m/s] for AEP evaluation
ws_sw = 1
# Specify wind direction step width in [deg] for AEP evaluation
wd_sw = 1

# --------------------
# CODE
import numpy as np
import os
import inspect
import yaml
import sys
from floris.tools import FlorisInterface, WindRose
import xarray as xr

# Function to write floris turbine file (yaml) if it does not yet exist
def WriteTurbine(turbine_name, data):
    # Read parameters
    HH = data['wind_farm']['turbines']['hub_height']
    RD = data['wind_farm']['turbines']['rotor_diameter']
    cp = data['wind_farm']['turbines']['performance']['Cp_curve']['Cp_values']
    cp_ws = data['wind_farm']['turbines']['performance']['Cp_curve']['Cp_wind_speeds']
    ct = data['wind_farm']['turbines']['performance']['Ct_curve']['Ct_values']
    ct_ws = data['wind_farm']['turbines']['performance']['Ct_curve']['Ct_wind_speeds']
    # Interpolate for one wind speed vector if they are not equal
    # note: could be removed with the updated turbine file
    if not cp_ws == ct_ws:
        int_speeds = np.linspace(np.min(np.min([cp_ws + ct_ws])), np.max(np.max([cp_ws + ct_ws])), 10000)
        cps_int = np.interp(int_speeds, cp_ws, cp)
        cts_int = np.interp(int_speeds, ct_ws, ct)
        # convert to list
        cps_int = cps_int.tolist()
        cts_int = cts_int.tolist()
        int_speeds = int_speeds.tolist()
    else:
        int_speeds = cp_ws
        cps_int = cp
        cts_int = ct
    # Dummy values for Floris Cp / Ct curves (necessary for interpolation)
    int_speeds[0:0] = [0,int_speeds[0]-0.00001]
    cps_int[0:0] = [0,0]
    cts_int[0:0] = [0,0]
    int_speeds.extend([int_speeds[-1]+0.00001, 100])
    cps_int.extend([0,0])
    cts_int.extend([0,0])
    # Create dict with input values
    dict_file = {'turbine_type' : turbine_name, 'generator_efficiency' : 1.0,
                 'hub_height' : HH, 'pP' : 1.88, 'pT': 1.88,
                 'rotor_diameter' : RD, 'TSR' : 7.0, 'ref_density_cp_ct': 1.225,
                 'power_thrust_table': {'power' : cps_int, 'thrust' : cts_int, 'wind_speed' : int_speeds}
                 }
    # Add turbine name as single-quoted (work around to include single quotation in yaml writer)
    class SingleQuoted(str):
        pass
    def single_quoted_presenter(dumper, data):
        return dumper.represent_scalar('tag:yaml.org,2002:str', data, style="'")
    yaml.add_representer(SingleQuoted, single_quoted_presenter)
    dict_file.update(turbine_type = SingleQuoted(turbine_name))
    # Write turbine yaml-file
    with open(tur_path, 'w') as file:
        yaml.dump(dict_file, file, default_flow_style=False,sort_keys=False)
        
# constructor for YAML "!include" command
# (this is included in windio.utils)
def include_constructor(loader, node):
    filepath = loader.construct_scalar(node)
    base_dir = os.path.dirname(loader.stream.name)
    abs_filepath = os.path.join(base_dir, filepath)
    
    with open(abs_filepath, 'r') as f:
        return yaml.safe_load(f)

def includeBathymetryNetCDF(loader, node):
    filepath = loader.construct_scalar(node)
    base_dir = os.path.dirname(loader.stream.name)
    abs_filepath = os.path.join(base_dir, filepath)
    dataset = xr.open_dataset(abs_filepath)
    bathymetry_data = {variable: list(dataset[variable].values.flatten()) for variable in dataset.variables}
    return bathymetry_data

# add new constructor
yaml.SafeLoader.add_constructor('!include', include_constructor)
yaml.SafeLoader.add_constructor('!includeBathymetryNetCDF', includeBathymetryNetCDF)

# Open, load and store the required WindIO file as defined above
data =  {}
file_path = sys.argv[1]
with open(file_path) as f:
     data = yaml.load(f, Loader=yaml.SafeLoader)

# Discretize Weibull
# Extract required parameters
wb_scale = data['site']['energy_resource']['wind_resource']['weibull_a']['data']
wb_shape = data['site']['energy_resource']['wind_resource']['weibull_k']['data']
wb_wd_freq = data['site']['energy_resource']['wind_resource']['sector_probability']['data']
wb_wd = data['site']['energy_resource']['wind_resource']['wind_direction']
# Format and reshape
wb_scale = np.reshape(np.array(wb_scale),(-1,1))
wb_shape = np.reshape(np.array(wb_shape),(-1,1))
wb_wd_freq = np.reshape(np.array(wb_wd_freq),(-1,1))
# Wind speed discretization
step = 1
wb_ws = np.arange(0, 51, step)
# Upper and lower boundaries of wind speed bins
ws_low = np.arange(np.min(wb_ws)-step/2,np.max(wb_ws)+step/2,step)
ws_high = ws_low + step
ws_low[ws_low<0] = 0
ws_high[ws_high<0] = 0
# Discretize distribution for each wind direction and store in list (Weibull CDF)
freq_grid_raw = wb_wd_freq * ((1 - np.exp(-(1 / wb_scale * ws_high) ** wb_shape)) -
              (1 - np.exp(-(1 / wb_scale * ws_low) ** wb_shape)))
    
# Extract farm information
x = data['wind_farm']['layouts']['initial_layout']['coordinates']['x']
y = data['wind_farm']['layouts']['initial_layout']['coordinates']['y']
turbine_name = data['wind_farm']['turbines']['name']
cut_in = data['wind_farm']['turbines']['performance']['cutin_wind_speed']
cut_out = data['wind_farm']['turbines']['performance']['cutout_wind_speed']
turbine_power = data['wind_farm']['turbines']['performance']['rated_power']

# Check if turbine already in library, otherwise write yaml file
path = os.path.dirname(inspect.getfile(FlorisInterface))
path = path.replace('tools','turbine_library\\')
tur_path = path + turbine_name + '.yaml'
if not os.path.isfile(tur_path):
    WriteTurbine(turbine_name,data)
else:
    print('Note: turbine name already exists in library. No changes applied.')

# Initialize floris and update values
fi = FlorisInterface("floris_ex_input.yaml")
fi.reinitialize(layout_x = x, layout_y = y, turbine_type=[turbine_name])
fi.assign_hub_height_to_ref_height()

# Initialize wind rose, feed it, and interpolate
wind_rose = WindRose()
wd_grid_raw, ws_grid_raw = np.meshgrid(wb_wd,wb_ws,indexing="ij")
wind_rose.make_wind_rose_from_user_dist(
    np.array([j for i in np.reshape(wd_grid_raw,(-1,1)).tolist() for j in i]),
     np.array([j for i in np.reshape(ws_grid_raw,(-1,1)).tolist() for j in i]),
     np.array([j for i in np.reshape(freq_grid_raw,(-1,1)).tolist() for j in i]),
     wd=np.array(wb_wd),
     ws=np.array(wb_ws))
wd_int = np.arange(0, 360, wd_sw)
ws_int = np.arange(0, 51, ws_sw)
wd_grid, ws_grid = np.meshgrid(wd_int,ws_int,indexing="ij")
freq_grid = wind_rose.interpolate(wd_grid,ws_grid)
# make sure sum is 1
freq_grid = freq_grid / np.sum(freq_grid)

# Update Floris and evaluate AEP
fi.reinitialize(wind_directions=wd_int,wind_speeds=ws_int)
aep, t_power, f_power = fi.get_farm_AEP_and_power(freq = freq_grid, cut_in_wind_speed = cut_in, cut_out_wind_speed = cut_out+1)
print('aep is %.2f GWh' % (aep/1e9))

## Uncomment to quantify wake losses:
# aep_nowake, t_power_nowake, f_power_nowake = fi.get_farm_AEP_and_power(freq = freq_grid, cut_in_wind_speed = cut_in, cut_out_wind_speed = cut_out+1, no_wake=True)
# print('(%.2f capcacity factor)' % ( aep / (turbine_power * 8760 * len(x))))
# print('(%.2f%% wake losses)' % (100 - aep / aep_nowake * 100))