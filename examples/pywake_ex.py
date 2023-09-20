import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import os
import yaml
from py_wake.site import XRSite
from py_wake.wind_turbines import WindTurbine
from py_wake.wind_turbines.power_ct_functions import PowerCtTabular
from py_wake.examples.data.hornsrev1 import Hornsrev1Site

# constructor for YAML !include command
# (this is included in windio.utils)
def include_constructor(loader, node):
    filepath = loader.construct_scalar(node)
    base_dir = os.path.dirname(loader.stream.name)
    abs_filepath = os.path.join(base_dir, filepath)
    
    with open(abs_filepath, 'r') as f:
        return yaml.safe_load(f)

def includeTimeseriesNetCDF(loader, node):
    filepath = loader.construct_scalar(node)
    base_dir = os.path.dirname(loader.stream.name)
    abs_filepath = os.path.join(base_dir, filepath)
    
    timeseries = xr.open_dataset(abs_filepath)
    timeseries_dicts = [{**{'time': str(time)}, **{var: float(data_vars[var].values) for var in data_vars.keys()}}
                    for time, data_vars in timeseries.groupby('time')]
    return timeseries_dicts


def includeBathymetryNetCDF(loader, node):
    filepath = loader.construct_scalar(node)
    base_dir = os.path.dirname(loader.stream.name)
    abs_filepath = os.path.join(base_dir, filepath)

    dataset = xr.open_dataset(abs_filepath)
    bathymetry_data = {variable: list(dataset[variable].values.flatten()) for variable in dataset.variables}
    return bathymetry_data


yaml.SafeLoader.add_constructor('!includeBathymetryNetCDF', includeBathymetryNetCDF)
yaml.SafeLoader.add_constructor('!includeTimeseriesNetCDF', includeTimeseriesNetCDF)
yaml.SafeLoader.add_constructor('!include', include_constructor)


system = sys.argv[1]
#system = 'examples/plant/wind_energy_system/IEA37_case_study_3_wind_energy_system.yaml'
with open(system, "r") as stream:
    try:
        system_dat = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

#farm = 'examples/plant/plant_wind_farm/IEA37_case_study_3_wind_farm.yaml'
#with open(farm, "r") as stream:
#    try:
#        farm_dat = yaml.safe_load(stream)
#    except yaml.YAMLError as exc:
#        print(exc)
farm_dat = system_dat['wind_farm']

#resource = 'examples/plant/plant_energy_resource/IEA37_case_study_4_energy_resource.yaml'
#with open(resource, "r") as stream:
#    try:
#        resource_dat = yaml.safe_load(stream)
#    except yaml.YAMLError as exc:
#        print(exc)
# construct site
resource_dat = system_dat['site']['energy_resource']
if 'timeseries' in resource_dat['wind_resource'].keys():
   timeseries = True
   wind_resource_timeseries = resource_dat['wind_resource']['timeseries']
   times = [d['time'] for d in wind_resource_timeseries]
   ws = [d['speed'] for d in wind_resource_timeseries]
   wd = [d['direction'] for d in wind_resource_timeseries]
   assert(len(times) == len(ws))
   assert(len(wd) == len(ws))
   site = Hornsrev1Site()
   TI = None
   #ite = XRSite(xr.Dataset(
   # data_vars={'P': (('time'), np.ones(len(ws)) / len(speeds)), },
   # coords={'time': range(len(times)),
   #         'ws': speeds,
   #         'wd': wd}))

elif 'weibull_k' in resource_dat['wind_resource'].keys():
   A = resource_dat['wind_resource']['weibull_a']
   k = resource_dat['wind_resource']['weibull_k']
   freq = resource_dat['wind_resource']['sector_probability']
   wd = resource_dat['wind_resource']['wind_direction']
   ws = resource_dat['wind_resource']['wind_speed']
   site = XRSite(
          ds=xr.Dataset(data_vars=
                           {'Sector_frequency': ('wd', freq['data']), 
                            'Weibull_A': ('wd', A['data']), 
                            'Weibull_k': ('wd', k['data']), 
                            'TI': (resource_dat['wind_resource']['turbulence_intensity']['dims'][0], resource_dat['wind_resource']['turbulence_intensity']['data'])
                            },
                         coords={'wd': wd, 'ws': ws}))
   
   timeseries = False
   TI =  resource_dat['wind_resource']['turbulence_intensity']['data']
else:
   timeseries = False
   ws = resource_dat['wind_resource']['wind_speed']
   wd = resource_dat['wind_resource']['wind_direction']
   P = np.array(resource_dat['wind_resource']['probability']['data'])
   site = XRSite(ds=xr.Dataset(data_vars={'P': (['wd', 'ws'], P)}, coords = {'ws': ws, 'wd': wd, 'TI': resource_dat['wind_resource']['turbulence_intensity']['data']}))
   TI = resource_dat['wind_resource']['turbulence_intensity']['data']

# get x and y positions
x = farm_dat['layouts']['initial_layout']['coordinates']['x']
y = farm_dat['layouts']['initial_layout']['coordinates']['y']

# define turbine
hh = farm_dat['turbines']['hub_height']
rd = farm_dat['turbines']['rotor_diameter']
cp = farm_dat['turbines']['performance']['Cp_curve']['Cp_values']
cp_ws = farm_dat['turbines']['performance']['Cp_curve']['Cp_wind_speeds']
ct = farm_dat['turbines']['performance']['Ct_curve']['Ct_values']
ct_ws = farm_dat['turbines']['performance']['Ct_curve']['Ct_wind_speeds']
int_speeds = np.linspace(np.min(np.min([cp_ws, ct_ws])), np.max(np.max([cp_ws, ct_ws])), 10000)
cps_int = np.interp(int_speeds, cp_ws, cp)
cts_int = np.interp(int_speeds, ct_ws, ct)
windTurbines = WindTurbine(name=farm_dat['turbines']['name'], diameter=rd, hub_height=hh, 
                      powerCtFunction=PowerCtTabular(int_speeds, 0.5 * cps_int * int_speeds ** 3 * 1.225 * (rd / 2) ** 2 * np.pi, power_unit='W', ct=cts_int))
from py_wake import NOJ, BastankhahGaussian
#noj = BastankhahGaussian(site, turbine, turbulenceModel=None)
site.interp_method = 'linear'
from py_wake.rotor_avg_models import RotorCenter
noj = NOJ(site, windTurbines, turbulenceModel=None, k=0.05, rotorAvgModel=RotorCenter())
sim_res = noj(x, y, time=timeseries, ws=ws, wd=np.arange(0, 360, 1), TI=TI)
aep = sim_res.aep(normalize_probabilities=False).sum()
#aep = sim_res.aep(normalize_probabilities=not timeseries).sum()
print('aep is ', aep, 'GWh')
#print('aep is ', sim_res.aep().sum(), 'GWh')
print('(%.2f capcacity factor)' % ( aep / (len(x) * windTurbines.power(10000) * 8760 / 1e9)))
