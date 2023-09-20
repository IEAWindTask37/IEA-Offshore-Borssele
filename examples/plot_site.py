import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from utils.yml_utils import validate_yaml, Loader, load_yaml
import os
import yaml
import xarray as xr


# ensure netcdf reader is present
def includeBathymetryNetCDF(self, node):
    filename = os.path.join(self._root,
                            self.construct_scalar(node))
    dataset = xr.open_dataset(filename)
    bathymetry_data = {variable:
                       list(dataset[variable].values.reshape(-1))
                       for variable in dataset.variables}
    return bathymetry_data


Loader.includeBathymetryNetCDF = includeBathymetryNetCDF
Loader.add_constructor('!includeBathymetryNetCDF',
                       Loader.includeBathymetryNetCDF)


# load wind farm systems
#with open('../IEA37_Borssele_Optimized_System.yaml', 'r') as f:
with open('../IEA37_Borssele_Regular_System.yaml', 'r') as f:
    regular_system = yaml.load(f, Loader)

with open('../IEA37_Borssele_irregular_System.yaml', 'r') as f:
    irregular_system = yaml.load(f, Loader)

# extract site and wind farm
b = regular_system['site']
regular = regular_system['wind_farm']
irregular = irregular_system['wind_farm']

# extract bathymetry data
# dim(x)!=dim(y)!=dim(z), array must be properly processed
X = np.array(regular_system['site']['Bathymetry']['x'])
Y = np.array(regular_system['site']['Bathymetry']['y'])
Z = np.array(regular_system['site']['Bathymetry']['depth'])
xi, yi = np.meshgrid(X, Y)
zi = griddata((xi.reshape(-1), np.flip(yi).reshape(-1)),
              Z, (xi, yi), method='linear')

# create plot
fig, ax = plt.subplots(1, 2, figsize=(6, 4), sharey=True)

# filled contour plots
zi[np.isclose(zi, 70)] = np.nan
CS = ax[0].contourf(xi, yi, zi, 300, cmap=plt.cm.binary,)
ax[1].contourf(xi, yi, zi, 300, cmap=plt.cm.binary,)
cb_ax = fig.add_axes([0.93, 0.2, 0.02, 0.6])
cb = fig.colorbar(CS, cax=cb_ax)
cb.set_label('Depth (m)')
cb.ax.invert_yaxis()
cb.set_ticks([20, 30, 40])

# plot boundaries
regx = regular['layouts']['initial_layout']['coordinates']['x']
regy = regular['layouts']['initial_layout']['coordinates']['y']
irrgx = irregular['layouts']['initial_layout']['coordinates']['x']
irrgy = irregular['layouts']['initial_layout']['coordinates']['y']
for ii in range(2):
    for jj in range(len(b['boundaries']['polygons'][0]['x']) - 1):
        ax[ii].plot(b['boundaries']['polygons'][0]['x'][jj:jj+2],
                    b['boundaries']['polygons'][1]['y'][jj:jj+2],
                    color='firebrick', ls='--')
    ax[ii].plot(np.array(b['boundaries']['polygons'][0]['x'])[[-1, 0]],
                np.array(b['boundaries']['polygons'][1]['y'])[[-1, 0]],
                color='firebrick', label='Site Boundaries', ls='--')
ax[1].scatter(irrgx, irrgy, c='darkorange', marker='2', zorder=3)
ax[0].scatter(regx, regy, c='darkorange', marker='2',
              zorder=3, label='Turbine Locations')

# decorate plot
ax[0].legend(prop={'size': 8})
ax[0].set_ylabel('Northing (m)')
for ii in range(2):
    ax[ii].set_aspect(1)
    ax[ii].set_xlabel('Easting (m)')

#ax[0].set_title('Optimized Layout, 3437.44 GWh', pad=10)
ax[0].set_title('Regular Layout\n3386.98 GWh', pad=10)
ax[1].set_title('Irregular Layout\n3437.19 GWh', pad=10)
plt.savefig('layouts.pdf', bbox_inches='tight')
plt.clf()


# plot turbine cp/ct curves
fig, ax = plt.subplots()
cp = regular['turbines']['performance']['Cp_curve']
ct = regular['turbines']['performance']['Ct_curve']
ax.plot(cp['Cp_wind_speeds'], cp['Cp_values'], label='$c_p$', c='k')
ax2 = ax.twinx()
ax2.plot(ct['Ct_wind_speeds'], ct['Ct_values'], label='$c_t$', c='k', ls='--')
ax.legend(loc='upper center')
ax2.legend(loc='upper right')
ax.set_ylabel('Power Coefficient')
ax2.set_ylabel('Thrust Coefficient')
ax.set_xlabel('Wind Speed (m/s)')
plt.savefig('Turbine.pdf')
plt.clf()
