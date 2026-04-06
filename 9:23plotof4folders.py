# updated on 09/23
# Concatenate all years
# plot 4 folders based on the time sequence (years) with a certain coordinate (lat, lon)
# inplement dictionary to make the script neat, the keys are names of sub_folders 
# and the values include time/pr data

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import glob
import os

# path
data_main_dir = '/Users/zhangyiyi/Desktop/datasets/model1'
sub_folders = ['historical', 'ssp245', 'ssp370', 'ssp585']

# define coordinate
defined_lat = 23
defined_lon = 113

# Function to load, concatenate, and extract precipitation data for a given sub-folder
def load_and_extract_pr(sub_folder):
    # access .nc
    file_paths = glob.glob(os.path.join(data_main_dir, sub_folder, '*.nc'))
    # initialize empty lists for the precipitation data and time data
    all_prdata = []
    all_timedata = []
    for file_path in file_paths:
        # open datasets
        dataset = nc.Dataset(file_path, mode='r')
        # and each variables inside
        precipitation = dataset.variables['pr'][:]  # pr
        time = dataset.variables['time'][:]  # time series
        lat = dataset.variables['lat'][:]  # lat
        lon = dataset.variables['lon'][:]  # lon
        # convert to datatime type
        time_units = dataset.variables['time'].units
        calendar = dataset.variables['time'].calendar
        time_dates = nc.num2date(time, units=time_units, calendar=calendar)
        
        # search the closest index of lat and lon
        lat_idx = np.abs(lat - defined_lat).argmin()
        lon_idx = np.abs(lon - defined_lon).argmin()
        # extract pr time series
        precip_time_series = precipitation[:, lat_idx, lon_idx]
        # convert units from kg/m^2/s to mm/day
        # 1 kg/m²/s is equivalent to 1 mm/s, 
        # because 1 mm of rain is defined as 1 liter of water (1 kg) falling on an area of 1 square meter.
        precip_time_series = precip_time_series * 86400 # There are 86400 seconds in a day
        print("Precipitation converted to mm/day")
        # deal with missing value, to NaN
        missing_value = dataset.variables['pr'].missing_value
        precipitation = np.where(precipitation == missing_value, np.nan, precipitation)
        # Append data
        all_prdata.extend(precip_time_series)
        all_timedata.extend(time_dates)
        
        # Close the dataset
        dataset.close()
        
    return np.array(all_timedata), np.array(all_prdata)

# Dictionary to store time series for each sub-folder (model/scenario)
precip_data_dict = {}

# Loop over each sub-folder (historical, ssp245, etc.) and load the data
for sub_folder in sub_folders:
    print(f"Loading data from {sub_folder}...")
    time_data, precip_data = load_and_extract_pr(sub_folder)
    precip_data_dict[sub_folder] = {'time': time_data, 'precipitation': precip_data}
    

# plot the time series fort each model/scenario
plt.figure(figsize=(20,9)) 
allyears=[]
for sub_folder in sub_folders:
    time_data = precip_data_dict[sub_folder]['time']
    precip_data = precip_data_dict[sub_folder]['precipitation']
    # convert time data to years
    years = np.array([t.year for t in time_data])
    # collect all years for the x-axis
    allyears.append(years)
    # check the min and max years 
    print(f"{sub_folder} years range:{min(years)} to {max(years)}")  
    # plot the time series for this model
    plt.plot(years, precip_data, label=sub_folder) 
# To show each year
allyears = np.concatenate(allyears)# concatenate all years into a single array after collecting them
xticks = np.arange(allyears.min(), allyears.max()+1, 5)
plt.gca().set_xticks(xticks)
plt.gca().set_xticklabels([str(x) for x in xticks], rotation = 45, ha = 'right')


# plot setting up
plt.title(f"Precipitation time series at (lat={defined_lat}, lon={defined_lon})")
plt.xlabel("Year")
plt.ylabel("Precipitation (mm/day)")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.15), ncol=len(sub_folders), fancybox=True)
plt.grid(True)

# Save the plot
plt.savefig('/Users/zhangyiyi/Desktop/datasets/Precipitation time series at Guangzhou.png')

# Show the plot
plt.tight_layout()
plt.show()

# Ensure all arrays have the same length
max_length=max(len(data['precipitation'])for data in precip_data_dict.values())
all_npydata=[data['precipitation']for data in precip_data_dict.values()]
all_npydata=[np.resize(arr, max_length) for arr in all_npydata]

# Convert the list of arrays to numpy array
all_data_array = np.concatenate(all_npydata)

# Save all .npy data in one file
np.save('/Users/zhangyiyi/Desktop/datasets/all precipitation npy data.npy', all_data_array)
print("Data has been saved")

# print all npy data
load_npy = np.load('/Users/zhangyiyi/Desktop/datasets/all precipitation npy data.npy', allow_pickle=True)
print(load_npy)