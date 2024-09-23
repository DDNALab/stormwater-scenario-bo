import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import glob
import os

# create a list of multi-models paths 
model_paths = [
    '/Users/zhangyiyi/Desktop/datasets/model1',
    '/Users/zhangyiyi/Desktop/datasets/model2' 
]

# list of each sub_folder
sub_folders = ['historical', 'ssp245', 'ssp370', 'ssp585']

# define coordinate
defined_lat = 23
defined_lon = 113

# Function to load, concatenate, and extract precipitation data for a given sub-folder
def load_and_extract_pr(model_dir, sub_folder):
    # access .nc
    file_paths = glob.glob(os.path.join(model_dir, sub_folder, '*.nc'))
    # initialize empty lists for the precipitation data and time data
    all_prdata = []
    all_timedata = []
    for file_path in file_paths:
        # open datasets
        dataset = nc.Dataset(file_path, mode='r')
        # and each variable inside
        precipitation = dataset.variables['pr'][:]
        time = dataset.variables['time'][:]
        lat = dataset.variables['lat'][:]
        lon = dataset.variables['lon'][:]
        # convert to datatime type
        time_units = dataset.variables['time'].units
        calendar = dataset.variables['time'].calendar
        time_dates = nc.num2date(time, units=time_units, calendar=calendar)
        
        # search the closest index of lat and lon
        lat_idx = np.abs(lat - defined_lat).argmin()
        lon_idx = np.abs(lon - defined_lon).argmin()
        # extract pr time series
        precip_time_series = precipitation[:, lat_idx, lon_idx] 
        
        # deal with NA
        missing_value = dataset.variables['pr'].missing_value
        precip_time_series = np.where(precip_time_series == missing_value, np.nan, precip_time_series)
        
        # convert unit from kg/m^2/s to mm/day
        # 1 kg/m²/s is equivalent to 1 mm/s, 
        # because 1 mm of rain is defined as 1 liter of water (1 kg) falling on an area of 1 square meter.
        precip_time_series = precip_time_series * 86400 # There are 86400 seconds in a day
        print("Precipitation converted to mm/day")
        # append data
        all_prdata.extend(precip_time_series)
        all_timedata.extend(time_dates)
        # close dataset
        dataset.close()
        
    return np.array(all_timedata), np.array(all_prdata)

# loop each model
for model_dir in model_paths:
    print(f"Processing model at {model_dir}...")
    # Dictionary to store time series for each sub-folder (model/scenario)
    precip_data_dict = {}
    
    # loop each sub_folder
    for sub_folder in sub_folders:
        print(f"Loading data from {sub_folder}...")
        time_data, precip_data = load_and_extract_pr(model_dir, sub_folder)
        precip_data_dict[sub_folder] = {'time': time_data, 'precipitation': precip_data}
    
    # plot the time series fort each model/scenario
    plt.figure(figsize=(20,9))
    allyears = []
    for sub_folder in sub_folders:
        time_data = precip_data_dict[sub_folder]['time']
        precip_data = precip_data_dict[sub_folder]['precipitation']
        # convert time data to years
        years = np.array([t.year for t in time_data])
        allyears.append(years)
        # check the min and max years 
        print(f"{sub_folder} years range:{min(years)} to {max(years)}")
        # plot the time series for this model  
        plt.plot(years, precip_data, label=f"{model_dir}/{sub_folder}")
    
    # collect data of each year
    allyears = np.concatenate(allyears)
    xticks = np.arange(allyears.min(), allyears.max()+1, 5)
    plt.gca().set_xticks(xticks)
    plt.gca().set_xticklabels([str(x) for x in xticks], rotation=45, ha='right')

    # plot setting up
    plt.title(f"Precipitation time series at (lat={defined_lat}, lon={defined_lon})")
    plt.xlabel("Year")
    plt.ylabel("Precipitation (mm/day)")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.15), ncol=len(sub_folders), fancybox=True)
    plt.grid(True)
    plt.tight_layout()
    #plt.show()
    
    # Save the plot
    plt.savefig(os.path.join(model_dir, 'Precipitation_time_series_at_Guangzhou.png'))
    plt.close()

    # save the data
    max_length = max(len(data['precipitation']) for data in precip_data_dict.values())
    all_npydata = [np.resize(arr, max_length) for arr in [data['precipitation'] for data in precip_data_dict.values()]]
    all_data_array = np.concatenate(all_npydata)
    np.save(os.path.join(model_dir, 'all_precipitation_npy_data.npy'), all_data_array)
    print(f"Data for {model_dir} has been saved")