"""
# This script is an attempt to view and extra info (variables) from ONE .nc dataset
# It will plot the images after running.
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cftime
from datetime import datetime


# path of the dataset
file_path = '/Users/zhangyiyi/Desktop/partdata/model1/historical/pr_day_ACCESS-CM2_historical_r1i1p1f1_gn_1950_v1.1.nc'
# access dataset
dataset = nc.Dataset(file_path, mode='r')

# view variables in the dataset
print("Variables in the dataset:", dataset.variables.keys())

# access each variable specifically
precipitation = dataset.variables['pr'][:]  # precipitation
time = dataset.variables['time'][:]  # time series
lat = dataset.variables['lat'][:]  # lat data
lon = dataset.variables['lon'][:]  # lon data

# print each data（shape）
print(f"Time dimension: {time.shape}")
print(f"Latitude dimension: {lat.shape}")
print(f"Longitude dimension: {lon.shape}")
print(f"Precipitation dimension: {precipitation.shape}")

# process time dimension to datatime type
time_units = dataset.variables['time'].units
calendar = dataset.variables['time'].calendar
print(f"Time units: {time_units}")
print(f"Calendar: {calendar}")

# transfer to time units
time_dates = nc.num2date(time, units=time_units, calendar=calendar)
time_dates_converted = [datetime(t.year, t.month, t.day) for t in time_dates]
print(f"First 5 time values: {time_dates_converted[:5]}")

# deal with missing value, to NaN
missing_value = dataset.variables['pr'].missing_value
precipitation = np.where(precipitation == missing_value, np.nan, precipitation)

# transfer pr from kg m-2 s-1 to mm/day
precipitation_mm_per_day = precipitation * 86400  # 86400 sec = 1 day
print("Precipitation converted to mm/day")

# plot pr distribution under a certain time
plt.figure(figsize=(10, 6))
plt.contourf(lon, lat, precipitation_mm_per_day[0, :, :], cmap='coolwarm')  # first time point 
plt.title(f"Precipitation Distribution (mm/day) at Time {time_dates_converted[0]}")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.colorbar(label='Precipitation (mm/day)')
plt.show()

# extra the time series data of a certain position
# set up the specific lat and lon
latitude = 23  # Guangzhou
longitude = 113  

# search the index close to the lat and lon
lat_idx = np.abs(lat - latitude).argmin()
lon_idx = np.abs(lon - longitude).argmin()

# extra the precipitation
precip_time_series = precipitation_mm_per_day[:, lat_idx, lon_idx]

# check if there are any valid data
if np.ma.is_masked(precip_time_series):
    print("Some or all of the data is masked (missing) in the precipitation time series.")
else:
    # print the first 5 values
    print(f"Precipitation time series: {precip_time_series[:5]}")

    # plot the pr time series
    plt.figure(figsize=(10, 6))
    plt.plot(time_dates_converted, precip_time_series)
    plt.title(f"Precipitation Time Series at lat={latitude}, lon={longitude} (mm/day)")
    plt.xlabel("Time")
    plt.ylabel("Precipitation (mm/day)")
    plt.grid(True)
    plt.show()

# close the dataset
dataset.close()
"""

# updated 09/20
# modified as using append() 
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import glob
import os

# path
data_dir = '/Users/zhangyiyi/Desktop/partdata/model1/historical/'

# access .nc
file_paths = glob.glob(os.path.join(data_dir, '*.nc'))

# initialize a empty list
all_results = []

# function of dealting with dataset and plot
def plot_time_series(lat, lon, precip, time_dates, latitude, longitude):
    # search index of lat and lon
    lat_idx = np.abs(lat - latitude).argmin()
    lon_idx = np.abs(lon - longitude).argmin()

    # extract pr time series
    precip_time_series = precip[:, lat_idx, lon_idx]
    
    # check data validity
    if np.ma.is_masked(precip_time_series):
        print("Some or all of the data is masked (missing) in the precipitation time series.")
    else:
        # plot pr time series
        plt.figure(figsize=(10, 6))
        plt.plot(time_dates, precip_time_series)
        plt.title(f"Precipitation Time Series at lat={latitude}, lon={longitude} (mm/day)")
        plt.xlabel("Time")
        plt.ylabel("Precipitation (mm/day)")
        plt.grid(True)
        plt.show()
    
    # save results as a dictionary and add them into the list
    result = {
        'latitude': latitude,
        'longitude': longitude,
        'time_dates': time_dates,
        'precip_time_series': precip_time_series.tolist()
    }
    all_results.append(result)  # by using append

# call all paths of files
for file_path in file_paths:
    # access datasets
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
    time_dates_converted = [datetime(t.year, t.month, t.day) for t in time_dates]

    # deal with NA value
    missing_value = dataset.variables['pr'].missing_value
    precipitation = np.where(precipitation == missing_value, np.nan, precipitation)

    # convert unit
    precipitation_mm_per_day = precipitation * 86400  # 86400 秒 = 1 天

    # plot pr distribution under a certain time
    plt.figure(figsize=(10, 6))
    plt.contourf(lon, lat, precipitation_mm_per_day[0, :, :], cmap='coolwarm')
    plt.title(f"Precipitation Distribution (mm/day) at Time {time_dates_converted[0]}")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.colorbar(label='Precipitation (mm/day)')
    plt.show()

    # call function to plot
    plot_time_series(lat, lon, precipitation_mm_per_day, time_dates_converted, 23, 113)
    
    # save all results into a .npy file
    np.save('all_results.npy', all_results)

    # close dataset
    dataset.close()

# print all results
#print(all_results)
load_npy = np.load('all_results.npy', allow_pickle=True)
print(load_npy)
