#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Yianna Bekris
Date: August 16th, 2023

This script calculates the following from Chirps, CPC, and gridMET data:
    
    - Average annual precipitation (average is the mean)
    - Average annual wet days; defined as days with accumulated precipitation over 10 mm
    - Average annual dry days; defined as days with accumulated precipitation under 0.1 mm
    - Annual precipitation range; defined as the maximum minus minimum accumulated daily precipitation
    
    All calculations are done for the years 1991-2020.
    
Chirps: https://www.chc.ucsb.edu/data/chirps
CPC: https://psl.noaa.gov/data/gridded/data.cpc.globaltemp.html
gridMET: https://www.climatologylab.org/gridmet.html

"""

## Import packages
import xarray as xr
import numpy as np
import glob
import xclim.indices as indices
import xclim.indicators as indicators
import xclim

## File paths for precipitation datasets
# chirps_path = "/Users/yiannabekris/Desktop/chirps/*nc"
# precip_path = "/Users/yiannabekris/Documents/butterfly_figs/cpc_precip/*nc"
precip_path = "/Users/yiannabekris/Desktop/gridmet/precip/*nc"

precip_files = sorted(glob.glob(precip_path))

## Define variables
concat_var = 'time'
precip_var = 'precip'

## Initiate empty lists to store data
mean_precip_list = []
wet_days_list = [] # days over 10mm
cons_dry_days_list = [] # days under 0.1mm
precip_range_list = [] # max minus min in a year

       
## Loop through and take a daily mean
for entry in precip_files:

    
    # min_lon = -125.25
    # min_lat = 50.25
    # max_lon = -68.75
    # max_lat = 24.75
            
    ## Keep attributes of original files
    xr.set_options(keep_attrs=True)
    
    ## Open file
    precip = xr.open_dataset(entry)
    
    precip = precip.rename(name_dict={'day' : 'time'})
    precip = precip.rename(name_dict={'precipitation_amount' : 'precip'})
    
    ## Convert longitude from 0-360 10 -180-180
    # precip = precip.assign_coords(lon=(((precip.lon + 180) % 360) - 180))
    
    # ## Mask outside US if dataset is global
    # mask_lon = (precip.lon >= min_lon) & (precip.lon <= max_lon)
    # mask_lat = (precip.lat >= min_lat) & (precip.lat <= max_lat)
    # precip = precip.sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))

    ## Average
    precip_avg = precip.mean(dim='time')
    
    ## Wet Days # "Manually" calculated but indices function is below
    wet_days_bool = xr.where(precip[precip_var] > 10, 1, 0)
    wet_days = wet_days_bool.sum(dim='time')
    # wet_days = indices.wetdays_prop(precip[precip_var], thresh="10 mm/d", freq="YS")
    
    ## Dry Days # "Manually" calculated but indices function is below
    dry_days_bool = xr.where(precip[precip_var] < 0.1, 1, 0)
    dry_days = dry_days_bool.sum(dim='time')
    # dry_days = indices.dry_days(precip[precip_var], thresh="0.1 mm/d", freq="YS")
    
    ## Precip range
    ## First find the minumum and maximum values at each grid cell
    max_precip = precip.max(dim="time")
    min_precip = precip.min(dim="time")
    precip_range = max_precip.copy() ## Make a copy of the dataset
    
    ## Calculate the range
    precip_range[precip_var] = max_precip[precip_var] - min_precip[precip_var]

    ## Append to empty lists
    mean_precip_list.append(precip_avg)
    wet_days_list.append(wet_days)
    cons_dry_days_list.append(dry_days)
    precip_range_list.append(precip_range)
    
    ## Close precip dataset
    precip.close()


## Concatenate along time dimension # mean of everything run gridmet and cpc
mean_precip_concat = xr.concat(mean_precip_list, dim=concat_var)
mean_precip_concat.to_netcdf('/Users/yiannabekris/Desktop/climpact/ncfiles/gridmet_mean_precip_1991_2020.nc')
mean_precip_trends = mean_precip_concat.polyfit(dim="time", deg=1, skipna=True)
mean_precip_trends_a = mean_precip_trends.precip_polyfit_coefficients.sel(degree=1)
mean_precip_trends_a.to_netcdf('/Users/yiannabekris/Desktop/climpact/ncfiles/gridmet_mean_precip_trends_1991_2020.nc')

wet_days_concat = xr.concat(wet_days_list, dim=concat_var)
wet_days_concat.to_netcdf('/Users/yiannabekris/Desktop/climpact/ncfiles/gridmet_wet_days_1991_2020.nc')
wet_days_trends = wet_days_concat.polyfit(dim="time", deg=1, skipna=True)
wet_days_trends_a = wet_days_trends.polyfit_coefficients.sel(degree=1)
wet_days_trends_a.to_netcdf('/Users/yiannabekris/Desktop/climpact/ncfiles/gridmet_wet_days_trends_1991_2020.nc')

cons_dry_days_concat = xr.concat(cons_dry_days_list, dim=concat_var)
cons_dry_days_concat.to_netcdf('/Users/yiannabekris/Desktop/climpact/ncfiles/gridmet_dry_days_1991_2020.nc')
cons_dry_days_trends = cons_dry_days_concat.polyfit(dim="time", deg=1, skipna=True)
cons_dry_days_trends_a = cons_dry_days_trends.polyfit_coefficients.sel(degree=1)
cons_dry_days_trends_a.to_netcdf('/Users/yiannabekris/Desktop/climpact/ncfiles/gridmet_dry_days_trends_1991_2020.nc')

precip_range_concat = xr.concat(precip_range_list, dim=concat_var)
precip_range_concat.to_netcdf('/Users/yiannabekris/Desktop/climpact/ncfiles/gridmet_mean_precip_range_1991_2020.nc')
precip_range_trends = precip_range_concat.polyfit(dim="time", deg=1, skipna=True)
precip_range_concat_trends_a = precip_range_trends.precip_polyfit_coefficients.sel(degree=1)
precip_range_concat_trends_a.to_netcdf('/Users/yiannabekris/Desktop/climpact/ncfiles/gridmet_precip_range_trends_1991_2020.nc')





