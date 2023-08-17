#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Yianna Bekris
Date: August 16th, 2023

This script calculates the following from gridMET data: 
    
    - Average annual mean temperature (average is the mean)
    - Average annual minimum temperature (Tmin)
    - Average annual maximum temperature (Tmax)
    - Average annual daily temperature range
    - Growing season; defined as number of days over 5 C -- at least 6 days after Jan 1st and 6 days < 5C after July 1st 
    - Frost days; defined as the annual count of days when Tmin is < 0C
    - Summer days; defined as the annual count of days when Tmax is > 25C
    - Average annual extreme monthly range; defined as the monthly Tmax - monthly Tmin
    
    All calculations are done for the years 1991-2020.
    
gridMET data can be downloaded from:
    
    https://www.climatologylab.org/gridmet.html

Abatzoglou, J. T. (2013), Development of gridded surface meteorological data for ecological applications and modelling. 
Int. J. Climatol., 33: 121–131.

"""

## Import packages
import xarray as xr
import numpy as np
import glob
import xclim.indices as indices
import xclim.indicators as indicators
import xclim


## Specify the paths where the files are located
tmmx_path = "/data/singh/data/gridmet/tmax/"
tmmn_path = "/data/singh/data/gridmet/tmmn/"

## Read in and sort the matching files
infiles_max = sorted(glob.glob(tmmx_path + "*199[1-9]*.nc")) + sorted(glob.glob(tmmx_path + "*20[0-1][0-9]*.nc"))
infiles_min = sorted(glob.glob(tmmn_path + "*199[1-9]*.nc")) + sorted(glob.glob(tmmn_path + "*20[0-1][0-9]*.nc"))

## Define variables
concat_var = 'time'
var_name = 'air_temperature'

## Initiate empty lists for each metric
avg_t_list = [] # max and min averaged
TNn_list = [] # mean of minimum temperature
TXx_list = [] # mean of maximum temperature
T_range_list = [] # max - min temperature (average of daily range)
grow_seas_list = [] # number of days over 5 C -- at least 6 days after Jan 1st and 6 days < 5C after July 1st 
frost_days_list = [] # annual count of days when Tmin is < 0C
summer_days_list = [] # annual count of days when Tmax is > 25C
x_mo_range_list = [] # monthly max - monthly min temp

       
## Loop through and compute variables to append to lists
for max_entry, min_entry in zip(infiles_max, infiles_min):
           
    ## Keep attributes of original files
    xr.set_options(keep_attrs=True)

    ## Open file
    ds_max = xr.open_dataset(max_entry)
    ds_min = xr.open_dataset(min_entry)
    
    ## Convert longitude to -180 to 180
    ds_max = ds_max.assign_coords(lon=(((ds_max.lon + 180) % 360) - 180))
    ds_min = ds_min.assign_coords(lon=(((ds_min.lon + 180) % 360) - 180))
 
    ## Rename time dimension from day to time to match other data (CPC, chirps)        
    ds_max = ds_max.rename(name_dict={'day' : 'time'})
    ds_min = ds_min.rename(name_dict={'day' : 'time'})
    
    ## Make a copy of the max dataset for creating the average
    ds_avg = ds_max.copy()
    
    ## Daily average
    ds_avg[var_name] = (ds_max[var_name] + ds_min[var_name]) / 2
    
    ## Find the monthly minimum of tmin and monthly maximum of tmax
    dsn_min = ds_min.groupby('time.month').min(dim=concat_var)
    dsx_max = ds_max.groupby('time.month').max(dim=concat_var)

    ## Find the extreme monthly range
    extreme_monthly_range = dsx_max - dsn_min
    
    ## Find the mean of the range
    x_mo_range = extreme_monthly_range.mean(dim='month')
    
    ## Air temperature
    airt_array = ds_avg[var_name]
    airt_mean = airt_array.mean(dim=concat_var)

    ## Growing season
    growing_season = indices.growing_season_length(tas=airt_array, mid_date="07-01", freq="AS")
    
    ## Frost days
    frost_t_array = ds_min[var_name]
    frost_count = indices.frost_days(frost_t_array)
    
    ## Summer days
    summer_t_array = ds_max[var_name]
    summer_count = indices.tx_days_above(summer_t_array, thresh='20.0 degC', freq='YS')
    
    ## Average of daily T range
    daily_T_r = ds_max[var_name] - ds_min[var_name]
    T_r = daily_T_r.mean(dim=concat_var)

    ## Annual verage of max and tmin
    TNn = ds_min[var_name].mean(dim=concat_var)
    TXx = ds_max[var_name].mean(dim=concat_var)

    ## Daily T range
    daily_T_r = ds_max[var_name] - ds_min[var_name]
    T_r = daily_T_r.mean(dim=concat_var)
    
    ## Daily mean
    avg_t_list.append(airt_mean)
    TNn_list.append(TNn)
    TXx_list.append(TXx)
    T_range_list.append(T_r)
    grow_seas_list.append(growing_season)
    frost_days_list.append(frost_count)
    summer_days_list.append(summer_count)
    x_mo_range_list.append(x_mo_range)
    
    ## Close datasets
    ds_max.close()
    ds_min.close()
    

## Concatenate along time dimension, calculate trends, and save netCDF files
avg_t_concat = xr.concat(avg_t_list, dim=concat_var)
avg_t_concat.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_avg_annual_T_1991_2020.nc')
avg_t_trends = avg_t_concat.polyfit(dim=concat_var, deg=1, skipna=True)
avg_t_trends_a = avg_t_trends.polyfit_coefficients.sel(degree=1)
avg_t_trends_a.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_avg_T_trends_1991_2020.nc')

TNn_concat = xr.concat(TNn_list, dim=concat_var)
TNn_concat.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_avg_annual_min_T_1991_2020.nc')
TNn_trends = TNn_concat.polyfit(dim=concat_var, deg=1, skipna=True)
TNn_trends_a = TNn_trends.polyfit_coefficients.sel(degree=1)
TNn_trends_a.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_TNn_trends_1991_2020.nc')

TXx_concat = xr.concat(TXx_list, dim=concat_var)
TXx_concat.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_avg_annual_max_T_1991_2020.nc')
TXx_trends = TXx_concat.polyfit(dim=concat_var, deg=1, skipna=True)
TXx_trends_a = TXx_trends.polyfit_coefficients.sel(degree=1)
TXx_trends_a.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_TXx_trends_1991_2020.nc')

T_range_concat = xr.concat(T_range_list, dim=concat_var)
T_range_concat.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_avg_annual_T_range_1991_2020.nc')
T_range_trends = T_range_concat.polyfit(dim=concat_var, deg=1, skipna=True)
T_range_concat_trends_a = T_range_trends.polyfit_coefficients.sel(degree=1)
T_range_concat_trends_a.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_T_range_trends_1991_2020.nc')

grow_seas_concat = xr.concat(grow_seas_list, dim=concat_var)
grow_seas_concat.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_avg_annual_grow_seas_1991_2020.nc')
grow_seas_trends = grow_seas_concat.polyfit(dim=concat_var, deg=1, skipna=True)
grow_seas_trends_a = grow_seas_trends.polyfit_coefficients.sel(degree=1)
grow_seas_trends_a.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_grow_seas_trends_1991_2020.nc')

frost_days_concat = xr.concat(frost_days_list, dim=concat_var)
frost_days_concat.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_avg_annual_frost_days_1991_2020.nc')
frost_days_trends = frost_days_concat.polyfit(dim=concat_var, deg=1, skipna=True)
frost_days_trends_a = frost_days_trends.polyfit_coefficients.sel(degree=1)
frost_days_trends_a.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_frost_days_trends_1991_2020.nc')

summer_days_concat = xr.concat(summer_days_list, dim=concat_var)
summer_days_concat.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_annual_summer_days_1991_2020.nc')
summer_days_trends = summer_days_concat.polyfit(dim=concat_var, deg=1, skipna=True)
summer_days_concat_trends_a = summer_days_trends.polyfit_coefficients.sel(degree=1)
summer_days_concat_trends_a.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_summer_days_trends_1991_2020.nc')

x_mo_range_concat = xr.concat(x_mo_range_list, dim=concat_var)
x_mo_range_concat.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_annual_ext_month_range_1991_2020.nc')
x_mo_range_trends = x_mo_range_concat.polyfit(dim=concat_var, deg=1, skipna=True)
x_mo_range_concat_trends_a = x_mo_range_trends.air_temperature_polyfit_coefficients.sel(degree=1)
x_mo_range_concat_trends_a.to_netcdf('/data/singh/yianna/butterfly_nc_files/gridmet_ext_month_range_trends_1991_2020.nc')



