#!/usr/bin/env python
# coding: utf-8

#=========================================
# Retrieve ECMWF forecast data from TIGGE
#=========================================
# Based on example: 
# * https://software.ecmwf.int/wiki/display/WEBAPI/TIGGE+retrieval+efficiency
#
# ------------------------------
# Description:
# ------------------------------
#   Period:    2015-2019
#   Init time: 00 UTC
#   Lead time: 24/36/240 h
#   Resolution: 0.5° 
#   Area: 2E, 30E; 42N, 55N 
#   Variables: 
#   * surface:  
#       Convective available potential energy
#       Surface pressure
#       Surface sensible heat flux
#       Surface latent heat flux
#       Mean sea level pressure
#       10 metre U wind component
#       10 metre V wind component
#       2 metre temperature
#       2 metre dewpoint temperature
#       Surface net solar radiation
#       Surface net thermal radiation
#       Skin temperature
#       Soil Moisture
#       Total Cloud Cover
#       Total Precipitation
#   * Level 500 hPa
#       UV wind,  
#       Geopotential 
#   * Level 850 hPa:
#       Temperature
#       UV wind 
#       Speciﬁc humidity
# ------------------------------
# Author: Patrik Benacek 
# Email: benacek.p@czechglobe.cz
# ------------------------------

from ecmwfapi import ECMWFDataServer
from calendar import monthrange
import os
from pathlib import Path

server = ECMWFDataServer()

data_dir = "/home/benacek.p/TIGGE/forecasts/"
years  = [2015, 2016, 2017, 2018, 2019]

def retrieve_tigge_data():
    # Settings
    # Make savepath directory
    for year in years:
        os.makedirs(os.path.join(data_dir, str(year)), exist_ok=True)
        for month in range(1, 13):
            days = range(1, monthrange(year, month)[1]+1)
            for day in days:
                # Date specification
                cyear  = str(year)
                cmonth = "%02d" % month
                cday   = "%02d" % day
                date   = "{}-{}-{}".format(cyear, cmonth, cday)
                # Downloading 
                # Surface parameters
                target = os.path.join(data_dir, cyear, "ecmwf_srf_{}{}{}.nc".format(cyear, cmonth, cday))
                if not os.path.exists(target):
                    print("Request for: ", target)
                    try:
                        tigge_request(date, target)
                        print("Finish.")
                    except:
                        print("MARS error downloading. Skip this term by touch.")
                        Path(target).touch()

def tigge_request(date, target):
    '''
       A TIGGE request for ECMWF perturbed forecasts of T2M.
    '''
    server.retrieve({
        'origin'    : "ecmf",
        'levtype'   : "sfc",
        'number'    : mem_numbers,
        'expver'    : "prod",
        'dataset'   : "tigge",
        'step'      : "24/36/240",
        'grid'      : "0.5/0.5",
        'param'     : "59/134/146/147/151/165/166/167/168/176/177/235/228039/228164/228228",
        'area'      : "55/2/42/30", # N/W/S/E
        'time'      : "00",
        'date'      : date,
        'type'      : "pf",
        'class'     : "ti",
        'format'    : "netcdf",
        'target'    : target
    })

mem_numbers = ''.join([''.join([str(i) + "/" for i in range(1,50)]),'50'])
retrieve_tigge_data()
