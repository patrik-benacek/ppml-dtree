import pandas as pd
import numpy as np
import time
import os
import sys

sys.path.append('config')
from config import (TARGET, LEADTIME, STATION, DATA_DIR)
                    
def load_dataset(file):
    """Read data from CSV file."""
    fname = os.path.join(DATA_DIR, file)
    try:
        return pd.read_csv(fname, parse_dates=True, index_col="date").drop(columns=['lon','lat','station'])
    except OSError:
        print("Could not open/read file:", fname)
        sys.exit()

def rename_station(df):
    """Transform station names."""

    df['station_names'] = df.station_names.apply(lambda x: x.rstrip().lstrip().replace(' / ','-'))
    return df

def read_dataset():
    """Get input dataset for particular experiment. 
    Output: 
        * dataset
        * expname: target_leadtime"""

    # Load dataset
    data = load_dataset(f"data_{TARGET}_ff{LEADTIME}.zip")

    # Experiment definition
    expname = f"{TARGET}_ff{LEADTIME}_{STATION}"

    # Rename station name
    data = rename_station(data)

    # Drop cases with nan observations
    data = data.dropna(subset=['obs'])

    # Get station data
    print(f"Run model for {STATION} station.")
    print(30*"-")
    if STATION=='all':
        return data, expname
    else:
        try:
            return data[data.station_names==STATION], expname
        except OSError:
            print("Station {} not exist. Available stations: \n{}".format(STATION, data.station_names.unique()))
            sys.exit()
