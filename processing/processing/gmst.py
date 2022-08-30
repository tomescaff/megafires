import pandas as pd
from os.path import join, abspath, dirname

currentdir = dirname(abspath(__file__))
relpath = '../../../megafires_data/GMST'

def read_gmst_mon():
    filename = 'GMST_mon_seas_year.csv'
    filepath = join(currentdir, relpath, filename)
    df = pd.read_csv(filepath, skiprows=1, parse_dates={'time': ['Year']})
    df = df.set_index('time')
    df = df.iloc[:-1,:]
    return df

def get_gmst_jan():
    df = read_gmst_mon()
    return df['Jan'].to_xarray().astype(float)

def get_gmst_annual():
    df = read_gmst_mon()
    return df['J-D'].to_xarray().astype(float)

def get_gmst_annual_5year_smooth():
    filename = 'GMST_year_smooth.csv'
    filepath = join(currentdir, relpath, filename)
    df = pd.read_csv(filepath, skiprows=1, parse_dates={'time': ['Year']})
    df = df.set_index('time')
    return df['Lowess(5)'].to_xarray().astype(float)