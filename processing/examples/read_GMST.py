import pandas as pd

df = pd.read_csv('../../../megafires_data/GMST/GMST_mon_seas_year.csv', skiprows=1, parse_dates={'time': ['Year']})
df = df.set_index('time')
df = df.iloc[:-1,:]

jan = df['Jan'].to_xarray().astype(float)
year = df['J-D'].to_xarray().astype(float)

dfs = pd.read_csv('../../../megafires_data/GMST/GMST_year_smooth.csv', skiprows=1, parse_dates={'time': ['Year']})
dfs = dfs.set_index('time')

years = dfs['No_Smoothing'].to_xarray().astype(float)
smooth = dfs['Lowess(5)'].to_xarray().astype(float)
