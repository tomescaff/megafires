import pandas as pd

f# define path of csv file with data
filepath = '../../../megafires_data/series/tmax_january_mon_mean_QN.csv'

# read the csv file using pandas
df = pd.read_csv(filepath, delimiter=",", decimal=",", parse_dates={'time': ['agno', 'mes', 'dia']})
df = df.rename({'valor':'tmax'}, axis='columns')
df = df.set_index('time')

# to xarray data array
da = df['tmax'].to_xarray()