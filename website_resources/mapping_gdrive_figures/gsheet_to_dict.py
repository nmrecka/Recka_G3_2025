import pandas as pd

def gsheet_to_df(sheet_id, sheet_name):
    url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}'
    df = pd.read_csv(url)
    return df

SHEET_ID = ''
SHEET_NAME = 'Sheet1'
SAVE_PATH = 'data.js'
TRIM_KEYS_RIGHT = '_dermisVlnPlot.png' #'_otherVlnPlot.png' '_EpidermisVlnPlot.png', '_dermisVlnPlot.png'

df = gsheet_to_df(SHEET_ID, SHEET_NAME)

#get the name of the columns
key_col_name = df.columns[0]
value_col_name = df.columns[1]

if TRIM_KEYS_RIGHT != '':
    df[key_col_name] = df[key_col_name].str.replace(TRIM_KEYS_RIGHT, '')

df = df.sort_values(by=key_col_name)

#create a new text file to save the data into as a javascript dictionary
with open(SAVE_PATH, 'w') as f:

    #write the start of the dictionary
    f.write('    const data = {\n')

    #loop through the rows of the dataframe
    for index, row in df.iterrows():
        f.write(f'      "{row[key_col_name]}": "{row[value_col_name]}",\n')

    #write the end of the dictionary
    f.write('};')

print(f'Data saved to: {SAVE_PATH}')

