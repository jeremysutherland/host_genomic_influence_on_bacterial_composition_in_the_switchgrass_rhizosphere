import pandas as pd
import numpy as np
import re

# x_1 - x_4 are subsets of the HapMapv2 SNP matrix found at https://datadryad.org/stash/dataset/doi:10.5061/dryad.mp6cp
# x_1 = 500000 positions
# x_2 = 500001 positions
# x_3 = 500001 positions
# x_4 = 378586 positions

df = pd.read_csv("x_4", sep = '\t')
cols = df.columns[4:]

#### ADD ALT ALLELE COLUMN ####
df['allele1'] = df.Alleles.str[-2]
df['allele2'] = df.Alleles.str[-1]

def label_alt (row):
   if row['allele1'] != row['Reference']:
       return row['allele1']
   else:
       return row['allele2']

df['Alt'] = df.apply (lambda row: label_alt(row), axis=1)


#### Reorder positions to Ref/Alt ####
def order_alt(row):
    if row['tmp3'] == True :
        return row['tmp1'] + str("/") + row['tmp2']
    else:
        return row['tmp2'] + str("/") + row['tmp1']

for name in cols:
    df[['tmp1', 'tmp2']] = df[name].str.split("/",expand=True,)
    df.tmp2 = np.where(df.tmp2.isnull(), df.tmp1, df.tmp2)
    df['tmp3'] = np.where(df.tmp1.str[0] == df['Reference'], True, False)
    df[name] = df.apply (lambda row: order_alt(row), axis=1)

# df.iloc[:10, -10:]

df.to_csv('x_4_reorder.csv', index=False)


#### Fill Homozygotes Ref ####
def find_number(row):
    num = re.findall(r'[0-9]+',row)
    return "".join(num)

def call_hom(row):
    if row['tmp3'] == True:
        return row['Reference']
    else:
        return row[name]


for name in cols:
    df[['tmp1', 'tmp2']] = df[name].str.split("/",expand=True,)
    df.tmp2 = np.where(df.tmp2.isnull(), df.tmp1, df.tmp2)
    df['tmp3'] = np.where(df.tmp2.str[0] == df['Reference'], True, False)
    df['tmp4'] = df['tmp1'].apply(lambda row: find_number(row))
    df['tmp5'] = df['tmp2'].apply(lambda row: find_number(row))
    df['tmp6'] = df['tmp4'].astype(int) / df['tmp5'].astype(int)
    df[name] = df.apply (lambda row: call_hom(row), axis=1)


df.to_csv('x_4_reorder_hom_ref.csv', index=False)


#### Fill Homozygotes Alt ####
def call_hom2(row):
    if row['tmp3'] == False and row['tmp6'] == 1:
        return row['Alt']
    else:
        return row[name]


for name in cols:
    df[['tmp1', 'tmp2']] = df[name].str.split("/",expand=True,)
    df.tmp2 = np.where(df.tmp2.isnull(), df.tmp1, df.tmp2)
    df['tmp3'] = np.where(df.tmp2.str[0] == df['Reference'], True, False)
    df['tmp4'] = df['tmp1'].apply(lambda row: find_number(row))
    df['tmp5'] = df['tmp2'].apply(lambda row: find_number(row))
    df['tmp4'] = pd.to_numeric(df['tmp4'])
    df['tmp5'] = pd.to_numeric(df['tmp5'])
    df.tmp4 = np.where(df.tmp4.isnull(), 1, df['tmp4'])
    df.tmp5 = np.where(df.tmp5.isnull(), 1, df['tmp5'])
    df['tmp6'] = df['tmp4'].astype(int) / df['tmp5'].astype(int)
    df[name] = df.apply (lambda row: call_hom2(row), axis=1)

df.to_csv('x_4_reorder_hom_ref_alt.csv', index=False)

#### Fill Should be Homozygotes Ref ####
def call_hom3(row):
    if row['tmp3'] == False and row['tmp6'] > 0.75:
        return row['Reference']
    else:
        return row[name]

for name in cols:
    df[['tmp1', 'tmp2']] = df[name].str.split("/",expand=True,)
    df.tmp2 = np.where(df.tmp2.isnull(), df.tmp1, df.tmp2)
    df['tmp3'] = np.where(df.tmp2.str[0] == df['Reference'], True, False)
    df['tmp4'] = df['tmp1'].apply(lambda row: find_number(row))
    df['tmp5'] = df['tmp2'].apply(lambda row: find_number(row))
    df['tmp4'] = pd.to_numeric(df['tmp4'])
    df['tmp5'] = pd.to_numeric(df['tmp5'])
    df.tmp4 = np.where(df.tmp4.isnull(), 1, df['tmp4'])
    df.tmp5 = np.where(df.tmp5.isnull(), 1, df['tmp5'])
    df['tmp6'] = df['tmp4'].astype(int) / (df['tmp5'].astype(int) + df['tmp4'].astype(int))
    df[name] = df.apply (lambda row: call_hom3(row), axis=1)

df.to_csv('x_4_reorder_hom_ref_alt_ref2.csv', index=False)

#### Fill Should be Homozygotes Alt ####

df_length = len(df)
to_append = ['/'] * len(df.columns)
df.loc[df_length] = to_append

def call_hom4(row):
    if row['tmp3'] == False and row['tmp6'] < 0.25:
        return row['Alt']
    else:
        return row[name]

for name in cols:
    df[['tmp1', 'tmp2']] = df[name].str.split("/",expand=True,)
    df.tmp2 = np.where(df.tmp2.isnull(), df.tmp1, df.tmp2)
    df['tmp3'] = np.where(df.tmp2.str[0] == df['Reference'], True, False)
    df['tmp4'] = df['tmp1'].apply(lambda row: find_number(row))
    df['tmp5'] = df['tmp2'].apply(lambda row: find_number(row))
    df['tmp4'] = pd.to_numeric(df['tmp4'])
    df['tmp5'] = pd.to_numeric(df['tmp5'])
    df.tmp4 = np.where(df.tmp4.isnull(), 1, df['tmp4'])
    df.tmp5 = np.where(df.tmp5.isnull(), 1, df['tmp5'])
    df['tmp6'] = df['tmp4'].astype(int) / (df['tmp5'].astype(int) + df['tmp4'].astype(int))
    df[name] = df.apply (lambda row: call_hom4(row), axis=1)

df.to_csv('x_4_reorder_hom_ref_alt_ref2_alt2.csv', index=False)

#### Fill Heterozygotes ####

def call_het(row):
    if row['tmp6'] <= 0.75 and row['tmp6'] >= 0.25:
        return row['Alleles']
    else:
        return row[name]

for name in cols:
    df[['tmp1', 'tmp2']] = df[name].str.split("/",expand=True,)
    df.tmp2 = np.where(df.tmp2.isnull(), df.tmp1, df.tmp2)
    df['tmp3'] = np.where(df.tmp2.str[0] == df['Reference'], True, False)
    df['tmp4'] = df['tmp1'].apply(lambda row: find_number(row))
    df['tmp5'] = df['tmp2'].apply(lambda row: find_number(row))
    df['tmp4'] = pd.to_numeric(df['tmp4'])
    df['tmp5'] = pd.to_numeric(df['tmp5'])
    df.tmp4 = np.where(df.tmp4.isnull(), 4, df['tmp4'])
    df.tmp5 = np.where(df.tmp5.isnull(), 1, df['tmp5'])
    df['tmp6'] = df['tmp4'].astype(int) / (df['tmp5'].astype(int) + df['tmp4'].astype(int))
    df[name] = df.apply (lambda row: call_het(row), axis=1)

df.drop(df.tail(1).index,inplace=True) # drop last n rows
df = df.iloc[: , :-9]

df.to_csv('x_4_processed.csv', index=False)

