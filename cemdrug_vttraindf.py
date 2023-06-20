import os
import time
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols

def save_protein_fasta(df, protein_list):
    '''
    将protein_accession集合A1对应的蛋白保存到FASTA文件，文件名称为protein_accession
    '''
    if not os.path.exists('key_protein'):
        os.makedirs('key_protein')
        
    protein_df = df[df['protein_accession'].isin(protein_list)][['protein_accession', 'protein_sequence']].drop_duplicates()
    for accession, sequence in protein_df.values:
        with open(f'key_protein/{accession}.fasta', 'w') as f:
            f.write(f'>{accession}\n{sequence}')

def save_compound_smiles(df, compound_list):
    comp_df = df[df['compound_chembl_id'].isin(compound_list)]
    comp_df[['canonical_smiles','compound_chembl_id']].drop_duplicates().to_csv('key_smiles.smi', sep='\t', index=False, header=False)

start_time = time.time()

# 设置工作路径
os.chdir('*****PATH***')

# 读取数据
data = pd.read_csv('protein_inhibitors.csv')
data = data.drop_duplicates(subset=['compound_chembl_id', 'protein_accession'])

# 剔除任何一列为空的数据
data.dropna(inplace=True)

# 筛选standard_units为nM的数据
data_nM = data[data.standard_units == 'nM']

# 根据standard_relation和standard_value对class列进行填充
def classify(row):
    if row['standard_relation'] == '=':
        if row['standard_value'] < 10000:
            return 'active'
        else:
            return 'inactive'
    elif row['standard_relation'] in ['>', '>=']:
        if row['standard_value'] >= 10000:
            return 'inactive'
        else:
            return 'unknown'
    elif row['standard_relation'] == '<':
        if row['standard_value'] <= 10000:
            return 'active'
        else:
            return 'unknown'
    else:
        return 'unknown'

data_nM['class'] = data_nM.apply(classify, axis=1)

# 提取需要的列
activity_df = data_nM[['molregno', 'compound_chembl_id', 'target_name', 'protein_accession', 'protein_sequence', 'target_type', 'canonical_smiles', 'class']]

# 去除 class 为 unknown 的结果
activity_df = activity_df[activity_df['class'] != 'unknown']

# 去除偏颇数据，假如都是active或者inactive则剔除。
count_df = activity_df.groupby(['protein_accession', 'class'])['class'].count().unstack(fill_value=0)
drop_proteins = count_df[(count_df['active'] == 0) | (count_df['inactive'] == 0) | (count_df['inactive'] < 52) | (count_df['active'] < 52) ].index
activity_df = activity_df[~activity_df['protein_accession'].isin(drop_proteins)]

# 保存化合物数量不少于200的蛋白 
count_df2 = activity_df.groupby('protein_accession').size().reset_index(name='count')
top_proteins = count_df2[count_df2['count'] >= 200]['protein_accession'].values  
top_df = activity_df[activity_df['protein_accession'].isin(top_proteins)]

# 保存结果
top_compounds = top_df['compound_chembl_id'].unique()
save_compound_smiles(top_df, top_compounds)
save_protein_fasta(top_df, top_proteins)
top_df.to_csv('chembl_cemdrug.csv', index=False)

end_time = time.time()
total_time = end_time - start_time
print(f"Total time elapsed: {total_time:.2f} seconds")

