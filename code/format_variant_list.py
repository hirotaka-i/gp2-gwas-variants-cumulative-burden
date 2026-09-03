# Make bed files for plink2 to slice the pfile + score file to calculate PRS
## Bed file format:
# 1	109345810	109345810
# 1	154925709	154925709
## Score file format:
# chr1:154925709:G:C	C	0.2812
# chr1:155162560:G:A	A	0.6068
# Output in the TEMP_DIR

import pandas as pd
temp_dir = 'temp'

# GP2_Nalls2019_hg38
d=pd.read_csv('variant_list/GP2_Nalls2019_hg38.txt', sep='\t', names=['ID', 'A1', 'Effect'])
d_score = d.to_csv(f'{temp_dir}/GP2_Nalls2019_hg38.score', index=False, header=False, sep='\t')
d1 = d['ID'].str.split(':', expand=True)
d1.columns = ['chr', 'pos', 'ref', 'alt']
d1['chr'] = d1['chr'].str.replace('chr', '').astype(int)
d1['pos'] = d1['pos'].astype(int)
# If duplicated chr:pos return the duplicated chr pos and reduce the duplicated chr:pos to one row
d1_uniq = d1.drop_duplicates(subset=['chr', 'pos'])
if len(d1) != len(d1_uniq):
    print(f'Warning: {len(d1) - len(d1_uniq)} duplicated chr:pos found in GP2_Nalls2019_hg38.txt. Only the first one will be used.')
d_bed = d1_uniq[['chr', 'pos', 'pos']]
d_bed.to_csv(f'{temp_dir}/GP2_Nalls2019_hg38.bed', index=False, header=False, sep='\t')


# Leonard_media-2.csv
d = pd.read_csv('variant_list/Leonard_media-2.csv').dropna(subset=['SNP ID'])
d = d.rename(columns={
    'SNP ID': 'SNP', 'Effect allele': 'Effect_Allele',
    'Beta, all studies': 'all',
    'Beta, case-control': 'all_cc',
    'Beta, biobank': 'all_biobank',
})
for score_col in ['all', 'all_cc', 'all_biobank']:
    d[['SNP', 'Effect_Allele', score_col]].to_csv(
        f'{temp_dir}/Leonard_{score_col}.score',
        index=False, header=False, sep='\t')
d1 = d['SNP'].str.split(':', expand=True)
d1.columns = ['chr', 'pos', 'ref', 'alt']
d1['chr'] = d1['chr'].str.replace('chr', '').astype(int)
d1['pos'] = d1['pos'].astype(int)
d1_uniq = d1.drop_duplicates(subset=['chr', 'pos']).sort_values(by=['chr', 'pos']).reset_index(drop=True)   
if len(d1) != len(d1_uniq):
    print(f'Warning: {len(d1) - len(d1_uniq)} duplicated chr:pos found in Leonard_media-2.csv. Only the first one will be used.')
d1_uniq[['chr', 'pos', 'pos']].to_csv(f'{temp_dir}/Leonard_risk.bed', index=False, header=False, sep='\t')


# combine two bed files
bed1 = pd.read_csv(f'{temp_dir}/GP2_Nalls2019_hg38.bed', sep='\t', names=['chr', 'start', 'end'])
bed2 = pd.read_csv(f'{temp_dir}/Leonard_risk.bed', sep='\t', names=['chr', 'start', 'end'])
bed_combined = pd.concat([bed1, bed2]).drop_duplicates().sort_values(by=['chr', 'start', 'end']).reset_index(drop=True) 
bed_combined.to_csv(f'{temp_dir}/GP2_combined.bed', index=False, header=False, sep='\t')