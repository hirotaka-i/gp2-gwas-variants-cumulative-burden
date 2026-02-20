# Usage example: python code/detect_dup_TYPED_IMPUTED.py --pvar working/hirotaka/prs/temp/AFR_interest_chrpos.pvar
# This script is used to resolve the TYPED and IMPUTED duplication issue such as following:
# #CHROM	POS	ID	REF	ALT	INFO
# 4 11023507    chr4:11023507:C:T   C   T   IMPUTED;AF=0.603904;MAF=0.396096;AVG_CS=0.981919;R2=0.94355
# 4 11023507    chr4:11023507:T:C   T   C   TYPED;AF=0.395802;MAF=0.395802;AVG_CS=1;R2=1
# We prioritize the TYPED over IMPUTED and remove the IMPUTED in these cases

## However - if it is a truly multi-allelic site, we should keep either IMPUTED or TYPED. No mix for consistency. 
## e.g. a multi allelic site with C 0.5, G 0.3, T 0.2
## REF(C) ALT(G) IMPUTED AF=0.3
## REF(C) ALT(T) IMPUTED AF=0.2
## REF(G) ALT(C) TYPED AF=0.5 --> This going to be REF(C) ALT(G) AF=0.5 in normalization process but that is wrong.


import pandas as pd
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Generate exclude list to resolve TYPED and IMPUTED duplication in a .pvar file")
parser.add_argument('--pvar', help='Input pvar file path', required=True)
args = parser.parse_args()
pvar_path = args.pvar


# Function to find the header row number in the pvar file
def find_header_row(filename):
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if line.startswith('#CHROM'):
                return i
    return None

# Read the input file
header_row = find_header_row(pvar_path)
if header_row is not None:
    pvar = pd.read_csv(pvar_path, sep='\s+', header=header_row)
    all_ids = pvar['ID']

    # detect TRUE multi-allelic sites which has multiple TYPE or IMPUTED variants
    ##  only keep the either (TYPED or IMPUTED) with highest allelic numbers. No mix for consistency.
    pvar['IMPUTED_flag'] = pvar['INFO'].str.contains('IMPUTED')*1 # TYPED:IMPUTED, IMPUTED is 1, TYPED is 0
    t = pvar.groupby(['#CHROM', 'POS', 'IMPUTED_flag']).size()
    # sort by the number of variants in each site and the IMPUTED_flag and then take the largest size and lowest IMPUTED_flag
    t = t.reset_index().sort_values(by=[0, 'IMPUTED_flag'], ascending=[False, True]) # if the same number of variants, keep the TYPED
    t = t.drop_duplicates(subset=['#CHROM', 'POS'], keep='first')
    pvar = pvar[pvar[['#CHROM', 'POS', 'IMPUTED_flag']].apply(tuple, axis=1).isin(t[['#CHROM', 'POS', 'IMPUTED_flag']].apply(tuple, axis=1))]
    left_ids = pvar['ID']

    # 2. Find the removed IDs
    rm_ids = set(all_ids) - set(left_ids)
    rm_ids = list(rm_ids)
    rm_ids.sort()

    # Save the removing list
    output_path = pvar_path.replace('.pvar', '_exclude.txt')
    with open(output_path, 'w') as f:
        for item in rm_ids:
            f.write("%s\n" % item)
    print(f'{len(rm_ids)} IMPUTED SNPs should be removed because a TYPED SNP is available -> IDs saved at {output_path}')

else:
    print("Failed to find header row starting with '#CHROM'. Please check the file format.")