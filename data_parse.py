"""
Gene names mapped to Uniprot IDs using https://www.uniprot.org/id-mapping, yielding data/genename_to_uniprotid.tsv.
Manually converted provided Screen_mutations_for_visualization_110223.xslx to Screen_mutations_for_visualization_110223.csv.
"""

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import seq1

# Read mapping data
df = pd.read_csv('data/genename_to_uniprotid.tsv', header=0, sep='\t')
df_reviewed = df[df['Reviewed'] == 'reviewed']
cols = list(df_reviewed.columns)
cols[0], cols[1], cols[-1] = 'Gene', 'Uniprot_ID', 'Sequence_wt'
df_reviewed.columns = cols

# Read mutation data
df_mut = pd.read_csv('data/Screen_mutations_for_visualization_110223.csv', header=0)
df_mut['aa_from'] = df_mut['Perturbation'].str[:3]
df_mut['aa_to'] = df_mut['Perturbation'].str[-3:]
df_mut['mut_pos'] = df_mut['Perturbation'].str[3:-3].astype(int)

cols_kept = ['Gene', 'Uniprot_ID', 'aa_from', 'aa_to', 'mut_pos', 'Sequence_wt']
df_mut = df_mut.merge(df_reviewed, how='left', on='Gene')[cols_kept]

# Add column for mutated sequences
seqs_mut = []
for idx, row in df_mut.iterrows():
    seq = [i for i in row['Sequence_wt']]
    mut_pos = row['mut_pos'] - 1
    aa_from = seq1(row['aa_from'])
    aa_from_seq = seq[mut_pos]

    if aa_from_seq == aa_from:
        print(f'Matching AA at mutated position {mut_pos}: {aa_from_seq} - AA mutated from: {aa_from}')
    else:
        print(f'Mismatch between AA at mutated position {mut_pos}: {aa_from_seq} - AA mutated from: {aa_from}')

    seq[mut_pos] = seq1(row['aa_to'])
    seq = ''.join([i for i in seq])

    seqs_mut.append(seq)

df_mut['Sequence_mut'] = seqs_mut

df_mut.to_csv('data/mut_seqs.tsv', header=True, sep='\t', index=False)

# Create DataFrame with combined mutations
genenames = list(set(df_mut['Gene']))

uniprotids = []
seqs_mut_combined = []
seqs_wt = []
for genename in genenames:
    df_gene = df_mut[df_mut['Gene'] ==  genename]
    seq = [i for i in df_gene['Sequence_wt'].iloc[0]]
    # only keeps PIK3CD Y524H and not Y524C !
    for idx, row in df_gene.iterrows():
        mut_pos = row['mut_pos'] - 1
        seq[mut_pos] = seq1(row['aa_to'])

    seqs_mut_combined.append(''.join([i for i in seq]))
    uniprotids.append(df_gene['Uniprot_ID'].iloc[0])
    seqs_wt.append(df_gene['Sequence_wt'].iloc[0])

df_mut_combined = pd.DataFrame({'Gene': genenames,
                                'Uniprot_ID': uniprotids,
                                'Sequence_wt': seqs_wt,
                                'Sequence_mut_combined':seqs_mut_combined })
df_mut_combined.to_csv('data/mut_seqs_combined.tsv', header=True, sep='\t', index=False)
