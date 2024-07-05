import gffpandas.gffpandas as gffpd
import pandas as pd
import sys
from collections import Counter
import Bio.SeqIO as IO
from itertools import product
from os.path import isfile, join

print(sys.argv)
cds_from_genomic_path = sys.argv[1]
gff_path = sys.argv[2]
codon_usage_path = sys.argv[3]

li = ['A', 'T', 'G' , 'C']
combi = ["combi"]
for comb in product(li, repeat=3):
    combi.append(''.join(comb))
combi = combi[1:66]


codon_usage_table = pd.DataFrame()
record_dict = IO.to_dict(IO.parse(cds_from_genomic_path, "fasta"))
for key in record_dict.items():
    if "protein_id=" in key[1].description:
        protein_id = key[1].description.split( "protein_id=" )[1].split("]")[0]
        length_cds = len(key[1].seq)
        vect_codon = Counter([key[1].seq[i:i + 3].upper() for i in range(0, len(key[1].seq)-2, 3)])
        freq_codon = [vect_codon[codon] for codon in combi]
        vect_pos3 = Counter([key[1].seq[i + 2].upper() for i in range(0, len(key[1].seq)-2, 3)])
        freq_gc3 = [vect_pos3[nucl] for nucl in li]
        start_codon = ''.join(key[1].seq[0:3].upper())
        end_codon = ''.join(key[1].seq[len(key[1].seq)-3:len(key[1].seq)].upper())
        data_prot = pd.DataFrame([[length_cds, protein_id] + freq_codon + freq_gc3 + [start_codon , end_codon]], columns=["length_cds", "protein_id"] + combi + ['A3', 'T3', 'G3' , 'C3' , 'start_codon','end_codon'] )
#        codon_usage_table = codon_usage_table.append(data_prot)
        codon_usage_table = pd.concat([codon_usage_table,data_prot],ignore_index=True)

annotation = gffpd.read_gff3( gff_path )
print( annotation.header )

cds_annotation = annotation.filter_feature_of_type(['CDS'])
cds_annotation = cds_annotation.attributes_to_columns()
cds_annotation = cds_annotation.loc[~cds_annotation.protein_id.duplicated()]
cds_annotation.set_index( cds_annotation.protein_id , inplace=True )

mRNA_annotation = annotation.filter_feature_of_type(['mRNA'])
mRNA_annotation = mRNA_annotation.attributes_to_columns()

cds_annotation["corresponding_gene"] = cds_annotation.Parent
if len(mRNA_annotation) != 0:
    mRNA_annotation.set_index( mRNA_annotation.ID , inplace=True )
    cds_annotation.loc[ cds_annotation.Parent.isin(mRNA_annotation.index),"corresponding_gene" ] = mRNA_annotation.loc[ cds_annotation.loc[ cds_annotation.Parent.isin(mRNA_annotation.index),"corresponding_gene"] ,"Parent" ].values


codon_usage_table["gene_id"] = cds_annotation.loc[codon_usage_table.protein_id,"corresponding_gene"].values
codon_usage_table["chromosome"] = cds_annotation.loc[codon_usage_table.protein_id,"seq_id"].values

codon_usage_table.to_csv(codon_usage_path, index=False, sep="\t")
