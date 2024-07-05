#!/bin/bash
# Download from NCBI genomes, annotations and protein sequences

# Paths
export file_path=$1
export genome_path=$2
export number_core=$4


nb_chevron=$(grep -c "^>" ${genome_path})
nb_by=$(echo "$((${nb_chevron}/1000))")
echo ${nb_by}
#nb_by=${nb_chevron}
#echo ${nb_by}
if [ ${nb_by} -eq 0 ]; then
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%3==0){file=sprintf("'${file_path}'split_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < ${genome_path}
else
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%'${nb_by}'==0){file=sprintf("'${file_path}'split_%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < ${genome_path}
fi


#echo ${file_path}/${file}

list_file=$(ls ${file_path})
echo ${list_file}
echo ${number_core}
(
for file in ${list_file}
do
((i=i%number_core)); ((i++==0)) && wait
echo ${file}
tRNAscan-SE -E -o ${file_path}/${file}_tRNA.tab -S off ${file_path}/${file} &
done
wait
)

find ${file_path} -type f -name *_tRNA.tab -exec cat {} + > $3

