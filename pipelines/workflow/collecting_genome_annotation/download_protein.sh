#!/bin/bash
# Download from NCBI protein sequences

# Paths
export genome_assembly=$1
export protOutput=$2
export symlink_directory=$3

echo "------------ ${genome_assembly} ------------"


if [[ ${genome_assembly} =~ GCA ]]; then
  PathLink=$(esearch -db assembly -query ${genome_assembly} | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank)
elif [[ ${genome_assembly} =~ GCF ]]; then
  PathLink=$(esearch -db assembly -query ${genome_assembly} | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq)
fi

PathLink=$(echo ${PathLink} | cut -d " " -f 1)

BASENAME=`basename ${PathLink}`

Path=https${PathLink:3}/${BASENAME}'_protein.faa.gz'

wget ${Path} -O ${symlink_directory}${BASENAME}'_protein.faa.gz'
gzip -d ${symlink_directory}${BASENAME}'_protein.faa.gz'

cd ${symlink_directory}
ln -s ${BASENAME}'_protein.faa' ${protOutput}
