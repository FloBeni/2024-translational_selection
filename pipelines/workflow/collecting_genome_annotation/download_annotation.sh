#!/bin/bash
# Download from NCBI annotations

# Paths
export genome_assembly=$1
export gffOutput=$2
export gtfOutput=$3
export symlink_directory=$4

echo "------------ ${genome_assembly} ------------"


if [[ ${genome_assembly} =~ GCA ]]; then
  PathLink=$(esearch -db assembly -query ${genome_assembly} | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank)
elif [[ ${genome_assembly} =~ GCF ]]; then
  PathLink=$(esearch -db assembly -query ${genome_assembly} | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq)
fi

PathLink=$(echo ${PathLink} | cut -d " " -f 1)

BASENAME=`basename ${PathLink}`

Path=https${PathLink:3}/${BASENAME}'_genomic.gff.gz'

wget ${Path} -O ${symlink_directory}${BASENAME}'_genomic.gff.gz'
gzip -d ${symlink_directory}${BASENAME}'_genomic.gff.gz'

cd ${symlink_directory}
ln -s ${BASENAME}'_genomic.gff' ${gffOutput}

/beegfs/banque/gtdrift/pipeline/logiciels/cufflinks-2.2.1.Linux_x86_64/gffread ${gffOutput} -T -o ${BASENAME}'_genomic.gtf'

ln -s ${BASENAME}'_genomic.gtf' ${gtfOutput}
