#!/bin/bash

# What do I have to set up?
# $1 = pathway to the ncbi tables and CDS/protein sequences
# $2 = species abbreviation, something like "dari", "dmoj"
# last command (egrep), I have to change the specie name

#seq_cds.fa = at FTP server "cds_from_genomic.fa" file
#seq_prot.fa = at FTP server "translated_CDS.fa" file

set -e

for tab in $1*.loc; do
	name="${tab%.*}"
	output="${name}"
	awk '{print $7}' "$tab" | grep 'LOC' > "$output"_ids.lst
done
	
for ids in $1*_ids.lst; do
	name="${ids%_*}"
	output="${name}"
	grep -f "$ids" $1*.fna | sed 's/>//g' | egrep -v 'X2|X3|X4|X5|X6|X7|X8|X9|X10|X11|X12|X13|X14|X15|X16|X17|X18|X19' > "$output"_CDS.lst
	grep -f "$ids" $1*.faa | sed 's/>//g' |  egrep -v 'X2|X3|X4|X5|X6|X7|X8|X9|X10|X11|X12|X13|X14|X15|X16|X17|X18|X19' > "$output"_PT.lst
	seqtk subseq "$1"*.fna "$output"_CDS.lst > "$output"_CDS.fa
	seqtk subseq "$1"*.faa "$output"_PT.lst > "$output"_PT.fa
done

mkdir "$1"PT "$1"CDS
mv "$1"*CDS.fa "$1"CDS/
mv "$1"*PT.fa "$1"PT/

rm "$1"*.lst

for seq in "$1"PT/*.fa; do
	name="${seq%.*}"
	output="${name}"
	grep '>' "$seq" | sed 's/lcl.*gene=/#/' | sed 's/db_xref.*protein=/#/' | sed 's/protein_id=.*CDS/#/' | tr -d \[]# | sed 's/.$//' | sed 's/$/_dhyd/' > "$output".txt
done

cp change_PT_id.py "$1"PT/

cd "$1"PT/

python change_PT_id.py

cd ../../

mv "$1"PT/*.txt "$1"CDS/

for fasta in "$1"PT/*tidy.fa; do
	name="${fasta%_tidy*}"
	output="${name}"
	seqkit rmdup "$fasta" > "$output"_nodup.fa
done

mkdir "$1"PT/temp
mv "$1"PT/*PT.fa "$1"PT/*tidy.fa "$1"PT/*.py "$1"PT/temp

cp change_CDS_id.py "$1"CDS/

cd "$1"CDS/

python change_CDS_id.py

cd ../../

echo "Detecting ID's duplication CDS"

for fasta in "$1"CDS/*tidy.fa; do
	name="${fasta%_tidy*}"
	output="${name}"
	seqkit rmdup "$fasta" > "$output"_nodup.fa
done

mkdir "$1"CDS/temp
mv "$1"CDS/*CDS.fa "$1"CDS/*tidy.fa "$1"CDS/*.py "$1"CDS/*.txt "$1"CDS/temp



#egrep -v 'X2|X3|X4|X5|X6|X7|X8|X9|X10|X11|X12|X13|X14|X15|X16|X17|X18|X19'
# And then change the IDs with "change_id.py"
# Remove LOC duplications with seqkit rmdup
