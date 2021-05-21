#module load blast
#MHCII
cd MHC2
#for black swan
#highest scaffolds selected as the genome
#makeblastdb -in blackswan.fasta -dbtype nucl -out bs
for i in mhc2*cons.fa
do
sp=bs
base=$(basename $i ".cons.fa")
echo "running blastn against $i"
blastn -evalue 1e-5 -out ${base}.${sp}.hits -num_threads 24 -query $i -outfmt 6 -perc_identity 0.6 -db ${sp}
cat ${base}.${sp}.hits
done

#for mute swan (chormosomes only)
#makeblastdb -in muteswan.fasta -dbtype nucl -out ms
for i in mhc2*cons.fa
do
sp=ms
base=$(basename $i ".cons.fa")
echo "running blastn against $i"
blastn -evalue 1e-5 -out ${base}.${sp}.hits -num_threads 24 -query $i -outfmt 6 -perc_identity 0.6 -db ${sp}
cat ${base}.${sp}.hits
done
