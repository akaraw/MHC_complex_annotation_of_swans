#Using these scripts one can do a batch processing of the exon specific sequences at order level
for i in mhc1_exon?;do echo $i; conda activate base; cd-hit-est -c 0.8 -i $i -o ${i}.rep -G 1; conda activate base; muscle -in ${i}.rep | tee  \
${i}.aln; conda activate r; Rscript --vanilla \ 
consensus.R ${i}.aln ${i}.cons; cat ${i}.cons ; done
