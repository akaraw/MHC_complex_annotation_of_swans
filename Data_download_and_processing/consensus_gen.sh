#Ref: https://academic.oup.com/gbe/article/13/2/evaa270/6050824?login=true#227709906
#We refere to the above paper to download MHC I and II representative sequences for each galliformes and anseriformes:
#Then we cluster them based on 80% identity using cd-hit as follows:
cd-hit-est -c 0.8 -i all_MHC1.fa -o /dev/stdout -G 1 #all_MHC?.fa contains all the repre. MHC I and II sequences from both galliformes and anseriformes seperately. 

#Then each clusert is seperated and the meber sequecnes for each cluser were aligned using muscle as follows:
muscle -in mhc2.clust0.fa  | tee mhc2.clust0.aln.fa | head
 
#then the aligned seq from each class of MHC in each sluter were used for consensus generation as follows:
# using the provided Rscript:

 
