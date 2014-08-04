
## Project : Identification_of_Novel_transcripts

### Retrival of Novel Transcripts:
- To retrive *true novel transcripts* from RNA-Seq dataset use the perl script `analyze_novel_tr.pl`
- `analyze_novel_tr.pl` requires TopHat executable `gtf_juncs`set in the enviornment 
- The three files which go as input to the script :

>      1. common_OP.bed	 ## common over predicted splice junctions (common_OP) between three aligners
	 2. cuff_u.gtf 	   ## cufflinks transcripts tagged as 'u' by cuffcompare
	 3. cdna.bed	   ## splice junctions of de Novo assembled transcripts

- Run the script as follows : 

>      ./analyze_novel_tr.pl cuff_u.gtf common_OP.bed cdna.bed

###### Outputs :

>      valid_novel_transcripts  ## cufflinks assembled novel transcripts that has evidence from common over predicted junctions aligners junctions and Oases.
>      valid_oases_transcripts  ## oases assembled de Novo transcripts that have evidence from common_OP but NOT cufflinks.

- - -
### Generating Input files for `analyze_novel_transcript.pl` :

#### 1.  common_OP.bed:
---

After aligning the reads, collect the splice junctions from TopHat, MapSplice and SOAPsplice in bed format and should contain NO headers; as shown below:

>   	Tophat => junctions.bed	
	  =========================
>       tophat -p 16 -G ref/500_transcripts.gtf -r 300 --mate-std-dev 20 -o sim_run/50x_75/tophat_1.4.1 --keep-tmp ref/1 sim_data_set/sim_50x_75_1.fastq sim_data_set/sim_50x_75_2.fastq
> 	
>     	head -4 junctions.bed
> 	    1	17020	17305	JUNC00000001	14	-	17020	17305	255,0,0	2	35,73	0,212
> 	    1	17295	17675	JUNC00000002	21	-	17295	17675	255,0,0	2	73,70	0,310
>       1	17655	17742	JUNC00000003	21	-	17655	17742	255,0,0	2	74,9	0,78
>       1	17733	17988	JUNC00000004	23	-	17733	17988	255,0,0	2	9,74	0,181
##  
	
>      MapSplice => best_remapped_junction.bed
	 ==========================================
>      python mapsplice_segments.py 50x_75_mapsplice_all.cfg
> 	
>   	 # 50x_75_mapsplice_all.cfg: 
>      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
> 	  reads_file = sim_data_set/sim_50x_75_1.fastq,sim_data_set/sim_50x_75_2.fastq
> 	  chromosome_files_dir = ref/
>  	  Bowtieidx = ref/1 
> 	  output_dir = sim_run/50x_75/mapsplice_all
> 	  reads_format = FASTQ
> 	  paired_end = yes
> 	  segment_length = 25
> 	  junction_type = canonica
> 	  fusion_junction_type = canonical
> 	  full_running = yes
> 	  anchor_length = 8
> 	  remove_temp_files = yes
> 	  segment_mismatches = 1
> 	  splice_mismatches = 2
> 	  remap_mismatches = 3
> 	  min_intron_length = 1
> 	  max_intron_length = 500000
> 	  threads = 4
> 	  max_hits = 4
> 	  max_insert = 3
> 	  search_whole_chromosome = no
> 	  map_segment_directly = no
> 	  run_MapPER = no
> 	  do_fusion = no
> 	  do_cluster = no
> 	  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
> 	  head -4 best_remapped_junction.bed
> 	  1	17020	17304	JUNC_1	15	-	17020	17304	255,0,0	2	35,72	0,212
> 	  1	17295	17667	JUNC_2	26	-	17295	17667	255,0,0	2	73,62	0,310
> 	  1	17731	17987	JUNC_3	7	-	17731	17987	255,0,0	2	11,73	0,183
> 	  1	17990	18339	JUNC_4	65	-	17990	18339	255,0,0	2	71,72	0,277
##
	
>     SOAPsplice => soap_all.junc
	===============================
>     soapsplice -d ref/1_soap.index -1 sim_data_set/sim_50x_75_1.fastq -2 sim_data_set/sim_50x_75_2.fastq -r 75 -I 300 -o sim_run/50x_75/soap_all -f 2 -q 1
> 	
> 	   **Convert soap_all.junc => soap_all.bed**

>     awk '{if($4 =="fwd"){print $1"\t"$2"\t"$3"\t+\t"$5}}' $1/soap_all/soap_all.junc > $1/soap_all/soap_all.bed
> 	  awk '{if($4 =="rev"){print $1"\t"$2"\t"$3"\t-\t"$5}}' $1/soap_all/soap_all.junc >> $1/soap_all/soap_all.bed
> 	
> 	  head -4 soap_all.bed 
> 	  1	89830	236395	+	1
> 	  1	334297	342392	+	23
> 	  1	342549	350525	+	12
> 	  1	342603	350525	+	35
	

  Junction predicted by three aligners will be compared against all possible splice junctions present in the reference genome; which can be collected from reference gtf using TopHat executable `gtf_juncs`
	
> 	  gtf_juncs ref/500_transcripts.gtf > 500_transcripts.juncs 
> 	
> 	  head 500_transcripts.juncs
> 	  1	17054	17232	-
> 	  1	17367	17605	-
> 	  1	17741	17914	-
> 	  1	18060	18267	-

Run the script `map_junctions.pl` as follows :
                  
>     ./map_junctions.pl -j 500_transcripts.juncs -t junctions.bed -m best_remapped_junction.bed -s soap_all.bed
	
	Output: common_OP.bed

	
#### 2.  cuff_u.gtf :
---

Execute cufflinks to generate transcripts.gtf followed by cuffcompare to generate cuff_u.gtf 
> 	
> 	  cufflinks -g ref/1_500_transcripts.gtf -o sim_run/50x_75/cufflinks_g_out -p 16 sim_run/50x_75/tophat_1.4.1/accepted_hits.bam
> 	  cuffcompare -r ref/500_transcripts.gtf transcripts.gtf 

>   	**Output by cuffcompare => cuffcmp.transcripts.gtf.tmap**

> 	  awk '{if ($3=="u") {print "transcript_id" " " "\""$5"\""}}' cuffcmp.transcripts.gtf.tmap >u_ids
> 	  grep -f u_ids transcripts.gtf >cuff_u.gtf

	Output: cuff_u.gtf

#### 3. cdna.bed :
---
Run Oases and Splign to generate `cdna.bed`
######     Oases :
 

>      velveth directory k1,k2 data/reads.fa

   Where 'k1' and 'k2' are the range of k-mers. Replace 'k1' and 'k2' with '21' and '43' for dataset with read    length 75 bp. Now,create oases directory for every kmer as follows :
>
>       velvetg directory_k* -read_trkg yes
>       oases directory_k*

After running the previous process for different values of k, merge the results of all non-redundant assemblies contained in transcripts.fa in directory 'merged' using optimum K value = 27 

>      velveth merged/ 27 -long directory*/transcripts.fa
>      velvetg merged/ -read_trkg yes -conserveLong yes
>      oases merged/ -merge


######     Splign :

>     splign -mklds fasta_dir

>     formatdb -pF -oT -i transcripts.fa 
>     formatdb -pF -oT -i 1.fa   ## 1.fa is reference genome file used for simulated dataset 
>     compart -qdb cdna.fa -sdb 1.fa > cdna.compartments

>     splign -ldsdir fasta_dir -comps cdna.compartments 

>     awk '{if($1>0 && $8 != "-"){ if($8 >$9){print $3"\t"$9"\t"$8"\t"$2"\t"1"\t+"}else{print $3"\t"$8"\t"$9"\t"$2"\t"1"\t+"}}}' cdna.compartments >cdna.bed
>     awk '{if($1<0 && $8 != "-"){ if($8 >$9){print $3"\t"$9"\t"$8"\t"$2"\t"1"\t-"}else{print $3"\t"$8"\t"$9"\t"$2"\t"1"\t-"}}}' cdna.compartments >>cdna.bed

	
	Output: cdna.bed

- - -


