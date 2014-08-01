#!/usr/bin/perl

################################
## input
my $cufflinks_gtf= $ARGV[0];## only "u" tagged transcripts; this input in gtf format' 	
my $common_splice_junc=$ARGV[1]; ## common_OP.bed	
my $oases_bed=$ARGV[2];		## oases or trinity output supported
print "Usage: ./analyze_novel_tr.pl cuff_u.gtf common_OP.bed cdna.bed\n";
### output
my $cufflinks_junc="cuff_u.junc.bed";	
my $oases_gtf="oases.gtf";		
my $oases_junc="oases.junc.bed";	
my $cuff_common_junc_oases_map="cuff_common_oases_junc_map";	
my $valid_cuff_list="valid_novel_transcripts";	
my $oases_tp_list="valid_oases_transcripts";	

##################################
## read gtf file in transcript and exon hash
my  ($cufflinks_transcript_hash, $cufflinks_exon_hash) = read_gtf($cufflinks_gtf);

## make junc_bed from gtf file
my $junc_flag=make_junc($cufflinks_gtf,$cufflinks_junc);

## read junc_bed in to hash
my %cufflinks_junc_hash = read_junc_bed($cufflinks_junc);

## read commonFP bed
my %common_splice_hash=read_junc($common_splice_junc);

## make gtf from oases cdna.bed
my ($oases_gtf_exon, $oases_gtf_junc, $oases_gtf_transcript_exon, $oases_gtf_transcript_junc, $oases_max_transcript_exon_hash)=make_gtf_oases($oases_bed,$oases_gtf,$oases_junc);

################ MAP CUFF ################
open(CUFF_MAP,">$cuff_common_junc_oases_map")||die "cannot create file";

my %oases_transcripts_mapped_with_cuff;

foreach my $transcript_id(keys %$cufflinks_transcript_hash )
{
	
	
	my @cuff_coord_arr= split(/_/, $$cufflinks_transcript_hash{$transcript_id});
	
#	print "$cuff_coord_arr[0]:$cuff_coord_arr[1]-$cuff_coord_arr[2]\n";
	
	
	my @cuff_exon;	### cufflinks exons for a transcript
	my @cuff_junc;	### cufflinks juncs for a transcript
	
	if(exists($$cufflinks_exon_hash{$transcript_id}))
	{
#		print $transcript_id,"\n";
		@cuff_exon=@{$$cufflinks_exon_hash{$transcript_id}}; ### exon coord array for a transcript
	}
	
	if(exists($cufflinks_junc_hash{$transcript_id}))
	{
#		print $transcript_id,"\n";
		@cuff_junc=@{$cufflinks_junc_hash{$transcript_id}}; ### junc coord array for a given transcript
		
	}
	
	my $common_junc_hit=0;			### count of cufflinks juncs have an exact match with common juncs for a transcript
	my %oases_junc_transcript_list;
	my %oases_junc_transcript_list_both_str;
	
	foreach my $cuff_junc_coords(@cuff_junc) ### search for match for each cuff junc for a transcript
	{
		
		$cuff_junc_coords=~ s/\t/_/g;  ### commom junc format: "chr_start_stop_strand"
#		print $cuff_junc_coords,"\n";
		my $cuff_junc_coords_opp_str=$cuff_junc_coords;  ### opposite strand for oases junc search
		
		#### search in common junc bed file####
		if(exists($common_splice_hash{$cuff_junc_coords})) ## exact match search in common junc 
		{
#			print $cuff_junc_coords,"\n";
			$common_junc_hit++;
		}
		####### search in oases junc ##########
		if(exists($$oases_gtf_junc{$cuff_junc_coords}))   ######"chr_start_stop_strand"
		{
			foreach(@{$$oases_gtf_junc{$cuff_junc_coords}})
			{
				$oases_junc_transcript_list{$_}++;
			}
		}
		
		if(exists($$oases_gtf_junc{$cuff_junc_coords_opp_str}))   ######"chr_start_stop_strand" for opposite strand match
		{
			foreach(@{$$oases_gtf_junc{$cuff_junc_coords_opp_str}})
			{
				$oases_junc_transcript_list_both_str{$_}++;
			}
		}
		
		$cuff_junc_coords_opp_str=change_coord_strand($cuff_junc_coords_opp_str);
		
		if(exists($$oases_gtf_junc{$cuff_junc_coords_opp_str}))   ######"chr_start_stop_strand" for opposite strand match
		{
			foreach(@{$$oases_gtf_junc{$cuff_junc_coords_opp_str}})
			{
				$oases_junc_transcript_list_both_str{$_}++;
				$oases_transcripts_mapped_with_cuff{$_}++;
			}
		}
	}
		print CUFF_MAP "Cuff:\t",$transcript_id,"\t";
		print CUFF_MAP "$cuff_coord_arr[0]:$cuff_coord_arr[1]-$cuff_coord_arr[2]\n";
		print CUFF_MAP "Junc:\t",scalar @cuff_junc,"\t",$common_junc_hit,"\nExon:\t",scalar @cuff_exon,"\n";
		
		########### Analyse Mapped Transcripts ################
		foreach $oases_transcripts(keys %oases_junc_transcript_list_both_str) ## FOR BOTH STRAND OASES
		{
			print CUFF_MAP "Oases:\t",$oases_transcripts,"\n";
			
			if(exists($$oases_gtf_transcript_junc{$oases_transcripts}))
			{
				
				my $oases_transcript_junc_count=scalar @{$$oases_gtf_transcript_junc{$oases_transcripts}};
				
				my $oases_transcript_junc_hit=0;
				my $oases_transcript_common_junc_hit=0;
				my $oases_transcript_junc_common_junc_hit=0;
			
				foreach my $oases_gtf_transcript_junc_coords(@{$$oases_gtf_transcript_junc{$oases_transcripts}}) ### list of juncs in a oases transcript
				{
					
					$oases_gtf_transcript_junc_coords=~ s/\t/_/g;

					my $oases_gtf_transcript_junc_coords_opp= change_coord_strand($oases_gtf_transcript_junc_coords);
					
#					print $oases_gtf_transcript_junc_coords,"\n";
#					print $oases_gtf_transcript_junc_coords_opp,"\n";
					
					foreach my $cuff_junc_coords(@cuff_junc)  ### for each cuff junc for a transcript
					{
#						print $cuff_junc_coords,"\n";
						
						
						if($cuff_junc_coords eq $oases_gtf_transcript_junc_coords)
						{
						
							
							$oases_transcript_junc_hit++;     ### cuff junc of a transcript map with a oases transcript junc
							
							if(exists($common_splice_hash{$oases_gtf_transcript_junc_coords})) 
							{
								$oases_transcript_junc_common_junc_hit++; ### mapping cuff junc, common junc and oases junc
							}
						}
						elsif($cuff_junc_coords eq $oases_gtf_transcript_junc_coords_opp) ## for opposite strand
						{
#							print $oases_gtf_transcript_junc_coords_opp,"\n";
							$oases_transcript_junc_hit++;
							
							if(exists($common_splice_hash{$oases_gtf_transcript_junc_coords_opp})) 
							{
								$oases_transcript_junc_common_junc_hit++;
							}
						}
					}
				
				
					if(exists($common_splice_hash{$oases_gtf_transcript_junc_coords}))
					{
						$oases_transcript_common_junc_hit++;
					}
					if(exists($common_splice_hash{$oases_gtf_transcript_junc_coords_opp}))
					{
						$oases_transcript_common_junc_hit++;
					}
				}
				print CUFF_MAP "Junc:\ttrascript: ",$oases_transcript_junc_count,"\tcuff_&_oases_junc: ",$oases_transcript_junc_hit,"\toases_&_common: ",$oases_transcript_common_junc_hit,"\tcuff_&_common_&_oases: ",$oases_transcript_junc_common_junc_hit,"\n";
			}
			
			if(exists($$oases_gtf_transcript_exon{$oases_transcripts}))
			{
				my $oases_transcript_exon_count=scalar @{$$oases_gtf_transcript_exon{$oases_transcripts}};
				my $oases_transcript_exon_hit=0;
				foreach my $oases_gtf_transcript_exon_coords(@{$$oases_gtf_transcript_exon{$oases_transcripts}})
				{
					$oases_transcript_exon_hit++;
				}
				print CUFF_MAP "Exon\t",$oases_transcript_exon_count,"\t",$oases_transcript_exon_hit,"\n";
			}		
		}
	
}
close(CUFF_MAP);
####################################################################

################## analyze cuff_common_oases_junc ##################
open(TP,">$valid_cuff_list");
open(OSTP,">$oases_tp_list");

my ($cuff_junc, $cuff_junc_common, $cuff_oases_hash)=read_cuff_common_oases_junc($cuff_common_junc_oases_map);

for my $cuff_id(keys %$cuff_junc)
{
#	print $cuff_id,"\t",$$cuff_junc{$cuff_id},"\t",$$cuff_junc_common{$cuff_id},"\n" unless(exists($$cuff_oases_hash{$cuff_id}));
	unless(exists($$cuff_oases_hash{$cuff_id}))
	{
		if($$cuff_junc{$cuff_id} == $$cuff_junc_common{$cuff_id} && $$cuff_junc{$cuff_id} != 0)
		{
			my $chr_line=chr_format($$cufflinks_transcript_hash{$cuff_id});
#			print TP $cuff_id,"\t",$$cufflinks_transcript_hash{$cuff_id},"\n"; ## for zero hit cuffs
			print TP $cuff_id,"\t",$chr_line,"\n"; ## for zero hit cuffs
		}
	}
}

my %cuff_with_two_oases_hit;
my %cuff_with_three_oases_hit;
my %cuff_with_more_oases_hit;

my ($oases_transcript_hash, $oases_exon_hash) = read_gtf($oases_gtf);

for my $cuff_id(keys %$cuff_oases_hash)
{
	my $cuff_junc_count=$$cuff_junc{$cuff_id};
	my $cuff_common_junc_count=$$cuff_junc_common{$cuff_id};
	
	if(scalar keys %{$$cuff_oases_hash{$cuff_id}} == 1) ## "1" Oases transcript mapped
	{
		for my $oases_id(keys %{$$cuff_oases_hash{$cuff_id}})
		{
#			print $oases_id,"\t",$$cuff_oases_hash{$pri}{$oases_id},"\n";

			my @arr=split(/\t/,$$cuff_oases_hash{$cuff_id}{$oases_id});
#			print $cuff_id,"\t",$oases_id,"\t",$arr[0],"\t",$arr[1],"\t",$arr[2],"\t",$arr[3],"\n";
			if($cuff_junc_count == $arr[3] && $cuff_junc_count == $arr[0])
			{
				
#				print TP $cuff_id,"\t",$$cufflinks_transcript_hash{$cuff_id},"\n";  ## for one hit exact match cuffs
				my $chr_line=chr_format($$cufflinks_transcript_hash{$cuff_id});
				print TP $cuff_id,"\t",$chr_line,"\n";
			}
			
			if($cuff_junc_count == $arr[3] && $arr[0] > $cuff_junc_count && $arr[2] == $arr[0]-1)
			{
#				print OSTP $oases_id,"\t",$$oases_transcript_hash{$oases_id},"\n";  ## for one hit correct oases
				my $chr_line1=chr_format($$cufflinks_transcript_hash{$cuff_id});
				my $chr_line2=chr_format($$oases_transcript_hash{$oases_id});
				print OSTP $oases_id,"\t",$chr_line2,"\n";  ## for one hit correct oases
				
			}
		}

	}
	
	if(scalar keys %{$$cuff_oases_hash{$cuff_id}} == 2) ## "2" Oases transcript mapped
	{
		for my $oases_id(keys %{$$cuff_oases_hash{$cuff_id}})
		{
#			print $$cuff_oases_hash{$cuff_id}{$oases_id},"\n";
			my @arr=split(/\t/,$$cuff_oases_hash{$cuff_id}{$oases_id});
			
			if($cuff_junc_count == $arr[3] && $arr[0] == $cuff_junc_count)
			{
#				print $cuff_id,"\t",$$cufflinks_transcript_hash{$cuff_id},"\n";  ## for two hit exact match cuffs
				$cuff_with_two_oases_hit{$cuff_id}=$$cufflinks_transcript_hash{$cuff_id};
			}
			
		}
		
	}
	
	if(scalar keys %{$$cuff_oases_hash{$cuff_id}} == 3) ## "3" Oases transcript mapped
	{
#		print $pri,"\n";
		for my $oases_id(keys %{$$cuff_oases_hash{$cuff_id}})
		{
			my @arr=split(/\t/,$$cuff_oases_hash{$cuff_id}{$oases_id});
			
			if($cuff_junc_count == $arr[3] && $arr[0] == $cuff_junc_count)
			{
#				print $cuff_id,"\t",$$cufflinks_transcript_hash{$cuff_id},"\n";  ## for three hit exact match cuffs
				$cuff_with_three_oases_hit{$cuff_id}=$$cufflinks_transcript_hash{$cuff_id};
			}
			
		}
	}
	
		if(scalar keys %{$$cuff_oases_hash{$cuff_id}} > 3) ## more than "3" Oases transcript mapped
	{
#		print $pri,"\n";
		for my $oases_id(keys %{$$cuff_oases_hash{$cuff_id}})
		{
			my @arr=split(/\t/,$$cuff_oases_hash{$cuff_id}{$oases_id});
			
			if($cuff_junc_count == $arr[3] && $arr[0] == $cuff_junc_count)
			{
#				print $cuff_id,"\t",$$cufflinks_transcript_hash{$cuff_id},"\n";  ## for three hit exact match cuffs
				$cuff_with_more_oases_hit{$cuff_id}=$$cufflinks_transcript_hash{$cuff_id};
			}
			
		}
	}
		

}
	
for my $cuff_id(keys %cuff_with_two_oases_hit)
{
#	print TP "$cuff_id\t$cuff_with_two_oases_hit{$cuff_id}\n"
	my $chr_line=chr_format($cuff_with_two_oases_hit{$cuff_id});
	print TP "$cuff_id\t$chr_line\n"
}

for my $cuff_id(keys %cuff_with_three_oases_hit)
{
#	print TP "$cuff_id\t$cuff_with_three_oases_hit{$cuff_id}\n"
	my $chr_line=chr_format($cuff_with_three_oases_hit{$cuff_id});
	print TP "$cuff_id\t$chr_line\n"
}

for my $cuff_id(keys %cuff_with_more_oases_hit)
{
#	print TP "$cuff_id\t$cuff_with_more_oases_hit{$cuff_id}\n"
	my $chr_line=chr_format($cuff_with_more_oases_hit{$cuff_id});
	print TP "$cuff_id\t$chr_line\n"
}
close(TP);;
close(OSTP);

#######################################################################################################################################

######################################################### subs #########################################################################

########################### read gtf in to transcript and exon hash ###########################
sub read_gtf ## input: gtf_file, output: %transcript_hash, %exon_hash
{
	my $gtf_file=$_[0];
	my  %gtf_transcript_hash;
	my  %gtf_exon_hash;
	
	open(F1,$gtf_file)||die "Provide input files $gtf_file\n";
	# open(F1,$gtf_file);
	while(<F1>)
	{
		chomp();
		my @arr=split(/\t/);
		my $transcript_id= $1 if(/transcript_id\s\"(\S+)\"/);
		my $exon_no= $1 if(/exon_number\s\"(\d+)\"/);
		if($arr[2]=~ /transcript/)
		{
			
			$gtf_transcript_hash{$transcript_id}=$arr[0]."_".$arr[3]."_".$arr[4]."_".$arr[6];
#			$gtf_transcript_hash{$transcript_id}=$arr[0].":".$arr[3]."-".$arr[4];
		}
		elsif($arr[2]=~ /exon/)
		{
			my $exon_value=$arr[0]."_".$arr[3]."_".$arr[4]."_".$arr[6];
			push @{$gtf_exon_hash{$transcript_id}},$exon_value;	
		}
	}
	close(F1);
	return (\%gtf_transcript_hash, \%gtf_exon_hash);
}


############## make junc from gtf ##################
sub make_junc ## input: cufflink_gtf, output cuflink_junc_bed file name
{
	my $gtf_file=$_[0];
	my $junc_file=$_[1];
	
	`rm -rf tmp`;
	`mkdir tmp`;
	
	open(GTF,$gtf_file)||die"no file";
	open(JUNC_BED,">$junc_file")||die "can't create $junc_file";
	my %hash;
	
	while(<GTF>)
	{
		chomp();
		push @{$hash{$1}},$_ if(/transcript_id \"(\S+)\"\;/);;
	}
	close(GTF);
	
	for my $id ( sort keys %hash) 
	{
		open(TMP_GTF_WRITE,">tmp/$id.gtf");
		
		my $transcript_start;
		my $transcript_end;
			
		print TMP_GTF_WRITE join "\n",@{ $hash{$id} };
		close(TMP_GTF_WRITE);
		
		my $pid=system(`gtf_juncs tmp/$id.gtf >tmp/$id.junc`);
		waitpid($pid,0);
	
		open(TMP_GTF,"tmp/$id.gtf");
		while(<TMP_GTF>)
		{		
			chomp();
			my @arr=split(/\t/);
			if($arr[2]=~ /transcript/)
			{
				$transcript_start=$arr[3];
				$transcript_end=$arr[4];
			}
		}
		close(TMP_GTF);
		
		open(TMP_JUNC,"tmp/$id.junc") ;
		while(<TMP_JUNC>)
		{
			chomp();
			my @arr=split(/\t/);
			print JUNC_BED $arr[0],"\t",$transcript_start,"\t",$transcript_end,"\t",$id,"\t",scalar @{$hash{$id}}-1,"\t",$arr[3],"\t",$arr[1],"\t",$arr[2],"\n";
		}
		close(TMP_JUNC);
	}

	close(JUNC_BED);
	return 1;
}

################# read junc.bed (transcript included) file in to junc hash ######################

sub read_junc_bed ## input: junc_bed file (transcript included), output: %junc_hash
{
	my $junc_file=$_[0];
	open(F2,$junc_file)||die "no $junc_file\n";

	my %junc_hash; ## Hash of Array of junc from cuff

	while(<F2>)
	{
		chomp();
		my @arr=split(/\t/);
		my $junc_value=$arr[0]."_".$arr[6]."_".$arr[7]."_".$arr[5];
		push @{$junc_hash{$arr[3]}},$junc_value;
	}
	close(F2);
	return %junc_hash;
}

############################## read junc file in to junc hash ######################################

sub read_junc ## input: junc file , output: %junc_hash
{
	my $junc_file=$_[0];
	open(F3,$junc_file)||die "no $junc_file\n";
	my %junc_hash;   ## Hash of common junc bed, "_" separated
	while(<F3>)
	{
		chomp();
		my @arr=split(/\t/);
		my $splice_value=$arr[0]."_".$arr[1]."_".$arr[2]."_".$arr[3];
		$junc_hash{$splice_value}++;
	}
	close(F3);
	return %junc_hash; 
}

############################# make gtf from oases cdna.bed ##########################################

sub make_gtf_oases ## input: cdna.bed, output oases.gtf file name, output oases_junc.bed file name
{
	my $oases_bed=$_[0];
	my $oases_gtf=$_[1];
	my $oases_junc=$_[2];
	
	
	my %oases_bed_hash;
	my %oases_max_length_bed_hash;
	my %oases_transcript_bed_hash;
	my %oases_max_transcript_exon_hash;
	
	open(CDNA_BED,$oases_bed)||die"no $oases_bed\n";
	
	while(<CDNA_BED>)
	{
		chomp();
		my @arr=split(/\t/);
		my $oases_key= $arr[0]."_".$arr[1]."_".$arr[2]."_".$arr[5];
		#print $arr[3],"\n";
		push @{$oases_bed_hash{$oases_key}},$arr[3];
		push @{$oases_transcript_bed_hash{$arr[3]}},$oases_key;
		push @{$oases_max_transcript_exon_hash{$arr[3]}},$oases_key; ## new
	}
	close(CDNA_BED);
	
	
	foreach my $oases_bed_coords(keys %oases_bed_hash)
	{
		my @oases_length;
		my $max_length=0;
		my $max_transcript;
	
#		foreach(@{$oases_bed_hash{$oases_bed_coords}})
#		{
#			if(/_Length_(\d+)$/)
#			{
#				if($1 > $max_length)
#				{
#					$max_length=$1;
#					$max_transcript= $_;
#				}
#			}
#		}
#		$oases_max_length_bed_hash{$max_transcript}=$oases_bed_coords;
	
#		if(exists($oases_transcript_bed_hash{$max_transcript}))
#		{
#			$oases_max_transcript_exon_hash{$max_transcript}= $oases_transcript_bed_hash{$max_transcript};
#			print join "\n",@{$oases_transcript_bed_hash{$max_transcript}},"\n";
#			print join "\n",@{$oases_max_transcript_exon_hash{$max_transcript}},"\n";
#		}
	}
	
	open(OASES_GTF,">$oases_gtf")||die "cant create $oases_gtf";
	open(OASES_JUNC,">$oases_junc")||die "cant create $oases_junc";
	
	##### New Oases Hashes #####
	my %oases_gtf_exon;
	my %oases_gtf_junc;
	my %oases_gtf_transcript_exon;
	my %oases_gtf_transcript_junc;
	#############################
	
	foreach my $keys(keys %oases_max_transcript_exon_hash)
	{
		my $plus;
		my $minus;
		
		foreach my $line(@{$oases_max_transcript_exon_hash{$keys}})
		{
			$line=~ s/_/\t/g;
			my @arr=split(/\t/,$line);
			
			$plus++ if($arr[3] eq "+");
			$minus++ if($arr[3] eq "-");
		}
		
		my @exon_arr;
		
		foreach my $line(@{$oases_max_transcript_exon_hash{$keys}})
		{
			$line=~ s/_/\t/g;
			my @arr=split(/\t/,$line);
			
			if($plus >= $minus)
			{
				push @exon_arr,$line if($arr[3] eq "+");
			}
			else
			{
			push @exon_arr,$line if($arr[3] eq "-");
			}
		}
		
		@exon_arr=sort_coord(@exon_arr);
		my @junc_arr=gtf_junc(@exon_arr);
		my $valid_junc=junc_check(@junc_arr);
		
		unless($valid_junc > 0)
		{
			my @first_line_arr=split(/\t/,$exon_arr[0]);
			my @last_line_arr=split(/\t/,$exon_arr[scalar @exon_arr - 1]);
			my $transcript_start;
			my $transcript_end;
			
			if($first_line_arr[3] eq "+")
			{
				$transcript_start=$first_line_arr[1];
				$transcript_end=$last_line_arr[2];
			}
			elsif($first_line_arr[3] eq "-")
			{
				$transcript_start=$last_line_arr[1];
				$transcript_end=$first_line_arr[2];
			}
			
			print OASES_GTF $first_line_arr[0],"\tOases\ttranscript\t",$transcript_start,"\t",$transcript_end,"\t1000\t",$first_line_arr[3],"\t.\tgene_id \"$keys\"; transcript_id \"$keys\";\n";
	
			for (my $i=0; $i< scalar @exon_arr;$i++)
			{
				my @exon_choord=split(/\t/,$exon_arr[$i]);
				
				print OASES_GTF $exon_choord[0],"\tOases\texon\t",$exon_choord[1],"\t",$exon_choord[2],"\t1000\t",$exon_choord[3],"\t.\tgene_id \"$keys\"; transcript_id \"$keys\"; exon_number \"",$i+1,"\";\n";

				my $gtf_exon_key=$exon_arr[$i];
				$gtf_exon_key=~ s/\t/_/g;
				
				push @{$oases_gtf_exon{$gtf_exon_key}},$keys;
				push @{$oases_gtf_transcript_exon{$keys}},$exon_arr[$i];
				
			}
			
			foreach(@junc_arr)
			{
				my $gtf_junc_key=$_;
				my @junc_coords=split(/\t/);
				print OASES_JUNC "$junc_coords[0]\t$transcript_start\t$transcript_end\t$keys\t",scalar @exon_arr,"\t$junc_coords[3]\t$junc_coords[1]\t$junc_coords[2]\n" ;
			
				push @{$oases_gtf_transcript_junc{$keys}},$gtf_junc_key;
			
				$gtf_junc_key=~ s/\t/_/g;
				push @{$oases_gtf_junc{$gtf_junc_key}},$keys;
			}
		}
	}
	close(OASES_JUNC);
	close(OASES_GTF);
	
	return (\%oases_gtf_exon,\%oases_gtf_junc,\%oases_gtf_transcript_exon,\%oases_gtf_transcript_junc,\%oases_max_transcript_exon_hash);
}

######################
sub sort_coord
{
	my @exon_un_sorted_arr=@_;
	my @exon_sorted_arr;
	# print $exon_un_sorted_arr[0],"\n";
	
	my @sort_order;
	my %exon_hash;
	my $chr;
	my $strand;
	
	
	foreach (@exon_un_sorted_arr)
	{	
		my @exon_line=split(/\t/);
		$chr=$exon_line[0];
		$strand=$exon_line[3];
		
		if( $strand eq "+")
		{
			# print $strand,"\n";
			push @sort_order,$exon_line[1];
			$exon_hash{$exon_line[1]}=$exon_line[2];
		}
		elsif($strand eq "-")
		{
			# print $strand,"\n";
			push @sort_order,$exon_line[2];	
			$exon_hash{$exon_line[2]}=$exon_line[1];
		}
	}
	if($strand eq "+")
	{
		@sort_order = sort {$a <=> $b} @sort_order;
	}
	elsif($strand eq "-")
	{
		@sort_order = sort {$b <=> $a} @sort_order;
	}
	
	# print $chr,"\n";
	
	if($strand eq "+")
	{
		foreach (@sort_order)
		{
			my $line_sorted1=$chr."\t".$_."\t".$exon_hash{$_}."\t".$strand;
			# print $exon_hash{$_},"\n";
			# print $line_sorted1,"\n";
			push @exon_sorted_arr, $line_sorted1;
		}
	}
	elsif($strand eq "-")
	{
		foreach (@sort_order)
		{
			my $line_sorted2=$chr."\t".$exon_hash{$_}."\t".$_."\t".$strand;
			# print $exon_hash{$_},"\n";
			# print $line_sorted2,"\n";
			push @exon_sorted_arr, $line_sorted2;
		}
	}
return @exon_sorted_arr;
}

##############################

sub gtf_junc
{
	my @exon_sorted_arr=@_;
	my $exon_length=scalar @exon_sorted_arr;
	
	my @value1;
	my @value2;
	my $chr;
	my $strand;
	my @junc_arr;
	
	foreach(@exon_sorted_arr)
	{
		my @exon_line=split(/\t/);
		push @value1,$exon_line[1];
		push @value2,$exon_line[2];
		$chr=$exon_line[0];
		$strand=$exon_line[3];
	}
	
	if($strand eq "+")
	{
		for(my $i=0;$i<$exon_length -1;$i++)
		{
			my $junc=$chr."\t".($value2[$i]-1)."\t".($value1[$i+1]-1)."\t".$strand;
			# print $junc,"\n";
			push @junc_arr, $junc;
		}
	}
	elsif($strand eq "-")
	{	
		for(my $i=0;$i<$exon_length -1;$i++)
		{
			my $junc=$chr."\t".($value2[$i+1]-1)."\t".($value1[$i]-1)."\t".$strand;
			# print $junc,"\n";
			push @junc_arr, $junc;
		}
	}
	return @junc_arr;
}

################################

sub junc_check
{
	my @junc_arr=@_;
	my $flag=0;
	foreach (@junc_arr)
	{
		my @junc_coords=split(/\t/);
		if(($junc_coords[2]-$junc_coords[1])>50000)
		{
			$flag++;
		}
	}
	return $flag;
}

#################################
sub change_coord_strand
{
	my $strand=$_[0];
	my $opp_strand=$strand;
	
	if($opp_strand=~/_\+$/)
	{
		$opp_strand=~ tr/+/-/;
	}
	elsif($opp_strand=~/_-$/)
	{
		$opp_strand=~ tr/-/+/;
	}
	return $opp_strand;
}

############################# read cuff_common_oases_junc_map file ####################################
sub read_cuff_common_oases_junc   ## input: cuff_common_oases_junc_map_file , output: %cuff_junc, %cuff_junc_common, %cuff_oases_hash
{
	my $cuff_common_oases_junc_map_file=$_[0];
	open(F4,$cuff_common_oases_junc_map_file)||die "no $cuff_common_oases_junc_map_file\n";
	
	my %cuff_junc;
	my %cuff_junc_common;
	my %cuff_oases_hash;
	
	my $cuff_id;
	my $oases_id;
	
	
	while(my $line=<F4>)
	{
		chomp($line);
		
		
		if($line=~ /^Cuff:\t(\S+)\t(\S+)/)
		{	
			$cuff_id=$1;
			
		}
		
		if($line=~ /Junc:\t(\S+)\t(\S+)/)
		{
			$cuff_junc{$cuff_id}=$1;
			$cuff_junc_common{$cuff_id}=$2;
		}
		
		if($line=~ /Oases:\t(\S+)/)
		{
			$oases_id=$1;
		}
		
		if($line=~ /^Junc:\ttrascript:\s(\d+)\tcuff_&_oases_junc:\s(\d+)\toases_&_common:\s(\d+)\tcuff_&_common_&_oases:\s(\d+)/)
		{
			my $value=$1."\t".$2."\t".$3."\t".$4;
			$cuff_oases_hash{$cuff_id}{$oases_id}=$value;
		}
	}
	return(\%cuff_junc, \%cuff_junc_common, \%cuff_oases_hash);
}

################################ chr format #############################

sub chr_format	## input: "_" separated line, output: chr format line
{
	my $input_line=$_[0];
	my @arr=split(/_/,$input_line);
	my $output_line=$arr[0].":".$arr[1]."-".$arr[2];
	return $output_line;
}
