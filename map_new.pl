#!/usr/bin/perl


use Getopt::Std;
our($opt_j, $opt_t, $opt_m, $opt_s);
getopts('j:t:m:s:');
(!$opt_j or !$opt_t or !$opt_m or !$opt_s) and die" -j juncs -t tophat_junc.bed -m mapsplice_junc.bed -s soapsplice_junc.bed\n";


#my $tophat_junction="tophat_novel/junctions.bed";
my $tophat_junction=$opt_t;

#my $mapsplice_junction="mapsplice_unmapped/best_remapped_junction.bed";
#my $soapsplice_junction="soapsplice_unmapped/run_unmapped.junc.bed";
#my $juncs_uniq="tophat_novel/novel_map_19_688.juncs_u";

my $mapsplice_junction=$opt_m;
my $soapsplice_junction=$opt_s;
my $juncs_uniq=$opt_j;

############## Read Input ###########################

open(F1,$tophat_junction); 
my %tophat_junctions;
my $tophat_total_j;

while(my $line=<F1>)
{
	
	chomp($line);
	my @arr=split(/\t/,$line);
	my $chr=$arr[0];
	my ($a1,$a2)=split(/\,/,$arr[10]);
	my $start=$arr[1]+$a1-1;
	my $stop=$arr[2]-$a2;
	my $strnd=$arr[5];
	my $key1=$chr."_".$start."_".$stop."_".$strnd;
	#print $key1,"\n";
	$tophat_junctions{$key1}++;
}
close(F1);

foreach my $t_k(keys %tophat_junctions)
{
	$tophat_total_j++;
}

open(F2,$juncs_uniq);
my %gtf_junctions;
my $known_junction;

while(my $line2=<F2>)
{
	
	chomp($line2);
	my @arr=split(/\t/,$line2);
	my $key2=$arr[0]."_".$arr[1]."_".$arr[2]."_".$arr[3];
	#print $key2,"\n";
	$gtf_junctions{$key2}++;
}
close(F2);

foreach my $gk(keys %gtf_junctions)
{
	$known_junction++;
}

open(F3,$mapsplice_junction);

my %mapsplice_junctions;
my $mapsplice_total_j;

while(my $line3=<F3>)
{
	
	chomp($line3);
	my @arr=split(/\t/,$line3);
	my $chr=$arr[0];
	my ($a1,$a2)=split(/\,/,$arr[10]);
	my $start=$arr[1]+$a1-1;
	my $stop=$arr[2]-$a2;
	my $strnd=$arr[5];
	my $key3=$chr."_".$start."_".$stop."_".$strnd;
	#print $key3,"\n";
	$mapsplice_junctions{$key3}++;
}
close(F3);

foreach my $m_k(keys %mapsplice_junctions)
{
	$mapsplice_total_j++;
}

open(F4,$soapsplice_junction);

my %soapsplice_junctions;
my $soapsplice_total_j;

while(my $line4=<F4>)
{
	chomp($line4);
	my @arr=split(/\t/,$line4);
	#print  $arr[0],"\t",$arr[1]-1,"\t",$arr[2]-1,"\t",$arr[3],"\n";
	my $chr=$arr[0];
	my $start=$arr[1]-1;
	my $stop=$arr[2]-1;
	my $strand=$arr[3];
	my $key4=$chr."_".$start."_".$stop."_".$strand;	
	#print $key4,"\n";
	$soapsplice_junctions{$key4}++;
}
close(F4);

foreach $s_k(keys %soapsplice_junctions)
{
	$soapsplice_total_j++;
}

###################################################
########### TN AND FP COUNT ################
print "Known Junction             :  $known_junction\n\n";
print "TopHat total junction      :  $tophat_total_j\n";
print "MapSplice total junction   :  $mapsplice_total_j\n";
print "SoapSplice total junction  :  $soapsplice_total_j\n\n";

#open(F5,">tophat_missed_true_negative.bed");
my $tophatTN;
my %tophatTN;

foreach my $key5(keys %gtf_junctions)
{
	$tophatTN++ unless(exists($tophat_junctions{$key5}));
	$tophatTP++ if(exists($tophat_junctions{$key5}));
	$tophatTN{$key5}++ unless(exists($tophat_junctions{$key5}));
#	print $key5,"\n"if(exists($tophat_junctions{$key5}));
	#print F5 $key5,"\n"unless(exists($tophat_junctions{$key5}));
}
#close(F5);


#open(F6,">tophat_over_predict_false_positive.bed");
my $tophatFP;
my %tophatFP;

foreach my $key6(keys %tophat_junctions)
{
	$tophatFP++ unless(exists($gtf_junctions{$key6}));
	$tophatFP{$key6}++ unless(exists($gtf_junctions{$key6}));
	#print F6 $key6,"\n"unless(exists($gtf_junctions{$key6}));
}
#close(F6);
print "\nTopHat:\n";
print "TopHat True Positive  : ",$tophatTP,"\n";
print "TopHat False Negative  : ",$tophatTN,"\n";
print "TopHat Over Prediction : ",$tophatFP,"\n";

#open(F7,">mapsplice_missed_true_negative.bed");

my $mapspliceTN;
my %mapspliceTN;

foreach my $key7(keys %gtf_junctions)
{
	$mapspliceTN++ unless(exists($mapsplice_junctions{$key7}));
	$mapspliceTP++ if(exists($mapsplice_junctions{$key7}));
	$mapspliceTN{$key7}++ unless(exists($mapsplice_junctions{$key7}));
	#print F7 $key7,"\n"unless(exists($mapsplice_junctions{$key7}));
}
#close(F7);



#open(F8,">mapsplice_over_predict_false_positive.bed");

my $mapspliceFP;
my %mapspliceFP;

foreach my $key8(keys %mapsplice_junctions)
{
	$mapspliceFP++ unless(exists($gtf_junctions{$key8}));
	$mapspliceFP{$key8}++ unless(exists($gtf_junctions{$key8})); 
	#print F8 $key8,"\n"unless(exists($gtf_junctions{$key8}));
}
#close(F8);

print "\nMapSplice:\n";
print "MapSplice True Positive  : ",$mapspliceTP,"\n";
print "MapSplice False Negative  : ",$mapspliceTN,"\n";
print "MapSplice Over Prediction : ",$mapspliceFP,"\n";

#open(F9,">soapsplice_missed_true_negative.bed");

my $soapspliceTN;
my %soapspliceTN;

foreach my $key9(keys %gtf_junctions)
{
	$soapspliceTN++ unless(exists($soapsplice_junctions{$key9}));
	$soapspliceTP++ if(exists($soapsplice_junctions{$key9}));
	$soapspliceTN{$key9}++ unless(exists($soapsplice_junctions{$key9}));
	#print F9 $key9,"\n"unless(exists($mapsplice_junctions{$key9}));
}
#close(F9);



#open(F10,">soapsplice_over_predict_false_positive.bed");

my $soapspliceFP;
my %soapspliceFP;

foreach my $key10(keys %soapsplice_junctions)
{
	$soapspliceFP++ unless(exists($gtf_junctions{$key10}));
	$soapspliceFP{$key10}++ unless(exists($gtf_junctions{$key10}));
	#print F10 $key10,"\n"unless(exists($gtf_junctions{$key10}));
}
#close(F10);

print "\nSoapSplice:\n";
print "SoapSplice True Positive  : ",$soapspliceTP,"\n";
print "SoapSplice False Negative  : ",$soapspliceTN,"\n";
print "SoapSplice Over Prediction : ",$soapspliceFP,"\n";
######################################################################
################  Minimize TN ########################################

my %tophatTN_mapspliceTN;
my $t_m_tn;
#my $t_m_s_tn;

foreach my $tophat_true_negative_coords(keys %tophatTN)
{
	$t_m_tn++ if(exists($mapsplice_junctions{$tophat_true_negative_coords}));
	$tophatTN_mapspliceTN{$tophat_true_negative_coords}++ unless(exists($mapsplice_junctions{$tophat_true_negative_coords}));
}
foreach my $tophatTN_mapspliceTN_coords(keys %tophatTN_mapspliceTN)
{
	$t_m_s_tn++ if(exists($soapsplice_junctions{$tophatTN_mapspliceTN_coords}));
}

#print "\nTN:\n";
#print "TopHat FN but mapped by MapSplice: ",$t_m_tn,"\n";
#print "TopHat and MapSplice FN but mapped by SoapSplice: ",$t_m_s_tn,"\n";
#print "Total FN reduced: ", $t_m_tn+$t_m_s_tn,"\/",$tophatTN,"\n";
#print "Total FN        : ", $tophatTN - ($t_m_tn+$t_m_s_tn),"\n";
#######################################################################
################# Count FP ############################################

open(F11,">common_OP.bed");
#open(F12,"tophat_novel/novel_juncs");
#my %novel_juncs;

#while(my $line12=<F12>)
#{
#	chomp($line12);
#	my @arr=split(/\t/,$line12);
#	my $key12=$arr[0]."_".$arr[1]."_".$arr[2]."_".$arr[3];
#	$novel_juncs{$key12}++;
#}
#close(F12);

my $s_fp;
my $t_m_fp;
my $t_s_fp;
my $m_s_fp;
my $t_m_s_fp;
#my $novel_t_m_s_fp;
#my $novel_t_m_fp;
#my $novel_t_s_fp;
#my $novel_m_s_fp;
#my $novel_s_fp;

foreach my $tophatFP_coords(keys %tophatFP)
{
	 #$t_m_fp++ if(exists($mapspliceFP{$tophatFP_coords}));
	if(exists($mapspliceFP{$tophatFP_coords}))
	{
		unless(exists($soapspliceFP{$tophatFP_coords}))
		{
			$t_m_fp++;
			#$novel_t_m_fp++ if(exists($novel_juncs{$tophatFP_coords}));
		}
	}
	 #$t_s_fp++ if(exists($soapspliceFP{$tophatFP_coords}));
	if(exists($soapspliceFP{$tophatFP_coords}))
	{
		unless(exists($mapspliceFP{$tophatFP_coords}))
		{
			$t_s_fp++;
			#$novel_t_s_fp++ if(exists($novel_juncs{$tophatFP_coords}));
		}
	}
}

foreach my $mapspliceFP_coords(keys %mapspliceFP)
{
	 #$m_s_fp++ if(exists($soapspliceFP{$mapspliceFP_coords}));
	if(exists($soapspliceFP{$mapspliceFP_coords}))
	{
		unless(exists($tophatFP{$mapspliceFP_coords}))
		{
			$m_s_fp++;
			#$novel_m_s_fp++ if(exists($novel_juncs{$mapspliceFP_coords}));
		}
	}
}

foreach my $tophatFP_coords(keys %tophatFP)
{
	 if(exists($mapspliceFP{$tophatFP_coords}))
	{
		#$t_m_s_fp++ if(exists($soapspliceFP{$tophatFP_coords}));
		if(exists($soapspliceFP{$tophatFP_coords}))
		{	#print $tophatFP_coords,"\n";
			$t_m_s_fp++;
			#$novel_t_m_s_fp++ if(exists($novel_juncs{$tophatFP_coords}));
		}
		print F11 join("\t",split(/\_/,$tophatFP_coords)),"\n" if(exists($soapspliceFP{$tophatFP_coords}));
	}
}

foreach my $soapspliceFP_coords(keys %soapspliceFP)
{
	unless(exists($mapspliceFP{$soapspliceFP_coords}))
	{
		unless(exists($tophatFP{$soapspliceFP_coords}))
		{
			$s_fp++;
			#$novel_s_fp++ if(exists($novel_juncs{$soapspliceFP_coords}));
		}
	}
}
close(F11);
print "\nOP:\n";
#print "TopHat,MapSplice,SoapSplice FP: $t_m_s_fp ($novel_t_m_s_fp)\n"; #harshal
print "TopHat,MapSplice,SoapSplice Common OP: $t_m_s_fp \n";
#print "TopHat,MapSplice FP           : $t_m_fp ($novel_t_m_fp)\n";
#print "TopHat,MapSplice FP           : $t_m_fp \n";
#print "TopHat,SoapSplice FP          : $t_s_fp ($novel_t_s_fp)\n";
#print "TopHat,SoapSplice FP          : $t_s_fp \n";
#print "MapSplice,SoapSplice FP       : $m_s_fp ($novel_m_s_fp) \n";
#print "MapSplice,SoapSplice FP       : $m_s_fp \n";
#print "SoapSplice FP                 : $s_fp ($novel_s_fp)\n";
#print "SoapSplice FP                 : $s_fp \n";

#############################################################################
##################### COMMON JUNCTIONS  #####################################
my $t_m_s_j;
my $t_m_j;
my $t_s_j;
my $m_s_j;
my $t_j;
my $m_j;
my $s_j;



foreach my $tophat_junctions_coords(keys %tophat_junctions)
{
	
	if(exists($mapsplice_junctions{$tophat_junctions_coords}))
	{
		$t_m_s_j++ if(exists($soapsplice_junctions{$tophat_junctions_coords}));
	}
	
	if(exists($mapsplice_junctions{$tophat_junctions_coords}))
	{
		$t_m_j++ unless(exists($soapsplice_junctions{$tophat_junctions_coords}));
	}
	
	unless(exists($mapsplice_junctions{$tophat_junctions_coords}))
	{
		$t_s_j++ if(exists($soapsplice_junctions{$tophat_junctions_coords}));
	}
	
	unless(exists($mapsplice_junctions{$tophat_junctions_coords}))
	{
		$t_j++ unless(exists($soapsplice_junctions{$tophat_junctions_coords}));
	}
}

foreach my $mapsplice_junctions_coords(keys %mapsplice_junctions)
{
	
	if(exists($soapsplice_junctions{$mapsplice_junctions_coords}))
	{
		$m_s_j++ unless(exists($tophat_junctions{$mapsplice_junctions_coords}));	
	}
	unless(exists($soapsplice_junctions{$mapsplice_junctions_coords}))
	{
		$m_j++ unless(exists($tophat_junctions{$mapsplice_junctions_coords}));
	}
}

foreach my $soapsplice_junctions_coords(keys %soapsplice_junctions)
{
	
		unless(exists($mapsplice_junctions{$soapsplice_junctions_coords}))
		{
			$s_j++ unless(exists($tophat_junctions{$soapsplice_junctions_coords}));	
		}
}


#print"\nJunctions:\n";

#print "TopHat,MapSplice,SoapSplice JN : $t_m_s_j\n";
#print "TopHat,MapSplice  JN           : $t_m_j\n";
#print "TopHat,SoapSplice JN           : $t_s_j\n";
#print "MapSplice,SoapSplice JN        : $m_s_j\n";
#print "TopHat JN                      : $t_j\n";
#print "MapSplice  JN                  : $m_j\n";
#print "SoapSplice JN                  : $s_j\n";

