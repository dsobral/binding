use strict;
use warnings;

use Bio::BindingMatrix;

my $matrix_file = $ARGV[0];
my $seq_file = $ARGV[1];
my $threshold = 0;
if($ARGV[2]){ $threshold = $ARGV[2]};

die "Usage: BindingMatrix_MotifFeature.pl matrix_file sequence_file [threshold]" if(!$matrix_file && !$seq_file);

open(FILE,$matrix_file);
my @mats;
my $mat;
my $mat_name;
my $i = 1;
while(<FILE>){	
	my $line1 = $_;
	if($line1 =~ /^>(\S+)*/){
		$mat_name=$1;
		$line1 = <FILE>;	
	} else {
		$mat_name = "Matrix"+$i;
	}
	my $line2 = <FILE>;
	my $line3 = <FILE>;
	my $line4 = <FILE>;
	$mat = $line1.$line2.$line3.$line4;
	push(@mats,Bio::BindingMatrix->new($mat_name, $mat));
	$i++;
}
close FILE;

#example of a matrix
#my $mat = ">CTCF_Test  
#A [ 87	167	281	56	8	744	40	107	851	5	333	54	12	56	104	372	82	117	402  ]
#C [ 291	145	49	800	903	13	528	433	11	0	3	12	0	8	733	13	482	322	181]
#G [ 76	414	449	21	0	65	334	48	32	903	566	504	890	775	5	507	307	73	266 ] 
#T [459	187	134	36	2	91	11	324	18	3	9	341	8	71	67	17	37	396	59 ] ";

open(FILE,$seq_file);
my %sequences;
my $sequence;
my $seq_name;
while(<FILE>){
	chomp;
        my $line = $_;
	$line =~ s/\s+$//g;
	next if(!$line);
        if($line =~ /^>(\S+)*/){
		my $new_name = $1;
		if($sequence){ 
			$sequences{$seq_name} = $sequence; 
			$sequence = "";
		}
        	$seq_name=$new_name;
        } else {
		if(!$seq_name){ $seq_name="SEQUENCE"; }	
		if(!$sequence){ $sequence = $line;} else { $sequence .= $line; }
	}
}
close FILE;
#don't forget the last
$sequences{$seq_name} = $sequence;

#tests with CTCF.matrix and CTCF_seqs.fa as input
#print $mats[0]->name."\n";
#print $mats[1]->name."\n";
#numeric tests do not work with decimals?
#print $mats[0]->relative_affinity($sequences{"site1"},1)."\n";
#print $mats[0]->relative_affinity($sequences{"site2"},1)."\n";
#print $mats[0]->relative_affinity($sequences{"site3"},1)."\n";
#print $mats[0]->relative_affinity($sequences{"site4"},1)."\n";
#print $mats[0]->relative_affinity($sequences{"site5"},1)."\n";
#print $mats[0]->relative_affinity($sequences{"site6"},1)."\n";

foreach my $seq_name (keys %sequences){
	foreach my $matrix (@mats){
		for(my $i=0; $i<= (length($sequences{$seq_name})-$matrix->length); $i++){
			if($matrix->relative_affinity(substr($sequences{$seq_name},$i,$matrix->length),1) >= $threshold){
				print $seq_name."\t".($i+1)."\t".($i+$matrix->length);
				print "\t".$matrix->name."\t".substr($sequences{$seq_name},$i,$matrix->length);
				print "\t".$matrix->relative_affinity(substr($sequences{$seq_name},$i,$matrix->length),1);
				print "\n";	
			}
		}
		#my $len = $matrix->length
		#print $seq_name."\t".length($sequences{$seq_name})."\t".$matrix->name."\t".$matrix->length."\n";
	}
}


1;
