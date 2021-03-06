#
# module for Bio::BindingMatrix
#

=head1 LICENSE

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  dsobral@igc.gulbenkian.pt

=head1 NAME

Bio::BindingMatrix - A module to represent a BindingMatrix. 
In EFG this represents the binding affinities of a Transcription Factor to DNA.

=head1 SYNOPSIS

use Bio::BindingMatrix;

my $matrix = Bio::BindingMatrix->new(
                                     -name  => "MA0122.1",
                                     );
$matrix->frequencies("A  [ 4  1 13 24  0  0  6  4  9 ]
                      C  [ 7  4  1  0  0  0  0  6  7 ]
                      G  [ 4  5  7  0 24  0 18 12  5 ]
                      T  [ 9 14  3  0  0 24  0  2  3 ]");

print $matrix->relative_affinity("TGGCCACCA")."\n";

print $matrix->threshold."\n";

=head1 DESCRIPTION

This class represents information about a BindingMatrix, containing the name 
(e.g. the Jaspar ID, or an internal name), and description. A BindingMatrix 
is always associated to an Analysis (indicating the origin of the matrix e.g. 
Jaspar) and a FeatureType (the binding factor).   

=cut


use strict;
use warnings;

package Bio::BindingMatrix;

=head2 new

  Arg [-name]: string - name of Matrix
  Arg [-frequencies]: (optional) string - frequencies representing the binding affinities of a Matrix
  Arg [-threshold]: (optional) float - minimum relative affinity for binding sites of this matrix
  Example    : my $matrix = Bio::BindingMatrix->new(
                                                   "MA0122.1",
                                                   "Jaspar Matrix",
                                                );
  Description: Constructor method for BindingMatrix class
  Returntype : Bio::BindingMatrix

=cut

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;

  my $self = bless({}, $class);

  my ( $name, $freq, $thresh ) = @_;
  
  if(! defined $name){
    print("Must supply a name\n");
    exit(1);
  }

  $self->name($name);
  $self->frequencies($freq) if $freq;
  $self->threshold($thresh) if $thresh;

  return $self;
}



=head2 name

  Arg [1]    : (optional) string - name
  Example    : my $name = $matrix->name();
  Description: Getter and setter of name attribute
  Returntype : string

=cut

sub name {
    my $self = shift;
    $self->{'name'} = shift if @_;
    return $self->{'name'};
}


=head2 threshold

  Arg [1]    : (optional) float - threshold
  Example    : my $thresh = $matrix->threshold();
  Description: Getter and setter of threshold attribute 
  Returntype : float

=cut

sub threshold {
    my $self = shift;
    $self->{'threshold'} = shift if @_;
    return $self->{'threshold'};
}




=head2 frequencies

  Arg [1]    : (optional) string - frequencies
  Example    : $matrix->frequencies($frequencies_string);
  Description: Getter and setter of frequencies attribute

  The attribute is a string representing the matrix binding
  affinities in the Jaspar format. E.g. 
  	">
  	[ ]
  	"
  
  Returntype : string

=cut

sub frequencies {
  my $self = shift;
 
  my $frequencies = shift if @_; 
  if($frequencies){
  	$self->_weights($frequencies);
  	$self->{'frequencies'} = $frequencies;  	
  }
  return $self->{'frequencies'};
}

=head2 frequencies_revcomp

  Example    : $matrix->frequencies_revcomp();
  Description: Getter for the reverse complement frequencies attribute

  The attribute represents the reverse complement of frequencies
  
  Returntype : string

=cut

sub frequencies_revcomp {
  my $self = shift;
 
  return $self->{'frequencies_revcomp'};
}


=head2 relative_affinity

  Arg [1]    : string - Binding Site Sequence
  Arg [2]    : (optional) boolean - 1 if results are to be in linear scale (default is log scale)
  Example    : $matrix->relative_affinity($sequence);
  Description: Calculates the binding affinity of a given sequence
	relative to the optimal site for the matrix
	The site is taken as if it were in the proper orientation
        Considers a purely random background p(A)=p(C)=p(G)=p(T)
  Returntype : double

=cut

sub relative_affinity {
  my ($self, $sequence, $linear) = (shift, shift, shift);
  $sequence =~ s/^\s+//;
  $sequence =~ s/\s+$//;
  
  if(!$sequence) { print "No sequence given\n"; exit(1); }
  $sequence = uc($sequence);
  if($sequence =~ /[^ACGT]/){
    print "Sequence $sequence contains invalid characters: Only Aa Cc Gg Tt accepted\n";
    exit(1);
  }
  
  my $weight_matrix = $self->_weights;
  my $matrix_length = scalar(@{$weight_matrix->{'A'}});
  if(length($sequence) != $matrix_length){
    print "Sequence $sequence does not have length $matrix_length\n";
    exit(1);
  }
  
  my $log_odds = 0;
  my @bases = split(//,$sequence);
  for(my $i=0;$i<$matrix_length;$i++){
    $log_odds += $weight_matrix->{$bases[$i]}->[$i];	
  }
  
  #This log scale may be quite unrealistic... but usefull just for comparisons...
  if(!$linear){
    return ($log_odds - $self->_min_bind) / ($self->_max_bind - $self->_min_bind);
  } else {
    return (exp($log_odds) - exp($self->_min_bind)) / (exp($self->_max_bind) - exp($self->_min_bind));
  }

}

=head2 is_position_informative

  Arg [1]    : Int - 1bp position within the matrix
  Arg [2]    : (optional) double - threshold [0-2] for information content [default is 1.5]
  Example    : $matrix->is_position_informative($pos);
  Description: Returns true if position information content is over threshold
  Returntype : boolean

=cut

sub is_position_informative {
  my ($self, $position, $threshold) = @_;
  
  if(!defined($position)) {  print "Need a position\n"; exit(1); }
  if(($position<1) || ($position > $self->length)) { print "Position out of bounds\n"; exit(1); }
  if(!defined($threshold)){ $threshold = 1.5; }
  if(($threshold<0) || ($threshold>2)){ print "Threshold out of bounds\n"; exit(1); }
  return ($self->{'ic'}->[$position-1] >= $threshold);
}



=head2 length

  Example    : $bm->length();
  Description: Returns the length of the the matrix (e.g. 19bp long)
  Returntype : int with the length of this binding matrix 
  Exceptions : none

=cut

sub length {
  my $self = shift;

  my $weight_matrix = $self->_weights;

  return scalar(@{$weight_matrix->{'A'}});
}

=head2 _weights

  Arg [1]    : (optional) string - frequencies
  Example    : _weights($frequencies);
  Description: Private Getter Setter for the weight matrix based on frequencies
  Returntype : HASHREF with the weights of this binding matrix 

=cut

sub _weights {
	my $self = shift;
	
 	#for the moment use equiprobability and constant pseudo-count
	my $pseudo = 0.1;

 	#TODO allow for it to be passed as parameters?
  	my $frequencies = shift if @_; 
  	if($frequencies){
  		$frequencies =~ s/^(>.*?\n)//;
		my $header = $1;

  		my ($a,$c,$g,$t) = split(/\n/,$frequencies);
  		my @As = split(/\s+/,_parse_matrix_line('[A\[\]]',$a));
  		my @Cs = split(/\s+/,_parse_matrix_line('[C\[\]]',$c));
  		my @Gs = split(/\s+/,_parse_matrix_line('[G\[\]]',$g));
  		my @Ts = split(/\s+/,_parse_matrix_line('[T\[\]]',$t));
		if((scalar(@As)!=scalar(@Cs)) || (scalar(@As)!=scalar(@Gs)) || (scalar(@As)!=scalar(@Ts)) ){
			print "Frequencies provided are not a valid frequency matrix\n";
			exit(1);
		}
		$self->_calc_ic(\@As,\@Cs,\@Gs,\@Ts,$pseudo);

		#Create the reverse complement
		my @revT = reverse(@As);
		my @revA = reverse(@Ts);
		my @revC = reverse(@Gs);
		my @revG = reverse(@Cs);
		my $revcomp = $header;
		$revcomp.= "A [ ".join("\t",@revA)." ]\n";
		$revcomp.= "C [ ".join("\t",@revC)." ]\n";
		$revcomp.= "G [ ".join("\t",@revG)." ]\n";
		$revcomp.= "T [ ".join("\t",@revT)." ]\n";
		$self->{'frequencies_revcomp'} = $revcomp;

  		my @totals;
  		for(my $i=0;$i<scalar(@As);$i++){ 
  			$totals[$i]=$As[$i]+$Cs[$i]+$Gs[$i]+$Ts[$i];
  		}
  		
		my %weights;			
		#We can allow distinct background per nucleotide, instead of 0.25 for all... pass as parameter
		#But if the matrix was obtained using in-vivo data, it shouldn't matter the organism nucleotide bias..
		#We're using 0.1 as pseudo-count... the matrix cannot have very few elements... (e.g. <30 not good)
		my @was; for(my $i=0;$i<scalar(@As);$i++){ $was[$i] = log((($As[$i] + $pseudo) / ($totals[$i]+(4*$pseudo))) / 0.25); };
		$weights{'A'} = \@was;
		my @wcs; for(my $i=0;$i<scalar(@Cs);$i++){ $wcs[$i] = log((($Cs[$i] + $pseudo) / ($totals[$i]+(4*$pseudo))) / 0.25); };
		$weights{'C'} = \@wcs;
		my @wgs; for(my $i=0;$i<scalar(@Gs);$i++){ $wgs[$i] = log((($Gs[$i] + $pseudo) / ($totals[$i]+(4*$pseudo))) / 0.25); };
		$weights{'G'} = \@wgs;
		my @wts; for(my $i=0;$i<scalar(@Ts);$i++){ $wts[$i] = log((($Ts[$i] + $pseudo) / ($totals[$i]+(4*$pseudo))) / 0.25); };	     
		$weights{'T'} = \@wts;
	
		$self->{'weights'} = \%weights;

		my $max = 0; my $min = 0;
		for(my $i=0;$i<scalar(@As);$i++){
			my $col = [ $was[$i], $wcs[$i], $wgs[$i], $wts[$i] ];
			$min += _min($col);
			$max += _max($col);
		}

		#Log scale
		$self->_max_bind($max);
		$self->_min_bind($min);
	}
	
	return $self->{'weights'};	

}

=head2 _calc_ic

  Example    : _calc_ic($as,$cs,$gs,$ts,$pseudo);
  Description: Private function to calculate the matrix information content per position

=cut

sub _calc_ic {
  my ($self,$as, $cs, $gs, $ts,$pseudo) = (shift,shift, shift, shift, shift, shift);
  my @ic = ();
  for (my $i=0;$i<scalar(@$as);$i++){
    my $total_i = $as->[$i] + $cs->[$i] + $gs->[$i] + $ts->[$i] + (4*$pseudo);
    my $fas = ($as->[$i] + $pseudo) / $total_i;
    my $fcs = ($cs->[$i] + $pseudo) / $total_i;
    my $fgs = ($gs->[$i] + $pseudo) / $total_i;
    my $fts = ($ts->[$i] + $pseudo) / $total_i;    
    my $ic_i = 2 + ($fas * log($fas)/log(2)) + ($fcs * log($fcs)/log(2)) + ($fgs * log($fgs)/log(2)) + ($fts * log($fts)/log(2));
    push @ic, $ic_i;
  }
  $self->{'ic'} = \@ic;
}

sub _parse_matrix_line {
	 my ($pat,$line) = (shift,shift);
	 $line=~s/$pat//g; $line=~s/^\s+//; $line=~s/\s+$//;	
	 return $line;
}

sub _max { return _min_max(shift, 0); }

sub _min { return _min_max(shift, 1); }

sub _min_max {
	my ($list,$min) = (shift, shift);
	my $min_max = $list->[0];
	map { if($min ? $_ < $min_max : $_ > $min_max){ $min_max = $_; } } @$list;
	return $min_max;
}


=head2 _max_bind

  Arg [1]    : (optional) double - maximum binding affinity
  Example    : $matrix->_max_bind(10.2);
  Description: Private Getter and setter of max_bind attribute (not to be called directly)
  Returntype : float with the maximum binding affinity of the matrix 

=cut

sub _max_bind {
  my $self = shift;
  
  $self->{'max_bind'} = shift if @_;

  return $self->{'max_bind'};
}

=head2 _min_bind

  Arg [1]    : (optional) double - minimum binding affinity
  Example    : $matrix->_min_bind(-10.2);
  Description: Private Getter and setter of min_bind attribute (not to be called directly)
  Returntype : float with the minimum binding affinity of the matrix 

=cut

sub _min_bind {
  my $self = shift;
  
  $self->{'min_bind'} = shift if @_;

  return $self->{'min_bind'};
}

1;
