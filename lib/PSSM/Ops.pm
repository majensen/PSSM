#$Id: Ops.pm 193 2013-07-14 03:03:28Z maj $
package PSSM::Matrix;
use List::MoreUtils qw/each_arrayref/;
use strict;
use warnings;


=head1 Matrix Operations

=head2 Add

 Title   : add
 Usage   : $P = add_mx($M, $N)
 Function: add $M->{_mx} and $N->{_mx} elt by elt, return 
           resulting op'd matrix
 Example :
 Returns : 
 Args    : two PSSM objects

=cut

sub Add {
    my ($M, $N) = @_;
    die "Requires two PSSM::Matrix objects"
	unless _chktypes($M, $N);

    die "Matrix dimensions do not conform" 
	unless $M & $N;
    my $m = $M->dta;
    my $n = $N->dta;
    my $mx = { 
	map {
	    my $res = $_;
	    ($res => [ 
		 map {
		     $m->{$res}->[$_] + $n->{$res}->[$_]
		 } (0..$M->length-1)
	     ])
	} @{$M->residues} };

    return $M->result_source->resultset->new_result({})->build($mx);
}

sub Subtract {
    my ($M, $N) = @_;
    die "Requires two PSSM::Matrix objects" unless _chktypes($M, $N);
    die "Matrix dimensions do not conform" unless $M & $N;
    return ($M + (-1)*$N);
}

sub Cmult {
    my ($c, $M) = @_;
    unless ($c && $M) { die "Requires a scalar and a PSSM object"; }
    if ( !ref($M) && ref($c) ) {
	my $tmp = $c;
	$c = $M;
	$M=$tmp;
    }
    die "Requires a scalar and a PSSM object" 
	unless (ref(\$c) eq 'SCALAR' && _chk_pssm_type($M));
    my $m = $M->dta;
    my $mx = { 
	map {
	    my $res = $_;
	    ($res => [ 
		 map {
		    $c * $m->{$res}->[$_]
		 } (0..$M->length-1)
	     ])
	} @{$M->residues} };
    return $M->result_source->resultset->new_result({})->build($mx);
}

sub Eq {
    my ($M, $N) = @_;

    die "Requires two PSSM objects"
	unless _chktypes($M,$N);
    unless ($M & $N) {
	warn "Matrix dimensions do not conform";
	return 0;
    }
    my $eq = 1;
    my $m = $M->dta;
    my $n = $N->dta;
  EQ:
    foreach my $res (@{$M->residues}) {
	foreach (0..$M->length-1) {
	    last EQ unless 
		$eq &&= ($m->{$res}->[$_] == $n->{$res}->[$_]);
	}
    }
    return $eq ;
}

sub Bool {
    my ($M) = @_;
    return ( defined $M );
}

sub Conform {
    my ($M, $N) = @_;
    die "Requires two PSSM objects"
	unless _chktypes($M,$N);
    return unless $M->length == $N->length;
    my @m = @{$M->residues};
    my @n = @{$N->residues};
    my $ea = each_arrayref( $M->residues, $N->residues );
    while (my ($m, $n) = $ea->()) {
      return unless $m eq $n;
    }
    my ($m, $n) = ($M->dta, $N->dta);
    foreach my $res (@{$M->residues}) {
      return unless @{$m->{$res}} == @{$n->{$res}};
    }
    return 1;
}

1;
