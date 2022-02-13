#$Id: NW.pm 363 2014-02-26 03:41:41Z maj $
package PSSM::NW;
use base qw/DBIx::Class::Core/;
use PSSM::NW::alg;
use feature qw/switch/;
# codes for incremental path vectors
use constant { II => 4, IZ => 2, ZI => 1 };
our $SEQ_SYMBOLS  = 'a-zA-Z-~\.\?';
our $GAP_SYMBOLS = '\\~\\-\\.\\?';
our %OBJTABLE = ();
use lib '../../lib';
use strict;
use warnings;

# Bioperlize beginning 12/10/08
# objectify beginning 12/10/08
#
# matrix interface: a stripped-down version of MX? 
# - include a matrix getter from a regular text file
#
# original version aligns a nt seq to an nt matrix >by codon<
# how to incorporate LocatableSeq mapping?
#
# oldspan => newspan

# mapping: align
# nt to nt
# nt to aa
# aa to nt
# aa to aa
# mapping of [1 => 1] means seq in aas/nts, coords in aas/nts, resp.
# mapping of [1 => 3] means sequence in aas, coords in nt
# mapping of [3 => 1] means sequence in nt, coords in aa
# mapping of [3 => 3] means sequence in nt, coords in nt, but treat in 
#                     groups of 3 (as codons)

# so spatyper mapping is align seq to matrix under a 3 => 3 mapping
# mapping to an alignment might be referred to as: 
#   [ sequence coords => matrix coords ]
#   or [$sspan => $mspan]

#
# cgi version: MAJ 12/25/07
#
# 2/17/08
# want to inhibit (not completely) the addition of gaps
# into the PSSM sequence compared to the target, but allow
# the pssm to slide more freely across the target. I.e., 
# want gaps before/after the PSSM motif, more rarely gaps
# within the motif. 
# For aureus.mx, the best score possible is 33, average pssm 
# best incremental score per site is 4.125. Of the next highest
# scoring codons, the avg inc score is -5.875.
# 
#
# use NW-type algorithm to identify
# repeat elements of different lengths
# with pssm
# codon-based; accepts sequences with length divisible by 3.


# slide 8, right panel shows the path vector possibilities;
# more than one vector can occupy a cell, so the path vectors
# are coded in binary, such that ORing them will allow representation
# of multiple vectors

# intentional globals...


#to prevent addition of gaps to sequence $s, set $first_seq_fixed:
# (not required in PSSM use)

#default aa align
# set NW parameters; these could be tweaked
our $GXPEN = 0; # gap extension penalty 
our $GPEN = 2; # gap penalty
our $GPEN_FL = 0; # flanking gap penalty
our $GPEN_IN = 8; # internal gap penalty

# restrict gaps to occur at codon boundaries (w/r to pssm), 
# then insert a 3-nt gap.

# constructor

=head2 build_new(), build()

 my $nw = PSSM::NW->build_new($schema, %args); # create new obj
 $nw->build(%args) # set internals on preexisting obj

 Args: MATRIX => PSSM::Matrix obj (underlying matrix)
       FSFIXED => 1|0 (fix the first sequence--the matrix--preventing gaps)
       GPEN => $gap_penalty
       GXPEN => $gap_extension_penalty
       GPEN_FL => $flanking_gap_penalty
       GPEN_IN => $internal_gap_penalty
       MAPPING => [$matrix_residue_number => $sequence_residue_number]

=cut

sub build_new {
  my $class = shift;
  my ($schema, @args) = @_;
  unless (ref $schema && $schema->isa('PSSM::Schema')) {
    die "Requires a PSSM::Schema object (arg 1)";
  }
  return $schema->source('NW')->resultset->new_result({})->build(@args);
}


sub build {
    my $self = shift;
    my %nw_args = @_;
    my ($MX, $FSFIXED, $GPEN, $GXPEN, $GPEN_FL, $GPEN_IN, $MAPPING) =
	@nw_args{qw( MATRIX FSFIXED GPEN GXPEN GPEN_FL GPEN_IN MAPPING)};

    foreach my $parm (qw( MATRIX FSFIXED GPEN GXPEN GPEN_FL GPEN_IN MAPPING)) {
      no warnings qw(experimental);
      given ($nw_args{$parm}) {
	when ($parm eq 'MATRIX') {
	  unless (!$nw_args{$parm} || $nw_args{$parm}->name) {
	    die 'PSSM::Matrix must have name() field set';
	  }
	  $self->matrix($_);
	}
	when (defined) {
	  my $method = lc $parm;
	  $self->$method($_);
	}
	default {
	  1;
	}
      }
    }
    $self->update_or_insert; # store
    $self->discard_changes; # pick up the database default values
    my $id = $self->id;
    return $self;
}

## accessors

=over

=item matrix

 Title   : matrix
 Usage   : $obj->mx;
 Function: 
 Example : 
 Returns : associated L<PSSM::Matrix> object

=item fsfixed()

 Title   : fsfixed
 Usage   : $obj->first_fixed();
 Function: if true, first sequence is fixed (no gaps inserted);
           if false, gaps are allowed in first sequence
 Example : 
 Returns : value of first_fixed (a scalar)


=item gpen()

 Title   : gpen
 Usage   : $obj->gpen;
 Function: gap penalty
 Example : 
 Returns : value of gpen (a scalar)

=item gxpen()

 Title   : gxpen
 Usage   : $obj->gxpen;
 Function: gap traversal penalty
 Example : 
 Returns : value of gxpen (a scalar)

=item gpen_fl()

 Title   : gpen_fl
 Usage   : $obj->gpen_fl
 Function: flanking gap penalty
 Example : 
 Returns : value of gpen_fl (a scalar)

=item gpen_in()

 Title   : gpen_in
 Usage   : $obj->gpen_in
 Function: internal gap penalty
 Example : 
 Returns : value of gpen_in (a scalar)

=item mapping()

 Title   : mapping
 Usage   : $arrayref = $obj->mapping
 Function: sequence 1 to sequence 2 residue number mapping
 Example : 
 Returns : value of mapping (an arrayref of two ints)

=back

=head2 nw_parms

 Title   : nw_parms
 Usage   : read alignment parameters as hash
 Function:
 Example :
 Returns : 
 Args    :

=cut

sub nw_parms{
    my $self = shift;
    return (
	'multiplier' => $self->mapping->[1],
	'gap penalty' => $self->gpen,
	'gap traversal penalty' => $self->gxpen,
	'flanking gap penalty' => $self->gpen_fl,
	'internal gap penalty' => $self->gpen_in,
	'matrix anchored' => $self->fsfixed ? 'yes' : 'no'
	);
}

=head2 num_alns

 Title   : num_alns
 Usage   : $obj->num_alns($newval)
 Function: number of isoscoring alignments
 Example : 
 Returns : value of num_alns (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub num_alns{
    my $self = shift;
    return $OBJTABLE{$self->id}{'num_alns'} = shift if @_;
    return $OBJTABLE{$self->id}{'num_alns'};
}

=head2 seq

 Title   : seq
 Usage   : $obj->seq($newval)
 Function: the sequence string of the aligned sequence
           (set by method align())
 Example : 
 Returns : value of seq (a sequence object or scalar string)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub seq{
    my $self = shift;
    if (@_) {
	my $s = shift;
	die "Bad char in seq" if !ref($s) && $s =~ /[^$SEQ_SYMBOLS]/;
	return $OBJTABLE{$self->id}{'seq'} = $s;
    }
    return  $OBJTABLE{$self->id}{'seq'};
}

# sequence string only
sub seqstr{
    my $self = shift;
    my $s = $self->seq;
    return $s->seq if ref($s);
    return $s;
}

=head2 align

 Title   : align
 Usage   :
 Function: align sequence arg to the matrix and set properties 
 Example :
 Returns : 
 Args    : a sequence

=cut

sub align { #align a seq to a pssm by Needleman-Wunsch
    my $self = shift;
    my ($s) = @_; # sequence
    ref($s) =~ 'Bio' && do {
	if ($s->can('seq')) {
	    $self->seq( $s );
	}
	else {
	    die "Require sequence or seq object in arg";
	}
    };
    !ref($s) && do {
	# validation in seq() method
	$self->seq($s);
    };
    $s = uc $self->seqstr; # up case
    $s =~ s/[$GAP_SYMBOLS]//g; # degap

    die "Matrix not yet set in object" unless defined $self->matrix;
    my ($mx, $dist, $pthm, $scr);
    $mx = $self->matrix->dta; # the raw matrix
    my $ws = $OBJTABLE{$self->id};
    for (qw/aln_score score paths num_alns ub lb/) { #reset
      undef $ws->{$_};
    }
    
    my $mlen = $self->matrix->length();
    my $parms = [ $mlen, $self->gpen_in, $self->gpen_fl, 0,
		 @{$self->mapping}];
    my @paths;
    ($ws->{aln_score},$ws->{score}, @paths) = nw_align($s, $mx, $parms);
    $ws->{paths} = \@paths;
    $ws->{num_alns} = scalar @paths;
    foreach my $pth (@paths) {
      my @wk = @$pth;
      my $i = 0;
      my $lb = my $ub = -1;
      while (my $c = pop @wk) {
	if ($c == 2) {
	  $i++;
	}
	if ($c == 3) {
	  $lb = ($i*$self->mapping->[0])+1 if $lb < 0;
	  $ub = (($i+1)*$self->mapping->[0]);
	  $i++;
	}
      }
      push @{$ws->{lb}},$lb;
      push @{$ws->{ub}},$ub;
    }
    return 1;
}

sub aln_score { $OBJTABLE{shift()->id}{aln_score} }
sub score { $OBJTABLE{shift()->id}{score} }
sub paths { $OBJTABLE{shift()->id}{paths} }
sub ub { $OBJTABLE{shift()->id}{ub} }
sub lb { $OBJTABLE{shift()->id}{lb} }


=head2 get_align

 Title   : get_align
 Usage   :
 Function: return a gap-inserted sequence and sequence alignment boundaries
 Example :
 Returns : (gap-inserted sequence, lower bound, upper bound) array of scalars
 Args    : ($seq, $idx) : $idx specifies the index of the alignment desired
           from among the isoscoring alignments (see method num_alns() )
           ($seq) : returns seqs from the next isoscoring alignment in 
           the queue (when all sequences have been retrieved, this form
           returns () )

=cut

sub get_align { 
# return a gap-inserted seq, using string of bases and binary gap template, last arg opt. num bases
    my $self = shift;
    my $idx = shift;
    if (!defined $idx) {
	$idx = each @{$self->paths};
	return unless defined $idx; # done with queue
    }
    my $seq = $self->seqstr || die "No sequence aligned yet";
    my $ws = $OBJTABLE{$self->id};
    my $pth = $self->paths->[$idx];
    # copy pth
    my ($sq, $mx, $lb, $ub) = apply_path([@$pth], $seq, $self->matrix->consensus,@{$self->mapping});
    return { idx => $idx,
	     seq => $sq,
	     mxc => $mx,
	     lb  => $lb,
	     ub  => $ub };
}

## internal accessors

=head2 _destroy_array

 Title   : _destroy_array
 Usage   :
 Function: unlink all deep references
 Example :
 Returns : 
 Args    :

=cut

sub _destroy_array{
    my $a = shift;
    return 1 unless ref($a) eq 'ARRAY';
    unless (ref($a->[0])) {
	$a = undef;
	return 1;
    }
    my @a = @$a;
    foreach my $i (@a) {
	_destroy_array($i);
    }
    return 1;
}

sub _clone_array{
    my $a = shift;
    return () unless ref($a) eq 'ARRAY';
    my $ret = [];
    unless (ref($a->[0])) {
	return [@$a];
    }
    my @a = @$a;
    foreach my $i (@a) {
	push @$ret, _clone_array($i);
    }
    return $ret;
}


__PACKAGE__->table('nwxs');
__PACKAGE__->add_columns(
  id => { data_type => 'INT',
	  is_auto_increment => 1 },
  name => { data_type => 'TEXT',
	    is_nullable => 1},
  matrix => { data_type => 'INT',
	      is_nullable=> 1,
	      is_foreign_key => 1 },
  gpen => { data_type => 'REAL', 
	    default_value => 2.0 },
  gpen_fl => { data_type => 'REAL',
	       default_value => 0.0 },
  gxpen => { data_type => 'REAL',
	     default_value =>  0.0 },
  gpen_in => { data_type => 'REAL', 
	       default_value => 8.0 },
  fsfixed => { data_type => 'INT', 
	       default_value => 1 },
  mapping => { data_type => 'TEXT',
	       default_value => '1:1' }
);

__PACKAGE__->add_unique_constraint(['name']);
__PACKAGE__->set_primary_key('id');
__PACKAGE__->inflate_column(
  mapping => {
    inflate => sub {
      my ($mapstr, $result_obj) = @_;
      $result_obj = [split /:/, $mapstr];
    },
    deflate => sub {
      my ($map, $result_obj) = @_;
      join(':',@$map);
    }
   }
);
__PACKAGE__->belongs_to( matrix => 'PSSM::Matrix' );

1;
