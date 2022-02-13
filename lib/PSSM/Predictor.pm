#$Id: Predictor.pm 366 2014-03-04 01:39:06Z maj $
package PSSM::Predictor;
use base qw/DBIx::Class::Core/;
#use base qw/PSSM::Schema::Result::Predictor PSSM/;
use lib '../../lib';
use JSON;
use List::MoreUtils qw(pairwise all);
use strict;
use warnings;

our $OBJTABLE = {};

=head2 build

 Title   : build
 Usage   : $pred->build(@args)
 Function: build (load) a predictor object as a linear combination of PSSM
 Example : create a predictor with a single PSSM matrix: 
           $pred->build( COEFFS=>[1], MXS=>[$M] );
 Returns : a PSSM::Predictor object
 Args    : key=>value pairs; 
           CUTOFF can be set to a real scalar (the prediction cutoff score)
           CALLBACK can be set to an anonymous sub {} that will
           be passed  the predictor and sequence, and should 
           return 0 or 1

=cut

sub build {
   my ($self,@args) = @_;
#   return $class->new_from_obj($args[0]) if $args[0] and ref($args[0]) =~ /Predictor/;
   my %args = @args;
   my ($coeffs, $mxs, $cutoff, $plbl, $scr95, $pred_callback) = 
       @args{qw( COEFFS MXS CUTOFF PLBL SCR95 CALLBACK)};
   #init
   $coeffs ||= [];
   $mxs ||= [];
   $cutoff && ($cutoff eq 'none') && ($cutoff=undef);

   #sanity
   die "Require equal numbers of coefficients and matrices" if ($coeffs && $mxs && scalar @$coeffs != scalar @$mxs);
   die "Requires a list of PSSM matrices for -MXS arg" 
       if ( $mxs && grep { ref($_) !~ /^PSSM/ } @$mxs );
   die "All matrices must conform in dimension" if ($coeffs && $mxs && !conform_ck($mxs));
   die "All matrices must be stored in DB (use \$mx->store(\$name))" unless all { defined $_->id } @$mxs;
   if ($pred_callback) {
       die "Predictor callback not a CODE ref" unless ref($pred_callback) eq 'CODE';
       $OBJTABLE->{$self->id}{_callback} = $pred_callback;
   }

   $self->cutoff($cutoff) if defined $cutoff;
   $self->plbl($plbl) if defined $plbl;
   $self->scr95($scr95) if defined $scr95;
   my $dta;
   foreach my $res (@{$$mxs[0]->residues}) {
     $dta->{$res} = [ (0) x $$mxs[0]->length ];
   }
   my $mx = $self->result_source->schema->source('Matrix')->resultset->new_result({})->build($dta);
   pairwise {
     $self->add_to_matrices($b, {coeff => $a});
     $mx = $mx + $a*$b;
   } (@$coeffs, @$mxs);
   
   $OBJTABLE->{$self->id}{_mx} = $mx;
   return $self;
}


sub build_new {
  my $class = shift;
  my ($schema, @args) = @_;
  unless (ref $schema and $schema->isa('PSSM::Schema')) {
    die 'build_new requires a PSSM::Schema obj (arg 1)';
  }
  my $self = $schema->source('Predictor')->resultset->new_result({});
  $self->build(@args);
}
sub new_from_obj {
  my $class = shift;
  my ($obj) = @_;
  die "new_from_obj requires a Predictor obj (arg 1)" unless (
    $obj->isa('PSSM::Predictor') || 
      $obj->isa('PSSM::Schema::Result::Predictor'));
  my $dta;
  my @mxs = $obj->matrices;
  my @coeffs;
  foreach my $mx (@mxs) {
    my ($item) = $obj->matrix_in_predictor({matrix_id => $mx->id});
    push @coeffs, $item->coeff;
  }
  foreach my $res (@{$mxs[0]->residues}) {
    $dta->{$res} = [ (0) x $mxs[0]->length ];
  }
  my $mx = PSSM::Matrix->build_new($obj->result_source->schema,$dta);
  pairwise {
    $mx = $mx + $a*$b;
  } (@coeffs, @mxs);
  
  $OBJTABLE->{$obj->id}{_mx} = $mx;
  return $obj;
}

=head2 score

 Title   : score
 Usage   : $pred->score(@seqs);
 Function: score a Bio::SimpleAlign/[array of] Bio::PrimarySeq/[array of]
           seqstr using the predictor, also returns a prediction array 
           if $pred->cutoff is set
 Example :
 Returns : arrayref to scores, arrayref to predictions
 Args    : a [collection of] sequence[s]

=cut

sub score {
    my $self = shift;
    my $seqs = shift;
    die "Predictor not initialized"    if (!defined $self->mx);
    return () unless $seqs;
    my @seqs;
    if (!ref($seqs)) {
	push @seqs, $seqs;
    }
    elsif ((ref($seqs) eq 'ARRAY')) {
	@seqs = ref($seqs[0]) ? @$seqs : map {$_->seq} @$seqs;
    }
    elsif ($seqs->isa('Bio::SimpleAlign')) {
	@seqs = map {$_->seq} $seqs->each_seq;
    }
    elsif ($seqs->isa('Bio::PrimarySeq')) {
	push @seqs, $seqs->seq;
    }
    else {
	die "Can't understand arg";
    }
    
    my @scr = $self->mx->simple_score(@seqs);
    my @pred;
    if (defined $OBJTABLE->{$self->id}{_callback}) {
	@pred = $self->_callback(\@seqs, \@scr);
    }
    else {
	@pred = defined $self->cutoff ? (map {$_ >= $self->cutoff ? 1 : 0} @scr) : ();
    }
    return [@scr], [@pred];
}

=head2 predict

 Title   : predict
 Usage   : @predictions = $pred->predict(@seqs)
 Function: score each sequence and return 1 if score > cutoff, 
           0 if score < cutoff
 Example :
 Returns : array 
 Args    : [an array of] sequence[s]

=cut

sub predict{
   my ($self,@seqs) = @_;
   die "Predictor not initialized"    if (!defined $self->mx);
   if (!defined $self->{_callback}) {
       die "Cutoff not set" unless defined $self->cutoff;
       return map { $self->score($_) >= $self->cutoff ? 1 : 0 } @seqs;
   }
   else {
       return $self->_callback(\@seqs, $self->mx->simple_score(@seqs));
   }
}

=head2 conform_ck

 Title   : conform_ck
 Usage   : conform_ck(\@mxs)
 Function: check that all mxs in list have same dimensions
 Example :
 Returns : 1 if all mxs conform, 0 otherwise
 Args    :

=cut

sub conform_ck{
    my @mxs = @_;
    @mxs = @{$mxs[0]} if ref($mxs[0]) eq 'ARRAY';
    return (undef) unless @mxs;
    return (1) if scalar @mxs == 1;
    my $mx = $mxs[0];
    foreach (@mxs[(1..$#mxs)]) {
	return 0 unless ($mx & $_);
    }
    return 1;
}

sub simple_score { shift->mx->simple_score(@_) }

=head2 cutoff

 Title   : cutoff
 Usage   : $obj->cutoff($newval)
 Function: access the cutoff used for prediction
 Example : 
 Returns : value of cutoff (a scalar)
 Args    : on set, new value (a scalar or undef, optional) (ok to provide 
           arrayref of scalar cutoffs for use in callback predictor
           function)

=cut

=head2 id

 Title   : id
 Usage   : $pred->id($newval)
 Function: gets/sets predictor db id
 Example : 
 Returns : value of id (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

=head2 name

 Title   : name
 Usage   : $obj->name($newval)
 Function: gets/sets predictor name
 Example : 
 Returns : value of name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

=head2 mx

 Returns the PSSM of the predictor.

=cut

sub mx {
  my $self = shift;
  return $OBJTABLE->{$self->id}{_mx} if $OBJTABLE->{$self->id}{_mx};
  $self->new_from_obj($self); # sets the PSSM of the predictor;
  return $OBJTABLE->{$self->id}{_mx};
}

sub length { shift->mx->length }
sub residues { shift->mx->residues }
sub dta { shift->mx->dta }


=head2 _callback

 Title   : _callback
 Usage   : $pred->_callback(@args)
 Function: do a custom prediction based on this predictor
           ( not accessed directly by user, but called via
            $pred->score if callback ref was set in new() )
 Example : 
 Returns : return value of &{$pred->{_callback}}
 Args    : arguments passed through to coderef

=cut

sub set_callback {
  my $self = shift;
  my ($code) = (@_);
  die 'Arg 1 must be a coderef' unless (ref $code eq 'CODE');
  $OBJTABLE->{$self->id}{_callback}=$code;
  return 1;
}

sub _callback{
    my $self = shift;
    return &{$OBJTABLE->{$self->id}{_callback}}($self,@_);
}

=head2 store

 Title   : store
 Usage   : $pred->store($name)
 Function: persist pred and associated mxs in db
 Example :
 Returns : 
 Args    :

=cut

sub store {
  my $self = shift;
  my ($name) = @_;
  unless ($name && !ref $name) {
    warn "Cannot store predictor without scalar name (arg 1)";
    return;
  }
  $self->name($name);
  $self->insert;
  $self->update;
}


__PACKAGE__->table('preds');
__PACKAGE__->add_columns(
  id => { data_type  => 'INT',
	  is_auto_increment => 1},
  name => { data_type => 'TEXT',
	   is_nullable => 1},
  label => {data_type => 'TEXT',
	    is_nullable => 1},
  cutoff => { data_type => 'TEXT',
	      is_nullable => 1 },
  scr95 => {data_type => 'TEXT',
	    is_nullable => 1},
  plbl => {data_type => 'TEXT',
	   is_nullable => 1},
  type => { data_type => 'TEXT',
	    is_nullable => 1},
  date_added => { data_type => 'TEXT',
		  is_nullable => 1},
  expires => { data_type => 'TEXT',
	       is_nullable => 1},
  is_deprecated => { data_type => 'INT',
		     is_nullable => 1}
);

__PACKAGE__->add_unique_constraint(['name']);

__PACKAGE__->set_primary_key('id');

__PACKAGE__->has_many('matrix_in_predictor' => 'PSSM::MatrixPredictor',
		     { 'foreign.predictor_id' => 'self.id'} );
__PACKAGE__->many_to_many('matrices' => 'matrix_in_predictor', 'matrix');

__PACKAGE__->inflate_column(
  'cutoff', {
    inflate => sub {
      my ($intvl, $result_obj) = @_;
      $result_obj = decode_json $intvl;
    },
    deflate => sub {
      my ($intvl, $result_obj) = @_;
      encode_json $intvl;
    }
   });
__PACKAGE__->inflate_column(
  'scr95', {
    inflate => sub {
      my ($intvl, $result_obj) = @_;
      $result_obj = decode_json $intvl;
    },
    deflate => sub {
      my ($intvl, $result_obj) = @_;
      encode_json $intvl;
    }
   });
__PACKAGE__->inflate_column(
  'plbl', {
    inflate => sub {
      my ($lbl_ary, $result_obj) = @_;
      $result_obj = decode_json $lbl_ary;
    },
    deflate => sub {
      my ($lbl_ary, $result_obj) = @_;
      encode_json $lbl_ary;
    }
   });

1;
