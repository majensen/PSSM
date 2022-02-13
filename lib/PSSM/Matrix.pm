# $Id: Matrix.pm 340 2014-02-02 01:23:07Z maj $ #
# pssm handling and building

package PSSM::Matrix;
use base qw/DBIx::Class::Core/;
use feature qw(switch);
use lib '../lib';
use PSSM::Ops;
use JSON;
use List::MoreUtils qw(all each_array);

use overload 
  '+' => \&Add,
  '*' => \&Cmult,
  '-' => \&Subtract,
  '==' => \&Eq,
  '&' => \&Conform,
  'bool' => \&Bool,
  'eq' => \&Eq;


use strict;
use warnings;

our $VERSION = '0.1';
our $SEQ_SYMBOLS  = 'a-zA-Z-~\.\?';
# our $err = new FB::ERR(-TITLE=>'spaTyper',-CSS=>"$HTDOCS/fb.css",-HREF=>"spaTyper.pl",-ICON=>"$HTDOCS/sT-icon.ico");

# build : initialize a new Matrix object, providing a matrix or a matrix 
# filename
# (new() is left to the parent dBIC class)
sub build{
   no warnings qw(experimental);
   my ($self, $mx) = @_;
   given (ref $mx) {
     when (/PSSM/) { $self->dta($mx->dta); }
     when (/HASH/) { $self->dta($mx); }
     when ('') {
       if (-e $mx) {
	 $self->load_matrix($mx);
       }
       else {
	 die "Matrix file error: ".$! ;
       }
     }
     default {
       die "new(): matrix arg type not supported";
     }
   }
   my $dta = $self->dta;
   if ($dta) {
     $self->residues( [sort keys %$dta] );
     $self->length( scalar @{$$dta{$self->residues->[0]}} );
   }
   return $self;
}

sub build_new {
  my $class = shift;
  my ($schema, @args) = @_;
  unless (ref $schema and $schema->isa('PSSM::Schema')) {
    die 'build_new requires a PSSM::Schema obj (arg 1)';
  }
  my $self = $schema->source('Matrix')->resultset->new_result({});
  return $self->build(@args);
}

=head2 simple_score

 Title   : simple_score
 Usage   : $score = $mx->simple_score(@seq)
 Function: score the input sequence using the matrix, without N-W alignment 
           to the matrix 
 Example :
 Returns : array
 Args    : [array of] string scalars

=cut

sub simple_score{
   my ($self,@seqs) = @_;
   my $mx = $self->dta or die "Matrix not yet initialized";
   my @ret;
   foreach my $s (@seqs) {
       my @s = split('',$s);
       
       push @ret, eval ( join('+',
			    map { 
				$$mx{$s[$_]}[$_] || 0 
			    } (0..$self->length-1)
			 ) );
       warn "in eval: $@" if $@;
   }
   return wantarray ? @ret : $ret[0];
}

=head2 as_text

 Title   : as_text
 Usage   : print $mx->as_text
 Function: return a tab-delimited version of the data matrix
 Example :
 Returns : array of text lines
 Args    : 'NOHDR' as arg 1 will suppress row/col headers

=cut

sub as_text{
    my $self = shift;
    my @args = @_;
    my @ret;
    my $nohdr = ($_[0] eq 'NOHDR');
    my $dta = $self->dta;
    unless ($nohdr) {
	push @ret, "\t", join("\t", (1..$self->length)),"\n";
    }
    foreach my $res (@{$self->residues}) {
	my $line = ($nohdr ? "" : "$res\t");
	$line .= join("\t", map {$dta->{$res}[$_]} (0..$self->length-1) )."\n";
	push @ret, $line;
    }
    return @ret;
}

=head1 Accessors

=over

=item dta()

Hashref of matrix. Keys are residues, values are arrayref of real values.

=item residues()

Arrayref of residues (keys of the matrix hash).

=item id()

DB id of matrix

=back

=head2 load_matrix

 Title   : load_matrix
 Usage   : called when new has a filename argument
 Function: load matrix from a file (as_text output works)
 Example :
 Returns : 
 Args    : matrix text filename

=cut

sub load_matrix { 
# get matrix from file

    my $self = shift;
    my ($fn) = @_;
    open my $mat, $fn || die "Matrix file problem: $!\n";
    my @infile = <$mat>;
    close $mat;

    my (%mx,$residues,$i);

    @infile = grep !/^#/, @infile; # get rid of comment lines
    foreach (@infile) {
	my(@a);
	@a = split /\s+/;
	if  (/^[${SEQ_SYMBOLS}]/) {
	    my($a) = shift(@a); 
	    # assume first string on line is one-letter 
	    # nt designator
	    push @{$residues}, $a;
	    $mx{$a} = \@a; # rest are site scores for that aa
	}
    }
    $self->dta( \%mx );
    $self->residues( [ sort @$residues ] );
    $self->length( scalar @{$mx{$$residues[0]}} );
    return 1;
}

sub store {
  my $self = shift;
  my ($name) = @_;
  unless ($name && !ref $name) {
    warn "Cannot store matrix without scalar name (arg 1)";
    return;
  }
  $self->name($name);
  $self->insert;
}

sub _schema {
  my $self = shift;
  return $self->result_source->schema;
}


=head2 consensus

 Title   : consensus
 Usage   : 
 Function: get consensus residues at each matrix site, according to 
           votes in the matrix, as a string
 Example :
 Returns : scalar string
 Args    : (1) : ignore gaps and X char

=cut

sub consensus{
   my ($self,$nogap) = @_;
   my $cons = "";
   my @res = $nogap ? grep(/[^X-]/, @{$self->residues} ) : @{$self->residues};
   my $dta = $self->dta;
   my @rows = map { $dta->{$_} } @res;
   foreach my $col (0..$self->length-1) {
       my @scr = map { $_->[$col] } @rows;
       my %scr;
       @scr{@res} = @scr;
       my @sres = sort {-$scr{$a} <=> -$scr{$b}} @res;
       $cons .= $sres[0];
   }
   return $cons;
}

sub mean_score_per_residue {
   my ($self,$nogap) = @_;
   my $acc = 0;
   my @res = $nogap ? grep(/[^X-]/, @{$self->residues} ) : @{$self->residues};
   my $dta = $self->dta;
   my @rows = map { $dta->{$_} } @res;
   foreach my $col (0..$self->length-1) {
       my @scr = map { defined $_ ? $_->[$col] : () } @rows;
       $acc += eval( join('+', @scr) )/@scr;
   }
   return $acc/$self->length;
}

sub _chk_pssm_type {
  my ($M) = @_;
  ref $M && all {$_} map {$M->can($_)} qw(dta residues length);
}

sub _chktypes  {
  return all {$_} map {_chk_pssm_type($_) } @_;
}

__PACKAGE__->table('mxs');
__PACKAGE__->add_columns(
  id => { data_type  => 'INT',
	  is_auto_increment => 1},
  name => { data_type => 'TEXT'},
  dta => { data_type => 'TEXT',
	   is_nullable => 1},
  residues => { data_type => 'TEXT',
		is_nullable => 1},
  length => { data_type => 'INT',
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

__PACKAGE__->inflate_column( 
  'dta', {
    inflate => sub {
      my ($mx_hash, $result_obj) = @_;
      $result_obj = decode_json $mx_hash;
    },
    deflate => sub {
      my ($mx_hash, $result_obj) = @_;
      encode_json $mx_hash;
    }
   }
);

__PACKAGE__->inflate_column( 
  'residues', {
    inflate => sub {
      my ($residues_array, $result_obj) = @_;
      $result_obj = decode_json $residues_array;
    },
    deflate => sub {
      my ($residues_array, $result_obj) = @_;
      encode_json $residues_array;
    }
   }
);

__PACKAGE__->has_many( 'mx_preds' => 'PSSM::MatrixPredictor',
		       { 'foreign.predictor_id' => 'self.id' } );
__PACKAGE__->many_to_many( 'predictors' => 'mx_preds', 'predictor' );


=head1 NAME

PSSM::Schema::Result::Matrix - DBIC source for PSSMs

=head1 SYNOPSIS

=head1 DESCRIPTION

Fields:

=over

=item id

Autoincrement integer primary key.

=item name

Unique name for the matrix.

=item dta

The matrix data. The form of the matrix is a hashref whose keys are the
residues in upper case, and whose values are arrayrefs of floats. The
arrayrefs are all ofthe same length; this is the length in residues of
the motif.

The class deflates the dta hashref to L<JSON> with L<JSON\encode_json>.
It inflates using L<JSON\decode_json>.

=item residues

A lexically sort array of residues, the keys of the dta hash.

The class deflates the residues arrayref to L<JSON> with L<JSON\encode_json>.
It inflates using L<JSON\decode_json>.

=item length

The length of the PSSM motif = the length of all dta arrayref values.

=item date_added, expires

=item deprecated_q

Boolean. TRUE if the matrix is 'deprecated'.

=back

=head1 AUTHOR

MAJ

=cut

1;


