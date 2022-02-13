#$Id: Freq.pm 192 2013-07-13 17:18:58Z maj $
#####################
# PSSM::Freq: objects containing unconditional residue frequency vectors
#####################

package PSSM::Freq;
use lib '../../lib';
use base 'PSSM';
use strict;
use warnings;

sub new {
    my $class = shift;
    my %args = @_;
    my ($freqs, $type, $name) = @args{qw( FREQS TYPE NAME )};

    my $self = $class->SUPER::new();
    die "Frequency hash required arg" unless $freqs;
    
    $self->freqs($freqs);
    $self->residues( [sort keys %$freqs] );
    $self->seq_type($type);
    $self->name($name);

    return $self;
}

=head2 residues

 Title   : residues
 Usage   : $obj->residues($newval)
 Function: accesses arrayref of residue symbols
 Example : 
 Returns : value of residues (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub residues{
    my $self = shift;

    return $self->{'residues'} = shift if @_;
    return $self->{'residues'};
}

=head2 freqs

 Title   : freqs
 Usage   : $obj->freqs($newval)
 Function: access hashref of frequency table
 Example : 
 Returns : value of freqs (a hashref
 Args    : on set, new value (a hashref or undef, optional)

=cut

sub freqs{
    my $self = shift;

    return $self->{'freqs'} = shift if @_;
    return $self->{'freqs'};
}
  
=head2 seq_type

 Title   : seq_type
 Usage   : $obj->seq_type($newval)
 Function: accesses sequence type for these residues
 Example : 
 Returns : value of seq_type (a scalar) ('NT' or 'AA')
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub seq_type{
    my $self = shift;
    if (@_) {
	my $t = shift;
	die "Requires 'NT' or 'AA' for arg" unless $t =~ /^NT|AA$/;
	return $self->{'seq_type'} = $t;
    }
    return $self->{'seq_type'};
}

1;
