package PSSM;
use strict;
use warnings;

our $VERSION = '0.21';

sub get_item_by_name {
  my $class = shift;
  my ($schema, $item, $name) = @_;
  unless ($schema && ref $schema && $schema->isa('PSSM::Schema')) {
    die 'Schema must be provided (arg 1)';
  }
  unless ($name && !ref $name) {
    die 'Name must be specified (arg 2)';
  }
  return $schema->source($item)->resultset->find( { name => $name } );
}

sub get_item_by_id {
  my $class = shift;
  my ($schema, $item, $id) = @_;
  unless ($schema && ref $schema && $schema->isa('PSSM::Schema')) {
    die 'Schema must be provided (arg 1)';
  }
  unless ($id && !ref $id) {
    die 'ID must be specified (arg 2)';
  }
  return $schema->source($item)->resultset->find( { id => $id } );
}

sub get_matrix_by_name { shift->get_item_by_name($_[0],'Matrix',$_[1]) }
sub get_predictor_by_name { shift->get_item_by_name($_[0],'Predictor',$_[1]) }
sub get_nw_by_name { shift->get_item_by_name($_[0],'NW',$_[1]) }

sub get_matrix_by_id { shift->get_item_by_id($_[0],'Matrix',$_[1]) }
sub get_predictor_by_id { shift->get_item_by_id($_[0],'Predictor',$_[1]) }
sub get_nw_by_id { shift->get_item_by_id
($_[0],'NW',$_[1]) }

1;
