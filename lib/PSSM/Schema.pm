package PSSM::Schema;
use base 'DBIx::Class::Schema';
use lib qw(../lib ..);
use strict;
use warnings;

 __PACKAGE__->load_classes( 
   { PSSM => [qw/Matrix MatrixPredictor Predictor NW/] }
  );
1;
