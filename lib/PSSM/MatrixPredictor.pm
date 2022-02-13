#$Id: MatrixPredictor.pm 192 2013-07-13 17:18:58Z maj $
package PSSM::MatrixPredictor;
use base qw/DBIx::Class::Core/;
use strict;
use warnings;

__PACKAGE__->table('mx_pred');
__PACKAGE__->add_columns(
  id => { data_type  => 'INT',
	  is_auto_increment => 1},
  matrix_id => { data_type => 'INT'},
  predictor_id => { data_type => 'INT'},
  coeff => { data_type => 'REAL',
	     default_value => 1.0 }
);

__PACKAGE__->set_primary_key('id');
__PACKAGE__->belongs_to('matrix' => 'PSSM::Matrix', 
			{'foreign.id' => 'self.matrix_id'});
__PACKAGE__->belongs_to('predictor' => 'PSSM::Predictor',
			{'foreign.id' => 'self.predictor_id'});

1;
