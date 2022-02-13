#$Id: 001_db.t 205 2013-07-29 04:07:41Z maj $
use Test::More;
use lib '../lib';
use PSSM::Schema;
use File::Spec;
use strict;
use warnings;
$SIG{__DIE__} = sub { print $_[0] unless $_[0] =~ /Fcntl|Unable to satisfy/ };
my $TESTDIR = -d 't' ? 't' : '.';
my $dbf = File::Spec->catfile($TESTDIR,'test.db');
ok my $schema = PSSM::Schema->connect("dbi:SQLite:$dbf"), 'connect to test.db';
$schema->deploy( {add_drop_table => 1} );
ok my ($obj) = $schema->populate(
  'Matrix',
  [ [ qw/name dta/ ],
    ['try', { A => [1,2,3,4], 
	      G => [5,6,7,8],
	      C => [9,10,11,12],
	      T => [13, 14, 15, 16] } ] ]
), 'add a matrix row';
is $obj->name, 'try', 'obj->name correct';
is_deeply $obj->dta->{C}, [9,10,11,12], 'C row of dta is correct (inflation worked)';
ok $obj->residues([sort keys %{$obj->dta}]), 'set residues';
is_deeply $obj->residues, [qw/A C G T/], 'residues set correctly';
$obj->update;
ok my $rs = $schema->resultset('Matrix'), 'get Matrix result set';
ok my $mx = $rs->find( { name => 'try'} ), 'find try object';
is_deeply $mx->residues, [qw/A C G T/], 'residues retrieved correctly';
is_deeply $mx->dta->{G}, [5,6,7,8], 'G row of dta retrieved correctly';


ok my $new_obj = $rs->new_result({}), 'create a new obj';
ok !$new_obj->id, 'has no id';
ok $new_obj->can('dta'), 'new obj can dta()';
ok $new_obj->name('sqelb'), 'name new obj';
ok $new_obj->insert, 'put in db';
ok $new_obj->id, 'now has id';
ok $rs->find({name =>'sqelb'}), 'found new obj';

ok my $prs = $schema->resultset('Predictor'), 'get Predictor result set';

ok my $pr = $prs->create( {name => 'flint'} ), 'new predictor';

ok $mx = $rs->create( { name => 'fred',
			dta => { A => [1,2,3], G => [4,5,6] },
			residues => [qw/A G/],
			length => 2 } ), 'new matrix 1';
ok my $nx = $rs->create( { name => 'wilma',
			dta => { A => [-1,-2,-3], G => [-4,-5,-6] },
			residues => [qw/A G/],
			length => 2 } ), 'new matrix 2';
ok $pr->add_to_matrices($mx, {coeff => 1}), 'add matrix 1 to predictor';
ok $pr->add_to_matrices($nx, {coeff => -1}), 'add matrix 2 to predictor';

ok my ($sx) = $pr->matrices({name => 'wilma'}), 'get matrix wilma from predictor';
is_deeply $sx->dta->{G}, [-4,-5,-6], 'got correct matrix';
ok my $item = $pr->matrix_in_predictor({matrix_id => $sx->id}), 'got linking row (for coeff)';

is $item->next->coeff, -1, 'correct coefficient';
done_testing;
1;

