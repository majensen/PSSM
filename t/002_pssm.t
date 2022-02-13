#$Id: 002_pssm.t 193 2013-07-14 03:03:28Z maj $
use Test::More;
use Test::Exception;
use File::Spec;
use lib '../lib';
use PSSM;
use PSSM::Schema;
use strict;
use warnings;

#$SIG{__DIE__} = sub { print $_[0] };

use_ok('PSSM::Matrix');
my $TESTDIR = -d 't' ? 't' : '.';
my $dbf = File::Spec->catfile($TESTDIR,'test.db');
ok my $schema = PSSM::Schema->connect("dbi:SQLite:$dbf"), 'connect to test.db';
$schema->deploy( {add_drop_table => 1} );
ok $schema->populate(
  'Matrix',
  [ [ qw/name residues length dta/ ],
    ['try',
     [qw/A C G T/],
     4,
     { A => [1,2,3,4], 
       G => [5,6,7,8],
       C => [9,10,11,12],
       T => [13, 14, 15, 16] } ] ]
), 'add a matrix row';

ok my $mx = PSSM->get_matrix_by_name($schema,'try'), 'got matrix by name';
isa_ok($mx, 'PSSM::Matrix');
is_deeply $mx->residues, [qw/A C G T/], 'residues correct';
ok my $nmx = PSSM::Matrix::Cmult(-1,$mx), 'cmult executed';
is_deeply $nmx->dta->{A},[-1,-2,-3,-4], 'operation correct';
ok $nmx = -1*$nmx, 'use overloaded *';
is $nmx, $mx, 'use overloaded ==, correct';
ok my $amx = $nmx + 5*$mx, 'linear combination executed';
is_deeply $amx->dta->{G}, [30, 36, 42, 48], 'operation correct';
ok $amx = 5*$mx - $nmx, 'use overloaded -';
is_deeply $amx->dta->{T}, [52, 56, 60, 64], 'operation correct';

ok my $bmx = $schema->source('Matrix')->resultset->new_result({}), 'new Matrix obj';
ok $bmx->build( { 'A' => [1,2],
		  'C' => [3,4],
		  'G' => [5,6] } ), 'PSSM::Matrix::build()';
is_deeply $bmx->residues, [qw/A C G/], 'build worked : residues()';
is_deeply $bmx->dta, { 'A' => [1,2],
		  'C' => [3,4],
		  'G' => [5,6] }, 'build worked : dta()';
is $bmx->length, 2, 'build worked : length()';
throws_ok { $bmx + $amx } qr/do not conform/, 'sum of nonconforming mxs barfs ok';
ok $bmx->store('little'), 'put matrix in db with store';
ok my $cmx = PSSM->get_matrix_by_name($schema,'little'), 'found it';
is $cmx, $bmx, 'and it\'s the same as what we put';

is $cmx->simple_score("AC"), 5, "score ok";
is_deeply [$cmx->simple_score('GG', 'GA', 'CG')], [11,7,9], "score array of seqs ok";
is $cmx->consensus, 'GG', 'consensus ok';

done_testing;
