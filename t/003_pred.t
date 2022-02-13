#$Id: 003_pred.t 344 2014-02-09 23:03:08Z maj $
use Test::More;
use Test::Exception;
use File::Spec;
use lib '../lib';
use PSSM;
use PSSM::Schema;
use PSSM::Predictor;
use PSSM::Matrix;
use strict;
use warnings;


my $TESTDIR = -d 't' ? 't' : '.';
my $dbf = File::Spec->catfile($TESTDIR,'test.db');
ok my $schema = PSSM::Schema->connect("dbi:SQLite:$dbf"), 'connect to test.db';
$schema->deploy( {add_drop_table => 1} );
ok my $m1 = PSSM::Matrix->build_new(
  $schema, 
  { A => [ 0.2, 0.3, 0.5, 0.1 ],
    G => [ 1.1, 0.5, 0.7, -0.3],
    C => [ -0.1, -0.7, 1.5, 2.0 ],
    T => [ 0.5, 0.5, 0.5, -1.0 ] } ), 'test mx 1';
ok my $m2 = PSSM::Matrix->build_new($schema, 0.5*$m1), 'test mx 2';
ok my $m3 = PSSM::Matrix->build_new(
  $schema, 
  { A => [ 0.2, 0.3, 0.5, 0.1 ],
    G => [ 1.1, 0.5, 0.7, -0.3] }
 ), 'test mx 3';
ok $m1->store('m1'), 'store mx 1';
ok $m3->store('m3'), 'store mx 3';
my $p1;
throws_ok { PSSM::Predictor->build_new($schema, MXS => [$m1, $m2],
		       COEFFS => [1, -2]) } qr/must be stored/;
ok $m2->store('m2'), 'store mx 2';
ok $p1 = PSSM::Predictor->build_new(
  $schema, 
  MXS => [$m1, $m2],
  COEFFS => [1, -2]),'create new predictor';

is_deeply $p1->residues, [qw/A C G T/], 'predictor residues correct';
is_deeply $p1->dta->{A}, [0,0,0,0];
is_deeply $p1->dta->{G}, [0,0,0,0];
is_deeply $p1->dta->{C}, [0,0,0,0];
is_deeply $p1->dta->{T}, [0,0,0,0], 'all predictor columns are correct (linear combination correct)';
throws_ok { PSSM::Predictor->build_new(
  $schema, 
  MXS => [$m1, $m3], COEFFS => [1,1]
 ) } qr/must conform/;

is $p1->simple_score("AGGC"), 0, "simple score (from PSSM) works";

ok $p1->label('lbl1'), 'set predictor label field';
is $p1->label, 'lbl1', 'label correctly set';
ok $p1->scr95([-10, 21.05]), 'set scr95 (95% score pctile) field';
is_deeply $p1->scr95, [-10,21.05], 'scr95 correctly set';
ok $p1->cutoff([-2.0,-2.0]), 'set cutoff interval';
is_deeply $p1->cutoff, [-2.0,-2.0], 'cutoff interval correctly set';
ok $p1->plbl([qw/R5 X4/]), 'set prediction label array field';
is_deeply $p1->plbl, [qw/R5 X4/], 'prediction label array correctly set';
ok $p1->store('p1'), 'store predictor in db';
ok my $p2 = PSSM->get_predictor_by_name($schema,'p1'), 'retrieved it';
is_deeply $p2->residues, [qw/A C G T/], 'predictor residues correct';
is_deeply $p2->dta->{A}, [0,0,0,0];
is_deeply $p2->dta->{G}, [0,0,0,0];
is_deeply $p2->dta->{C}, [0,0,0,0];
is_deeply $p2->dta->{T}, [0,0,0,0], 'all predictor columns are correct (linear combination correct)';
is $p2->label, 'lbl1', 'retrieved label';
is_deeply $p2->scr95, [-10,21.05], 'retrieved scr95';
is_deeply $p1->cutoff, [-2.0,-2.0], 'retrieved cutoff interval';

done_testing;
