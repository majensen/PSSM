#$Id: 004_nwx.t 205 2013-07-29 04:07:41Z maj $
use Test::More;
use Test::Exception;
use File::Spec;
use lib '../lib';
use PSSM::Schema;
use PSSM::Matrix;
use PSSM::NW;
use Bio::SeqIO;
use strict;
use warnings;
#$SIG{__DIE__} = sub { print $_[0] };
my $TESTDIR = -d 't' ? 't' : '.';
my $dbf = File::Spec->catfile($TESTDIR,'test.db');
my $seqf = File::Spec->catfile($TESTDIR,'samples','simple.fasta');
ok my $schema = PSSM::Schema->connect("dbi:SQLite:$dbf"), 'connect to test.db';
$schema->deploy( {add_drop_table => 1} );
my $mxf = File::Spec->catfile($TESTDIR, qw/samples aureus.mx/);

ok my $mx = PSSM::Matrix->build_new($schema, 
				 { A => [.5, .5, .5],
				   C => [.5, .5, .5],
				   G => [.5, .5, 1],
				   T => [1, 1, .5] }), 'build NW pssm';
$mx->name('simple');
ok my $nwx = PSSM::NW->build_new($schema, MATRIX => $mx), 'create NW matrix using simple mx';
isa_ok $nwx, "PSSM::NW";
is $nwx->gpen, 2, 'gpen default';
is $nwx->gxpen, 0, 'gxpen default';
is $nwx->gpen_fl, 0, 'gpen_fl default';
is $nwx->gpen_in, 8, 'gpen_in default';
ok $nwx->fsfixed, 'fsfixed default';
is $nwx->fsfixed(0),0, 'clear fsfixed for test';
is $nwx->gpen_in(2), 2, 'set gpen_in for test';
is_deeply $nwx->mapping, [1,1], 'mapping default';
ok my $seqio = Bio::SeqIO->new(-file => $seqf), 'seqio';
ok my $seq = $seqio->next_seq;
ok $nwx->align($seq), 'aligned seq';
is $nwx->lb->[0], 2, 'lower coord correct';
is $nwx->ub->[0], 4, 'upper coord correct';
my @info;
while (my $info = $nwx->get_align ) {
  push @info, $info;
}
is scalar @info, 1, 'got 1 alignment';
is $info[0]->{mxc}, '-TTG-', 'alignment correct';
is $info[0]->{lb}, 2, 'lower coord correct';
is $info[0]->{ub}, 4, 'upper coord correct';

# is_deeply $PSSM::NW::OBJTABLE{$nwx->id}{_dist},
#   [ [0, 2, 4, 6],
#     [0,-.5,1.5,3.5],
#     [0,-.5,-1,1],
#     [0,-1,-1.5,-1.5],
#     [0,-.5,-1.5,-2.5],
#     [0,-.5,-1,-2.5] ], 'nw dist table correct';
# is_deeply $PSSM::NW::OBJTABLE{$nwx->id}{_pthm},
#   [ [0,1,1,1],
#     [2,4,5,5],
#     [2,4,4,5],
#     [2,4,4,4],
#     [2,4,4,4],
#     [2,4,4,2] ], 'nw path table correct';

is $nwx->fsfixed(1), 1, 'set fsfixed for test';
ok $nwx->align($seq), 'aligned seq';
# is_deeply $PSSM::NW::OBJTABLE{$nwx->id}{_pthm},
#   [ [0,1,1,1],
#     [2,4,4,4],
#     [2,4,4,4],
#     [2,4,4,4],
#     [2,4,4,4],
#     [2,4,4,2] ], 'nw path table correct w/fs fixed';

$seqf = File::Spec->catfile($TESTDIR,'samples','sequence.fasta');
$seqio = Bio::SeqIO->new(-file => $seqf);
$seq = $seqio->next_seq;
$mx = PSSM::Matrix->build_new($schema, $mxf);
$mx->name('aureus');
ok $nwx = PSSM::NW->build_new($schema, MATRIX => $mx), 'create NW matrix using sample mx file';
ok $nwx->align($seq), 'larger sequence transits algorithm';
@info = ();
while (my $info = $nwx->get_align ) {
  push @info, $info;
  diag $info->{lb}." ".$seq->subseq($info->{lb}, $info->{ub});
}
is scalar @info, 4, 'found 4 exact repeats';
done_testing;
1;
