#$Id: 0040_alg.t 205 2013-07-29 04:07:41Z maj $
use Test::More;
use lib '../lib';
use PSSM::NW::alg;
use strict;
use warnings;

my $s = 'AATGC';
my $mx = { A => [.5, .5, .5],
	   C => [.5, .5, .5],
	   G => [.5, .5, 1],
	   T => [1, 1, .5] };

ok my ($aln_score, $score, @paths) = nw_align($s, $mx, [3, 2, 0, 0]),'run aligner';
is $aln_score, -2.5, 'alignment score correct';
is $score, 2.5, 'pssm score correct';
is @paths, 1, 'got 1 isoscoring alignment';

ok my ($gs, $gm, $lb, $ub) = apply_path($paths[0], $s, 'TTG'), 'apply path to seqs';
is $gs, 'AATGC';
is $gm, '-TTG-';
is $lb, 2;
is $ub, 4;

my $pthm =   
  [ [0,1,1,1],
    [2,4,5,5],
    [2,4,4,5],
    [2,5,4,4],
    [2,4,4,4],
    [2,4,4,6] ];

@paths = fpaths(5, 3, $pthm);
is @paths, 3, "got 3 paths";
foreach (@paths) {
  my ($gs, $gm) = apply_path($_, $s, 'TTG');
  1
}

$s = 'AAATTTGGGCCC';
ok( ($aln_score, $score, @paths) = nw_align($s, $mx, [3, 2, 0, 0,1,1]),'run aligner : mapping 1:1' );
ok (($gs, $gm, $lb, $ub) = apply_path($paths[0], $s, 'TTG',1,1), 'apply path');
is $score, 3, 'score correct';
is $lb, 5, 'lower boundary correct';
is $ub, 7, 'upper boundary correct';
  

ok (($aln_score, $score, @paths) = nw_align($s, $mx, [3, 2, 0, 0,3,3]),'run aligner : mapping 3:3');
foreach (@paths) {
  ($gs, $gm, $lb, $ub) = apply_path($_, $s, 'TTG',3,3);
  1;
}

done_testing;
1;
