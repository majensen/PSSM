#$Id: alg.pm 205 2013-07-29 04:07:41Z maj $
package PSSM::NW::alg;
use base 'Exporter';
use feature qw/switch/;
use constant { II => 4, ZI => 1, IZ => 2 };
use strict;
use warnings;

our @EXPORT = qw/nw_align fpaths apply_path/;

my $GAPCHAR = '-';

sub nw_align {
 my ($s, $mx, $parms) = @_;
 my ($mlen, $gpen_in, $gpen_fl, $fsfixed, $sspan, $mspan) = @$parms;
 my ($dist,$pthm,$scr);

 my ($slen,$c,$i,$j,$k);
 my ($h11, $h10, $h01,$min);
 $mspan ||=1;
 $sspan ||=1;
 $mlen /= $mspan;
 $slen = length($s)/$sspan;
 $$dist[0][0] = $$pthm[0][0] = 0;
 for ($c=0,$i=1;$i<=$mlen;$i++) {
   $c += $gpen_in*$mspan;
   $$dist[0][$i]= $c;
   $$pthm[0][$i] = ZI;
 }
 for ($c=0,$i=1;$i<=$slen;$i++) {
   $c += $gpen_fl*$sspan;
   $$dist[$i][0] = $c;
   $$pthm[$i][0] = IZ;
 }

 $$scr[0][0]=0;

 for ($i=1;$i<=$slen;$i++) {
   for ($j=1; $j<=$mlen; $j++) {
     for ($c=0,$k=0;$k<$mspan;$k++) {
       $c += $mx->{substr($s,$sspan*($i-1)+$k,1)}[$mspan*($j-1)+$k];
     }
     $h01 = $$dist[$i][$j-1] + $gpen_in*$mspan;
     $h10 = $$dist[$i-1][$j] + (($j==$mlen) ? $gpen_fl : $gpen_in)*$mspan;
     $h11 = $$dist[$i-1][$j-1] - $c;
     $min = ($h01 > $h10) ? $h10 : $h01;
     $min = ($min > $h11) ? $h11 : $min;
     $$pthm[$i][$j] = 0;
     $$dist[$i][$j] = $min;
     if ($min == $h01) {
       $$scr[$i][$j] = $$scr[$i][$j-1];
       unless ($fsfixed) {
	 $$pthm[$i][$j] |= ZI;
       }
     }
     if ($min == $h10) {
       $$scr[$i][$j] = $$scr[$i-1][$j];
       unless ($fsfixed && ($i<=$j)) {
	 $$pthm[$i][$j] |= IZ;
       }
     }
     if ($min == $h11) {
       $$scr[$i][$j] = ($$scr[$i-1][$j-1]||0)+$c;
       $$pthm[$i][$j] |= II;
     }
   }
 }
 return $$dist[$slen][$mlen], $$scr[$slen][$mlen], fpaths($slen, $mlen, $pthm);
 1;
}

sub apply_path {
  my ($pth, $s, $m, $sspan, $mspan) = @_;
  my ($i, $j) = (0,0);
  my ($lb,$ub) = (-1,-1);
  my ($gs, $gm) = ('','');
  $sspan ||=1;
  $mspan ||=1;
  while (my $c = pop @$pth) {
    no warnings qw(experimental);
    given ($c) {
      when (1) {
	$gs .= $GAPCHAR x $sspan;
	$gm .= substr($m, ($mspan*$j++), $mspan);
      };
      when (2) {
	$gs .= substr($s, ($sspan*$i++), $sspan);
	$gm .= $GAPCHAR x $mspan;
      };
      when (3) {
	if ($lb < 0) { #not set
	  $lb = ($i*$sspan)+1;
	}
	$ub = ($i+1)*$sspan;
	$gs .= substr($s, ($sspan*$i++), $sspan);
	$gm .= substr($m, ($mspan*$j++), $mspan);
      };
      default {
	die 'apply_path - unknown path element';
      };
    }
  }
  return ($gs,$gm,$lb,$ub);
}

sub fpaths {
  my ($slen,$mlen,$pthm) = @_;
  my @stacks = ([$mlen,$slen]);
  my ($i, $j, @paths);
  while ( my $stk = pop @stacks ) {
    $i = pop @$stk;
    $j = pop @$stk;
    my $pth = $$pthm[$i][$j];
    unless ($pth) {
      push @paths, $stk;
      next;
    }
    my $ct = (($pth & 4) && 1) + (($pth & 2) && 1) + ($pth & 1);
    my @stks = ($stk);
    while (--$ct) {
      # copy stk
      push @stks, [@$stk];
    }
    foreach my $stkc (@stks) {
      if ($pth & IZ) {
	push @$stkc, 2, $j, $i-1;
	$pth ^= IZ; # clear flag
	next;
      }
      if ($pth & ZI) {
	push @$stkc, 1, $j-1, $i;
	$pth ^= ZI;
	next;
      }
      if ($pth & II) {
	push @$stkc, 3, $j-1, $i-1;
	$pth ^= II;
	next;
      }
    }
    push @stacks, @stks;
  }
  return @paths;
}

=head1 NAME

PSSM::NW::alg - Needleman-Wunsch alignment to PSSM

=head1 SYNOPSIS

See L<METHODS>.

=head1 DESCRIPTION

=head1 METHODS

=over

=item nw_align()

 $parm_array = [ $matrix_len, $gpen_in, $gpen_fl, $fs_fixed, $sspan, $mspan ];
 ($aln_score, $pssm_score, @paths) = nw_align($sequence_string,
   $matrix_hash, $parm_array);

=item fpaths()

 @paths = fpaths($sequence_len, $matrix_len, $path_matrix);

=item apply_path()

 ($seq_with_gaps, $matrix_seq_with_gaps, $lower_seq_coord, $upper_seq_coord) = 
  apply_path( $path, $sequence_string, $matrix_sequence_string,
              $seq_span, $matrix_span);

=back

=cut
