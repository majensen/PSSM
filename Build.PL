#$Id: Build.PL 335 2014-01-27 05:23:35Z maj $
use Module::Build;
use strict;
use warnings;

my $build = Module::Build->new(
  dist_name => 'PSSM',
  dist_version_from => 'lib/PSSM.pm',
  dist_abstract => 'DBIx::Class for Position-specific scoring matrices',
  perl => '5.010',
  license => 'perl',
  config_requires => {
    'Module::Build' => 0
   },
  build_requires => {
    'Test::More' => 0,
    'Test::Exception' => 0,
    'File::Spec' => 0,
    'Bio::SeqIO' => 0
   },
  requires => {
    'DBIx::Class' => 0,
    'SQL::Translator' => 0,
    'List::MoreUtils'=> 0,
    'JSON' => 2
  },
);

$build->create_build_script;
1;
