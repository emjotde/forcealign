use strict;
use FindBin;
use lib "$FindBin::Bin/build";
use perlforcealigner;

my $model = $ARGV[0];
my $srcL  = $ARGV[1];
my $trgL  = $ARGV[2];

my $srcFile = $ARGV[3];
my $trgFile = $ARGV[4];

open(SRC, $srcFile) or die;
open(TRG, $trgFile) or die;

my $mode = $perlforceAligner::SymForceAligner::GrowDiagFinal;
if($ARGV[3]) {
    $mode = eval( "\$perlforcealigner::SymForceAligner::" . $ARGV[3] );
}

my $fa = new perlforcealigner::SymForceAligner($srcL, $trgL, $model);
$fa->setMode($mode);

my $c = 0;
while(defined(my $src = <SRC>) and defined(my $trg = <TRG>)) {
    chomp($src, $trg);
    print $fa->addSentence($src, $trg);
    last if $c++ > 10000;
}
$fa->alignCorpus();