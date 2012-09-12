use strict;
use FindBin;
use lib "$FindBin::Bin/build";
use perlforcealigner;

my $model = $ARGV[0];
my $srcL  = $ARGV[1];
my $trgL  = $ARGV[2];

my $mode = $perlforceAligner::SymForceAligner::GrowDiagFinal;
if($ARGV[3]) {
    $mode = eval( "\$perlforcealigner::SymForceAligner::" . $ARGV[3] );
}

my $fa = new perlforcealigner::SymForceAligner($srcL, $trgL, $model);
$fa->setMode($mode);

my ($src, $trg);
while(<STDIN>) {

    if($. % 2) {
	$src = $_;
    }
    else {
	$trg = $_;
	
	print $fa->alignSentenceStr($src, $trg), "\n";
    }

}