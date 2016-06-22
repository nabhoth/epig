#!/usr/bin/perl

$c = 0;
$indata;
@args;
foreach$arg (@ARGV)
{
	if ($c == 0){
		$indata=$ARGV[$c];
	} elsif ($c == 1){
		$runs=$ARGV[$c];
	} elsif ($c == 2){
		$output_path=$ARGV[$c];
			
	}
	$c++;
}

if ($c != 3) {print "$#ARGC usage: ./runall.pl input.file generations> <output_path>\n";exit;}
else {print "Runs: $runs, Output Folder: $output_path\n";}
`echo "" > runall.log`;
for($r=0; $r<$runs;$r++){
	print "Startin $r -(st/th)\n";
	`./epig $indata $output_path >> $output_path.runall.log`;
	print "Ending $r -(st/th)\n";
}	
