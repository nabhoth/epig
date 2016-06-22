#!/usr/bin/perl
# requires the path as argument
use strict;

# Parameters fitype, population, generation are substitutions
# Parameters mincost, segs:minsegs, measus are using a formula to be switched on/off or scaled 
# formula is new_val = (oldVal + param*oldVal). To leave in a given state put 0, to zero put -param and so on
# any double parameter such as segs:minsegs is calculated as max param - (param*input_coef)
my %counter = (#"mincost",[0,0,0.5,1],
		"segs:minsegs",[0, "0:0", "0.1:0.2", "0.2:0.4", "0.3:0.5"],
		"erra",[0,0.5,0.6,0.7,0.8,0.9,1]);#,
		#"measu",[0,0],
		#"fitype",[0,1,2,3],
		#"pop",[0, 100],
		#"gener",[0, 500,900],
		#"alph:beta",[0, "0.01:-0.01"]);

my $error = 0.01;
my @results;
my $c = 0;
my $file_to_open;
my $runs;
my (@avg_fitness,@best_fitness, @avg_error, @avg_cost);
my @dir_contents;
my $HEADER = 0;
my $STATES = 1;
my $DONE = -1;
my @target_content;
my $out_dir;
my $out_prefix;

foreach my $arg (@ARGV)
{
	if ($c == 1){
		$file_to_open=$ARGV[$c];
	} elsif ($c == 0){
		$runs=$ARGV[$c];
	} elsif ($c == 2){
		$out_prefix=$ARGV[$c];
	}
	$c++;
}

if ($c != 3) {
	print "usage: ./runall.pl <number of runs>  <input.file> <output_path_prefix>\n";
	exit;
} else {
	print "Runs: $runs, Output Folder Prefix: $out_prefix\n";
}

my $out_file = trim($file_to_open);

#my @nicks = ("mincost","segs:minsegs","measu","fitype","pop","gener","alph:beta");
#my @nicks = ("erra");
my @nicks = ("segs:minsegs","erra");
#sets the total number of iterations
my $cc = 0;
my @max_count;
my @tem_arr = %counter;
#for (my $h = 0; $h < @tem_arr; $h++){
for (my $h = 0; $h < @nicks; $h++){
	print "'".$nicks[$h]."', ";
	my $nck = $counter{$nicks[$h]};
	my $tu = 0;
	#print $h.": '".$nck."', ";
	for (my $p = 0; $p < 15; $p++){
		#$tu+= 1 if defined $tem_arr[$h+1]->[$p];
		$tu+= 1 if defined $nck->[$p];
	}
	$max_count[$cc++] = $tu-1;
	#$h = $h + 1;
}
#my @max_count = (length($counter{$nicks[0]})-1,length($counter{$nicks[1]})-1,length($counter{$nicks[2]})-1,length($counter{$nicks[3]})-1,length($counter{$nicks[4]})-1,length($counter{$nicks[5]})-1);
print "\n\n";
print @max_count."\n";
print " ".$max_count[0]." ".$max_count[1]." ".$max_count[2]." ".$max_count[3]." ".$max_count[4]." ".$max_count[5]."\n";

my %hashes_in_file;
# open the template file
open(FILE, "assigntemplate.txt") or die("Unable to open file");
# read template file into an array

my @raw_array = <FILE>;
# close template file 
close(FILE);
my %template_array;
my $seppat = '(\w+)(:)(.+)$';
my $i = 0;
foreach my $LINE (@raw_array){
	my $line = trim($LINE);
	if ($line =~ m/$seppat/){
		$template_array{$1} = "$3";
		$hashes_in_file{$1} = $i;
	}	
	$i++;

}

$i = 1;
for (my $g = 0; $g < @max_count; $g++){
	$i *= $max_count[$g];
}

print "Running $i cofnigurations\n";
for(my $j = 0; $j < $i; $j++){

	$out_dir = $out_prefix;
	$out_file = loadTarget();
	nextInput();
	for(my $r=0; $r<$runs;$r++){
		print "Startin $r -(st/th)\n";
		print "File: ".$out_file.", Out Dir: ".$out_dir."/\n";
		`./epig $out_file $out_dir/ >> $out_dir/.runall.log`;
		print "Ending $r -(st/th)\n";
	}	
	`rm $out_file`; 
}


sub loadTarget(){
use vars qw($file_to_open);
# open file with handle DIR
open(TARGET,$file_to_open) || die("Cannot open file !\n");
# Get contents of directory
@target_content= <TARGET>;
# Close the directory
close(TARGET);

my $out_file;
if ($file_to_open =~ /(.+)\.(.+)/m){
	$out_file = trim($1);
} else {
	$out_file = trim($file_to_open);
}
return $out_file;
}

#generate input file for the next counter
sub nextInput (){

	use vars qw($out_dir $out_file %counter @nicks @target_content);
	my $ext = 0;
	my $carry = 1;
	my $count = 0;
	print "\n";

	foreach my $hashname(@nicks){
		print "Hash: ".$hashname." <> ";
		my $param_row = $counter{"$hashname"};
		my $index = $param_row->[0];
		my $param = $param_row->[$index+1];

		print $param." <> ";#.$param0." <> ".$param1."\n";
		my ($minsegval, $maxsegval, $segval);
		if ($param =~/(.+):(.+)/m){
			#print "Found: ".$param."\n";
			$minsegval = trim($1);	
			$maxsegval = trim($2);	
		} else {
			$segval = int($param);
		}
		my $patt_main = "(fitype|pop|gener)";
		my $patt_sec = "(measu|mincost|segs:minsegs)";
		my $dir_patt = "(measu|pop|gener)";
		my $double_num_patt = "(alph:beta)";
		my $entang_patt = "(erra)";
		if ($hashname =~/(\w+):(\w+)/m){
			my $maxseg = trim($1);
			my $minseg = trim($2);
			my $maxindex = $hashes_in_file{"$maxseg"};
			my $minindex = $hashes_in_file{"$minseg"};
			$maxseg = trim($target_content[$maxindex]);
			$minseg = trim($target_content[$minindex]);
			print $minindex." <> ".$maxindex." :: ".$minsegval." <> ".$maxsegval." ::old  ".$minseg." <> ".$maxseg." ::new ";
			if ($hashname =~/$double_num_patt/m){
				$maxseg = $maxseg + $maxsegval;
				$minseg = $minseg + $minsegval;
			} else {
				$maxseg = $maxseg - round(eval(trim($target_content[$maxindex]))*$maxsegval);
				$minseg = $minseg - round(eval(trim($target_content[$minindex]))*$minsegval);
			}
			print $minseg." <> ".$maxseg."\n";
			$target_content[$maxindex] =  $maxseg."\n"; 
			$target_content[$minindex] =  $minseg."\n"; 
			$out_file .="_".$hashname."_";
			$out_file .=$minseg.":".$maxseg."_";
		} else{
			if ($hashname =~ /$patt_main/m){
				my $lineindex = $hashes_in_file{"$hashname"};
				print $lineindex." <>Old: ".trim($target_content[$lineindex])." ::new  ".$segval."\n";
				$target_content[$lineindex] =  $segval."\n"; 
				if ($hashname =~ /$dir_patt/m){
					$out_dir .="_".$hashname."_";
					$out_dir .=$segval."_";
				} else {
					$out_file .="_".$hashname."_";
					$out_file .=$segval."_";
				}
			} else {
				my $lineindex = $hashes_in_file{"$hashname"};
				my $val = trim($target_content[$lineindex]);
				print "Val: ".$val.",  ".$segval.",  ".$param."\n";
				print $hashname."\n";
				my $segvalloc = 0.0;
				if ($hashname =~/$entang_patt/m){
					$segvalloc = $val - $param;
				} else {
					$segvalloc = abs($val + ($segval*$val));
				}
				print $lineindex." <>Old ".trim($target_content[$lineindex])." ::new  ".$segvalloc."\n";
				$target_content[$lineindex] =  $segvalloc."\n"; 
				if ($hashname =~ /$dir_patt/m){
					$out_dir .="_".$hashname."_";
					$out_dir .=$segvalloc."_";
				} else {
					$out_file .="_".$hashname."_";
					$out_file .=$segvalloc."_";
				}
		}
			
		}
		if ($carry > 0){
			if ($max_count[$count] > 0){
				$index++;
				if ($index > $max_count[$count]-1){
					$param_row->[0] = 0;
					$carry = 1;
				} else {
					$param_row->[0] = $index;
					$carry = 0;
				}
			}
			$count++;
		}
	}

	$out_file = substr($out_file, 0, -1).".inn";
	$out_dir = substr($out_dir, 0, -1);
	print "Dir: ".$out_dir.", File: ".$out_file."\n";
	`mkdir $out_dir`;
	my $output = $out_dir."/".$out_file;
	open (TARGET,">$out_file") || die("Cannot open file !\n");
	foreach my $line(@target_content){
		print TARGET $line;
	}
	close(TARGET);

}
# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
# Left trim function to remove leading whitespace
sub ltrim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($)
{
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}

# Split String with a given character
sub rtrim($)
{
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}
sub round {
	my($number) = shift;
	return int($number + .5);
}
