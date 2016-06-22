#!/usr/bin/perl
# requires the path as argument
#$dir_to_open=$ARGV[0];

my $c = 0;
@result_contents;
$conversion, $file_source;
$HEADER = 0;
$DATA = 1;


foreach my $arg (@ARGV)
{
	if ($c == 1){
		$file_source=$ARGV[$c];
	} elsif ($c == 0){
		$conversion=$ARGV[$c];
	}
	$c++;
}

if ($c != 2) {print "usage: ./convrevlibga.pl <measure=1/matrix=0> <input_file> \n";exit;}

# open the template file
open(FILE, $file_source) or die("Unable to open file");
# read template file into an array
@source_array = <FILE>;
# close template file 
close(FILE);

$mode = $HEADER;
$variables = 0;
@garbage;
$count = 2;

$patt_var = '\.numvars (\d+)';
$patt_garbage = '\.garbage (.+)';
$begin_patt = '\.begin';
$end_patt = '\.end';

foreach $source_line (@source_array){
	if ($mode == $HEADER){
		print ($source_line."\n");
		if ($source_line =~m/$patt_var/){
			$variables = $1;
		}
		if ($source_line =~m/$begin_patt/){
			$mode = $DATA;
		}
		if ($source_line =~m/$patt_garbage/){
			$garb = $1;
			$garb_length = length $garb;
			for ($r = 0; $r < $garb_length-1; $r++){
				$char = substr($garb, $r, 1);
				push(@garbage, $char);
				print @garbage[$r].", ";
			}
		}
	} elsif ($mode == $DATA) {
		if ($source_line =~m/$end_var/){
			break;
		}
		if ($conversion == 1){#build a measurement based specification
			$length = @garbage;
			if (substr($source_line, $length, 1) != '-') 
				print substr($source_line, $length, 1)."\n";

		} else {#build a matrix based specification




		}	




	}



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
