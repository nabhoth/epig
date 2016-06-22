#!/usr/bin/perl
# requires the path as argument
#$dir_to_open=$ARGV[0];

my $c = 0;
@avg_fitness,@best_fitness, @avg_error, @avg_cost;
@dir_contents;
$dir_to_open, $file_prefix,$rough_mode;
$HEADER = 0;
$GATES = 1;
$RESULTS = 2;


foreach my $arg (@ARGV)
{
	if ($c == 2){
		$rough_mode=$ARGV[$c];
	} elsif ($c == 1){
		$file_prefix=$ARGV[$c];
	} elsif ($c == 0){
		$dir_to_open=$ARGV[$c];
	}
	$c++;
}

if ($c != 3) {print "usage: ./runall.pl <output_path_prefix>  <input_file_prefix> <rough_mode:yes|no>\n";exit;}
else {print "Input File Prefix: $file_prefix, Output Folder Prefix: $dir_to_open\n";}
`echo "" > runall.log`;

	if($file_prefix =~ m/(.*(_[^_]+_).*)*/){
		print $file_prefix.": ".$1."\n";
		$1 =~ m/(.*(_[^_]+_).*)*/;
		print ":: $1 ::\n";
	}
	#exit;

# open the template file
open(FILE, "template.txt") or die("Unable to open file");
# read template file into an array
@template_array = <FILE>;
# close template file 
close(FILE);

# open file with handle DIR
opendir(DIR,$dir_to_open) || die("Cannot open directory !\n");
# Get contents of directory
@dir_contents= readdir(DIR);
# Close the directory
closedir(DIR);

$mode = $HEADER;
$walking = 0;
$maxgen = 0;
$count = 2;
$filecount = 0;
$clock = 0;
@results;
@best_ones;


# Now loop through array and print file names
foreach $file (@dir_contents)
{

	if(!(($file eq ".") || ($file eq "..")))
	{
		if ($file =~ m/$file_prefix/){
			$size = @results;
			$rescounter = 0;
			$gatecounter = 0;
			if($file =~ m/.out$/){ 
				print $dir_to_open.$file."\n";
				open FILE, $dir_to_open."$file" or
					die "Cant open file";
				@lines = <FILE>;
				close FILE;
				$i =0, $j = 0;
				$best_counter = 0;
				foreach $line(@lines)
				{
					if ($mode == $HEADER){
#					print $i."  ".$patt."\n";
						if ($rough_mode =~ m/no/i){
							if (length($template_array[$i]) > 2){
								$patt =$template_array[$i];
								$patt = trim($patt);
								$line = trim($line);
								$patt =~ s/\(.*$//;
								if ($line =~ m/$patt/i){
#						print $j."  ".$line."\n";
									$numpatt = '((\d+)|(\d+\.\d+))$';
									$line =~ m/$numpatt/;
									$nures = $1;
#						print " Found: $nures \n";
									if ($line =~ m/Initialization done/)
									{
										$mode = $GATES;
									} elsif ($line =~ m/Total number of GA generations/){
										$maxgen = $nures;
											print "Maxgen: $maxgen\n";
									}
									$i++;
								}
							} else {$i++;}
						} else {
							if ($line =~ m/Initialization done/){
								$mode = $GATES;
							} elsif ($line =~ m/Total number of GA generations/){
								$numpatt = '((\d+)|(\d+\.\d+))$';
								$line =~ m/$numpatt/;
								$nures = $1;
								$maxgen = $nures;
								print "Maxgen: $maxgen\n";

							}
						}
					} elsif ($mode == $GATES){
						$line = trim($line);
					#$patt1 = '(generation max reached|Solution found)';
						$patt2 = 'Generation at';
						$patt3 = '<- Fitness';
						$patt4 = '<- Error';
						$patt1 = '<- Cost';
						if ($line =~ m/.*-----Results------.*/)
						{
#						print " Gates: $line \n";
							$mode = $RESULTS;
						} elsif ($line =~ m/$patt1/i){
#						print "Found: ". $line."\n";
							if($line =~ m/([\d\.]*)/){
								$best_ones[$gatecounter][2] = $1;
#						print "just walking ".$best_ones[$gatecounter][2]."  \n";
							}
							$walking =1;
						} elsif ($line =~ m/$patt3/i){
#						print "Found: ". $line."\n";
							if($line =~ m/([\d\.]*)/){
								$best_ones[$gatecounter][0] = $1;
#						print "just walking ".$best_ones[$gatecounter][0]."  \n";
							}
							$walking =1;
						} elsif ($line =~ m/$patt4/i){
#						print "Found: ". $line."\n";
							if($line =~ m/([\d\.]*)/){
								$best_ones[$gatecounter][1] = $1;
#						print "just walking ".$best_ones[$gatecounter][1]."  \n";
							}
							$walking =1;
						} elsif($line =~ m/$patt2/i){
							if($line =~ m/(\d+)/){
#						print "Found: '". $1."', ".$line."\n";
								$gatecounter = $1;
							}
						}
#						print " GateCounter: $gatecounter \n";
						if ($walking == 0){
						


						}
						$i++;
					} elsif ($mode == $RESULTS){
#					print " Results: $line \n";
						$line = trim($line);
						$numpatt = '(([\d\.]+)(\s+))(([\d\.]+)(\s+))(([\d\.]+)(\s+))(([\d\.]+))$';
						$charpatt = '.+,.+';
					#print "Size: $size\n";
						if ($line =~ m/$numpatt/){
#					print " Results: 1:|$1| 5:|$5| 7:|$7| $10\n";
							$one = trim($1);
							$five = trim($5);
							$seven =trim($7);
							$ten = trim($10);
#							print "Size: $size  Results: 1:|$one| 5:|$five| 7:|$seven| $ten\n";
							if ($size > 0){
								if ($one == 0 || $five == 0 || $seven == 0){
									$results[0][$rescounter]= $results[0][$rescounter]+1;
									$results[1][$rescounter]= $results[1][$rescounter]+1 ;
									$results[2][$rescounter]= $results[2][$rescounter]+1 ;
									$results[3][$rescounter]= $results[3][$rescounter]+1 ;
								} else {
									$results[0][$rescounter] += $1;
									$results[1][$rescounter] += $5;
									$results[2][$rescounter] += $7;
									$results[3][$rescounter] += $10;
								}
							} else {
								if ($one == 0 && $five == 0){
									push(@results, (1, 1, 1, 1));
								} else { 
									push(@results, ($one, $five, $seven, $ten));
								}
							}
							$rescounter++;
						} elsif ($line =~ m/.*\-\-\-Results\-\-\-.*/)
						{
							$mode = $HEADER;
						}
					}
					$j++;
				}
			$filecount++;
			}
			while ($rescounter < $maxgen){
				$results[0][$rescounter]= $results[0][$rescounter]+1;
				$results[1][$rescounter]= $results[1][$rescounter]+1 ;
				$results[2][$rescounter]= $results[2][$rescounter]+1 ;
				$results[3][$rescounter]= $results[3][$rescounter]+1 ;
				$rescounter++;
			}
		}
	}
}
for ($r = 0; $r < $maxgen; $r++){
	print $best_ones[$r][0].", ";
}
print "Results: \n";
print "Number fo files: $filecount\n";
open FILE, ">".$dir_to_open."stats_out_$file_prefix.txt" or
	die "Cant open file";
$bt_f_one = 0;
$bt_e_one = 1;
$bt_c_one = 0;
for ($r = 0; $r < $maxgen; $r++){
	if($best_ones[$r][0] == 0){
		$best_ones[$r][0] = $bt_f_one;
		$best_ones[$r][1] = $bt_e_one;
		$best_ones[$r][2] = $bt_c_one;
	} else {
		$bt_f_one = $best_ones[$r][0];
		$bt_e_one = $best_ones[$r][1];
		$bt_c_one = $best_ones[$r][2];
	}
}
print FILE "'Gen.', 'Avg. Fit.', 'Best Avg. Fit.',  'Avg. Error',  'Avg. Cost', 'Best Fit.', 'Best Err.',  'Best Cost' "." \n";
for ($r = 0; $r < $maxgen; $r++){
		print FILE $r." ".($results[0][$r]/$filecount)." ".($results[1][$r]/$filecount)." ".(1-($results[2][$r]/$filecount))." ".($results[3][$r]/$filecount)." ".$best_ones[$r][0]." ".$best_ones[$r][1]." ".$best_ones[$r][2]." \n";
	}


close FILE;
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
