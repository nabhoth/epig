#!/usr/bin/perl
# requires the path as argument
#$dir_to_open=$ARGV[0];

my $c = 0;
$file_prefix;
$only_best;
$HEADER = 0;
$GATES = 1;
$RESULTS = 2;
$gates_start = "Number of Input Gates:";

$gate_rep = 'Representation:\s+(.+)$';
$gate_cost = 'Cost of the gate:\s+(.+)$';
$gate_name = 'The name of this gate is:\s+(.+)$';
$gate_io = 'IO of the gate:\s+(.+)$';

$result_best = "Solution found";
$result_start = "Output of synthesis:";
$result_init = "----- Best Individual -----";
$result_rep = "<- Representation";
$result_fitness = "<- Fitness";
$result_error = "<- Error";
$result_cost = "<- Cost";
$result_io = "<- number of inputs";

$replace_I = '& \qw';
$replace_H = '& \gate{H}';
$replace_T = '& \gate{T}';
$replace_Tt = '& \gate{T^\dagger}';
$replace_NOT = '& \targ';
$replace_V = '& \gate{V}';
$replace_Vt = '& \gate{V^\dagger}';
$replace_U = '& \gate{U}';
$replace_Y = '& \gate{Y}';
$replace_Z = '& \gate{Z}';
$replace_C = '& \ctrl{1}';
$replace__C = '& \ctrl{-1}';
$replace_CI = '& \ctrl{2}';
$replace__IC = '& \ctrl{-2}';
$replace_CII = '& \ctrl{3}';
$replace__IIC = '& \ctrl{-3}';
$replace_CIII = '& \ctrl{4}';
$replace__IIIC = '& \ctrl{-4}';
$replace_CIIII = '& \ctrl{5}';
$replace__IIIIC = '& \ctrl{-5}';
$replace_CIIIII = '& \ctrl{6}';
$replace__IIIIIC = '& \ctrl{-6}';

$patt_control = '(I|C)';
$patt_data = '(T|T\+|H|V|V\+|Z|Y|X|NOT|WIRE)';


foreach my $arg (@ARGV)
{
	if ($c == 0){
		$file_prefix=$ARGV[$c];
	} elsif ($c == 1){
		if ($ARGV[$c] > 0){
			$only_best = -1;
		} else {
			$only_best = 1;
		}
	}
	$c++;
}

if ($c != 2) {print "usage: ./translate_circuits.pl <input_file> <mode>\n";exit;}
else {print "Input File: $file_prefix $only_best\n";}


#`echo "" > runall.log`;

# open the template file
open(FILE, "<$file_prefix") or die("Unable to open file");
# read template file into an array
@file_array = <FILE>;
# close template file 
close(FILE);

$mode = $HEADER;
our %gates = ();
@results;
@result;
$j = 0;
$gate_n,$gate_r,$gate_i,$gate_c;
$res_n,$res_r,$res_i,$res_c,$res_f,$res_e;

# Now loop through array and print file names
foreach $line (@file_array)
{
		if ($mode == $HEADER){
			if ($line =~ /$gates_start/m){
				$mode = $GATES;
				print "Mode: GATES\n";
			}
		} elsif ($mode == $GATES){
			if ($line =~ m/$gate_rep/){
				$gate_r = $1;
			} elsif ($line =~ m/$gate_cost/){
				$gate_c = $1;
			} elsif ($line =~ m/$gate_name/){
				$gate_n  = $1;
				print "$gate_r  $gate_n \n";
				$gates->{$gate_r}->{rep} = $gate_r;
				$gates->{$gate_r}->{name} = $gate_n;
				$gates->{$gate_r}->{cost} = $gate_c;
				$gates->{$gate_r}->{io} = $gate_i;
			} elsif ($line =~ m/$gate_io/){
				$gater_i = $1;
			} elsif ($line =~ m/$result_start/){
				$mode = $RESULTS;
			}
		} elsif ($mode == $RESULTS){
			$line = trim($line);
			if ($line =~ m/.*---Results---.*/)
			{
				@out_results = %gates;
				print " Decoded Circuits: out_results \n";
				print @out_results;
				exit;
			} elsif ($line =~ m/$result_best/i){
				$only_best = -1;
			} elsif ($line =~ m/$result_rep/i && $only_best == -1){
				if($line =~ m/([^\s]+)/){
					print "Found: ". $1.", ".$line."\n";
					$res_r = trim($1);
					decode($res_r, $res_i);
					$best_ones[$gatecounter][0] = $1;
				}
			} elsif ($line =~ m/$result_fitness/i && $only_best == -1){
				if($line =~ m/([\d\.]*)/){
					print "Found: ". $line."\n";
					$res_f = trim($1);
					$j++;
				}
			} elsif($line =~ m/$result_error/i && $only_best == -1){
				if($line =~ m/([\d\.]*)/){
					print "Found: '". $1."', ".$line."\n";
					$res_e = trim($1);
				}
			} elsif($line =~ m/$result_cost/i && $only_best == -1){
				if($line =~ m/(\d+)/){
					print "Found: '". $1."', ".$line."\n";
					$res_c = trim($1);
				}
			} elsif($line =~ m/$result_io/i && $only_best == -1){
				if($line =~ m/(\d+)/){
					print "Found: '". $1."', ".$line."\n";
					$res_i = trim($1);
				}
			}
		} 
}

sub decode
{
	my $rep = shift;
	my $io = shift;
	my $result = '';
	my $octave_result = '';
	my @qcircuit = ();
	$qcircuit[0] = '\[\Qcircuit @C=1em @R=1em {';
	$count = 0;$pos = 0;
	$ncount = 0; $npos = 0;$nfunct = 0;
	print "Decoding \n";
	while ( ($pos = index($rep, "p", $pos)) != -1) {
		$count++;
		$pos++;
	} 
	for ($r = 0; $r < $count/2; $r++){
		if ($rep =~ m/(p([^p]+)p)/g)
		{
			$str = $2;
			$result .= '(';
			if ($r == 0){
	#			$octave_result .= 'kron(';
			} else {
				$octave_result .= '*';
			}
			$i_count = 0;
			$kroneckers = 0;
			for ($t = 0; $t < length $str; $t++){
				if ($t > 0) {$result .= ',';}
				$ch = substr($str, $t,1);
				$result .= $gates->{$ch}->{'name'};
				my $octave_gate = $gates->{$ch}->{'name'};
				$octave_gate =~ s/\+/t/;
				if ($t > 0) {
					if ($t < (length $str)-1 ){
						$octave_result .= ',kron('.$octave_gate;
						$kroneckers++;
					}else{
						$octave_result .= ','.$octave_gate;
					}
				} else {
					$octave_result .= 'kron('.$octave_gate;
				}
				#Encode the Qcircuit
				$gate_name = $gates->{$ch}->{'name'};
				$nfunct = 0;
				$wires = 0;
				$npos = 0;
				while ( $npos < length $gate_name ) {
					$next = substr ($gate_name, $npos, length ($gate_name) - $npos);
					print "$i_count, $wires: ".$next.",   ";
					if  ($next =~ /^NOT/m){
						$qcircuit[$i_count++] .= $replace_NOT;
						$nfunct++;
						$npos += 3;
                                        } elsif ($next =~ /^V\+/m){
						$qcircuit[$i_count++] .= $replace_Vt;
						$nfunct++;
						$npos += 2;
                                        } elsif ($next =~ /^V/m){
						$qcircuit[$i_count++] .= $replace_V;
						$nfunct++;
						$npos++;
					} elsif ($next =~ /^t/m){
						$qcircuit[$i_count++] .= $replace_T;
						$nfunct++;
						$npos++;
					} elsif ($next =~ /^t\+/m){
						$qcircuit[$i_count++] .= $replace_Tt;
						$nfunct++;
						$npos++;
					} elsif ($next =~ /^H/m){
						$qcircuit[$i_count++] .= $replace_H;
						$nfunct++;
						$npos++;
                                        } elsif ($next =~ /^U/m){
						$qcircuit[$i_count++] .= $replace_U;
						$nfunct++;
						$npos++;
                                        } elsif ($next =~ /^Z/m){
						$qcircuit[$i_count++] .= $replace_Z;
						$nfunct++;
						$npos++;
                                        } elsif ($next =~ /^Y/m){
						$qcircuit[$i_count++] .= $replace_Y;
						$nfunct++;
						$npos++;
                                        } elsif ($next =~ /^CIIIII/m){
						$qcircuit[$i_count++] .= $replace_CIIIII;
						$npos++;
                                        } elsif ($next =~ /^CIIII/m){
						$qcircuit[$i_count++] .= $replace_CIIII;
						$npos++;
                                        } elsif ($next =~ /^CIII/m){
						$qcircuit[$i_count++] .= $replace_CIII;
						$npos++;
                                        } elsif ($next =~ /^CII/m){
						$qcircuit[$i_count++] .= $replace_CII;
						$npos++;
                                        } elsif ($next =~ /^CI/m){
						$qcircuit[$i_count++] .= $replace_CI;
						$npos++;
                                        } elsif ($next =~ /^C/m){
						if ($wires > 0){
							$wires++;
							$qcircuit[$i_count++] .= '& \ctrl{-'.$wires.'}';
							$wires = 0;
						} elsif ($nfunct > 0){
							$qcircuit[$i_count++] .= '& \ctrl{-1}';
							$wires = 0;
						} else {
							$qcircuit[$i_count++] .= $replace_C;
						}
						$npos++;
                                        } elsif ($next =~ /^I/m){
						$qcircuit[$i_count++] .= $replace_I;
						$npos++;
						if ($nfunct > 0){
							$wires++;
						}
					} else {
						#print "Did not match: $next \n";
						$npos++;
					}
				} 
			}
			$result .= ')';
			$octave_result .= ')';
			while ($kroneckers-- > 0){
				$octave_result .= ')';
			}
			print "\n";
		}
	}
	for ($u = 0; $u < $io-1; $u++){
		$qcircuit[$u] .= ' \\\\';
	}
	$qcircuit[$io-1] .= '}\]';

	print "\n\n";
	for ($u = 0; $u < 10; $u++){
		print $qcircuit[$u] ."\n\n";
	}
	print $result."\n";
	$octave_result =~ s/kron(\([^,]+\))/$1/g;
	print $octave_result."\n";

	
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
