#!/usr/bin/perl

#this code is designed to count the number of destabilizing mutations produced from a file 
#outputted by the calculaterobustness program

#the format should be an initial line with the initial refined energy
#thereafter all lines starting with a digit or negative symbol

open(M,$ARGV[0]) || die("could not open file");



@file1 = <M>;

close M;

foreach$line(@file1){if($line =~ /initial refined energy = (-.+)/)
                       {push @array, $1;
                       }
		     if($line =~ /robustness = 0\.+/)
			{print $line;
			}
                     if($line =~ /^-\d/ || $line =~ /^\d/)
                       {push @array, $line;
                       }
                     }

$size = scalar@array-1;

print "initial energy = ", $array[0], "\n";

$x = 0; $y = 0; $z = 0;

foreach$line1(@array){if($line1 < $array[0])
                        {$x = $x + 1;
                        }
                      if($line1 > $array[0])
                        {$y = $y + 1;
                        }
		      if($line1 ==  $array[0])
			{$z = $z + 1;
			}
                      }

$final = ($x / $size) * 100;
$final2 = ($y / $size) * 100;


print "this is the percentage that are stabilizing = ", $final, "\n";
print "this is the percentage that are destabilizing = ", $final2, "\n";
print "proportion of stabilizing to destabilizing = ", $final/$final2, "\n";
print "this is the number of neutral neighbours = ", $z, "\n";
