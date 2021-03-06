#!/usr/bin/perl

#This script will find the most common sequence in a list of DNA sequences
#and is integrated into the neutralrobustness program as a system call


#The argument is the input file

open(M,$ARGV[0]) || die("could not open file");

@file1=<M>;
	
$x=0;

foreach$line1(@file1){
                     if($line1 =~ /^[A,T,C,G]+/){
		     #if($line1 =~ /^[A,T,C,G]{168}/){
                                             $homologs[$x]=$line1;
                                             #print $homologs[$x]; print "\n";  
                                             $x = $x + 1; next;}
		      else{print "no sequence in file popndnaseqs, reported by mostcommonDNA.pl";}
                      }
		    


$numlines=scalar(@homologs);

# there should be a thousand lines
# the following code calculates the number of duplicates each line possesses

#print "this is the number of lines: ";
#rint $numlines; print"\n";

for($y=0;$y<$numlines-1;$y++){$duplicates[$y]=0; $z=1;
                           for($a=0;$a<($numlines-1)-$y;$a++){
                            if($homologs[$y] eq $homologs[$y+$z]){
                            
                            $duplicates[$y]=$duplicates[$y] + 1;
                               $z=$z + 1; next;}
                             
                             else{$z=$z + 1;next;}
                                           
                               }
                           #print$duplicates[$y]; print"\n";   
                            }


# the next lines will find line with the largest number of duplicates


$largest=0;

for($y=0;$y<$numlines-1;$y++){$z=1;
                             for($a=0;$a<($numlines - 1) -$y;$a++){
                             if($duplicates[$y] > $largest)
                               {$largest=$duplicates[$y]; 
                                #print"yes";$z=$z+1;
                                next;}

                            # else{$largest=$duplicates[$y+$z]; $z=$z+1;print"no";next;}
                                                     }
                              }

#print $largest + 1; print "\n";

# the next lines will extract the sequence number that corresponds to the highest
#number of duplicates


for($y=0;$y<$numlines-1;$y++){if($duplicates[$y]==$largest){
                              $number=$y; 
                              #print $number + 1; print "\n"; 
                              next;}
                             }


#print"\n";
print $homologs[$number];
