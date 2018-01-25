use strict;
use warnings;

### 10/11/2016
### Input FastQ Trim output files
### Read key trimming statistics
### Output to file 


###
### Setup. Create output file, print headers.
###

my $count = 1;

if (-e "./trim/trimOutput/outputTrimStats.txt") {
	system("rm ./trim/trimOutput/outputTrimStats.txt");
	print "Previous output removed";
}
open (my $ofh, ">>", "./trim/trimOutput/outputTrimStats.txt") 
	or die "Cannot open output : $!\n";

print $ofh "ProcessedReads\tProcessedBases\tTrimmedReads\tQuality-trimmed\tTrimmedBases\tTooShortReads\tTooLongReads\tTotalTime\tTimePerRead\n";


###
### Open Sample.txt, loop through file lines. Extract stats between lines 21 & 29 and print to file
###

foreach my $file (@ARGV) {
	open (my $fh, $file) 
		or die "Cannot open $file : $!\n";
	$count = 1;
	while (<$fh>){
		if ($count > 29){
			close($fh);
			last;
		}
		
		if ($count >= 21){
		
			# Regex to extract data
			$_ =~ s/ //g;
			$_ =~ s/[(].*[)]//g;
			$_ =~ s/bp$|ms$|s$|m$//g;
			my @heads = $_ =~ /(^[^:]*:)/g;
			$_ =~ s/^[^:]*://g;
			$_ =~ s/\n//g;
			
			# Print file name
			if ($count == 21){
				print $ofh $file."\t";
			}
			
			# Tab seperate entries
			if ($count < 29){
				$_ .= "\t";
				print $ofh $_;
			}
			
			# line 29 must end with new line
			if ($count == 29){
				$_ .= "\n";
				print $ofh $_;
			}
		}
		$count++;
	}
	close($fh);
	
}
close($ofh);
