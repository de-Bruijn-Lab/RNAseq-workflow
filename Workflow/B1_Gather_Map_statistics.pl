use strict;
use warnings;

### 10/11/2016
# Input Mapping output files (The final.out versions)
# Read key mapping statistics
# Output to file 


### Setup. Create output file, print headers.

my $count = 1;

if (-e "./StarBAM/Log/outputMappingStats.txt") {
	system("rm ./StarBAM/Log/outputMappingStats.txt");
	print "Previous output removed";
}

open (my $ofh, ">>", "./StarBAM/Log/outputMappingStats.txt") 
	or die "Cannot open output : $!\n";

print $ofh "Filename\tTotal_Input_Reads\tAvg_read_length\tUniquely_Mapped_Reads\tAvg_Mapped_Length\tMismatch_rate_perBase\tMultimapped_reads\tMultimapped_tooManyLoci\tPercUnmapped_Mismatches\tPercUnmapped_tooShort\tPercUnmapped_other\n";


### Open output file, loop through file lines. Extract stats between lines 21 & 29 and print to file

foreach my $file (@ARGV) {
	open (my $fh, $file) 
		or die "Cannot open $file : $!\n";
		
	$count = 1;
	print $ofh $file."\t";
	
	while (<$fh>){		
		if ($count == 6 || $count == 7 || $count == 9 || $count == 11 || $count == 18 || $count == 24 || $count == 26 || $count == 29 || $count == 30 || $count == 31){
		
			# Regex to extract data
			$_ =~ s/[^\d.]//g;
			
			# Tab seperate entries
			if ($count < 31){
				$_ .= "\t";
				print $ofh $_;
			}
			
			# line 31 must end with new line
			if ($count == 31){
				$_ .= "\n";
				print $ofh $_;
				close($fh);
				last;
			}
		}
		$count++;
	}
	
	
}
close($ofh);
