#!/usr/bin/perl

# runs analysis/AnalysisTree.C macro in each runX directory to make the Tree ROOT file
# runs analysis/MergeOutput.C macro to merge ROOT & TXT files in final output

use strict;

if(!defined($ARGV[0])) {die("Please give me the number of jobs!\n");}
my $nJobs  =  $ARGV[0];

my $basePath = $ENV{IonCatchersPath};
my $baseDir = $basePath.'/test/';
my $execMac = $basePath.'/analysis/AnalysisTree.C';
my $mergMac = $basePath.'/analysis/MergeOutput.C';
for(my $iJob=0; $iJob<$nJobs; $iJob++){
    my $dirJob = $baseDir.'run'.$iJob.'/';
    chdir $dirJob;                                                     # GO TO JOB DIRECTORY
    system 'pwd';
    `root -b $execMac > Output.txt`;                   # EXECUTE ANALYSIS MACRO
    print("root -b $execMac > Output.txt\n");
    `sleep 60`;                                                        # WAIT 1 MINUTE
}
chdir $baseDir;                                                     # GO TO BASE DIRECTORY
system 'pwd';
`root -b $mergMac`;                                             # EXECUTE MERGE MACRO
print("root -b $mergMac\n");
