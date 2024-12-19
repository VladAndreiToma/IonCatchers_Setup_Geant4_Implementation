#!/usr/bin/perl

# SETS UP AND LAUNCHES THE NUMBER OF JOBS N GIVEN AS FIRST PARAMETER
# SECOND PARAMETER IS THE RUN MODE IN RunInput.txt: 0 or 1

use strict;

if(!defined($ARGV[0])) {die("Please give me the number of jobs!\n");}
my $nJobs  =  $ARGV[0];
if(!defined($ARGV[1])) {die("Please give me the run step!\n");}
my $runStep  =  $ARGV[1];

my $basePath = $ENV{IonCatchersPath};
my $baseDir = $basePath.'/build/';
for(my $iJob=0; $iJob<$nJobs; $iJob++){
    my $dirJob = $baseDir.'run'.$iJob.'/';
    unless(-d $dirJob) {                                        # CREATE JOB DIRECTORY
	mkdir($dirJob);                                           # IF IT DOES NOT EXIST
    }
    chdir $dirJob;                                                  # GO TO JOB DIRECTORY
    if($runStep == 0) {
	`rm -f *`;                                                     # CLEAN UP
	`cp ../IonCatchers .`;                                 # COPY EXECUTABLE
	`./IonCatchers >Output_Mode0.txt 2>Error_Mode0.txt&`;   # LAUNCH THE JOB
	print("Running ./IonCatchers >Output_Mode0.txt 2>Error_Mode0.txt in:\n");
	system 'pwd'; 
	`sleep 180`;                                               # WAIT 3 MINUTES, IMPORTANT FOR RANDOM SEED!
    } elsif ($runStep == 1) {
	`./IonCatchers >Output_Mode1.txt 2>Error_Mode1.txt&`;   # LAUNCH THE JOB
	print("Running ./IonCatchers >Output_Mode1.txt 2>Error_Mode1.txt in:\n");
	system 'pwd'; 
	`sleep 10`;                                                 # WAIT 10 SECONDS	
    } else {die("Unknown run step!");}
}
