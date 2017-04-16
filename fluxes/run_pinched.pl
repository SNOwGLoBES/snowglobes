#!/usr/bin/perl  
# Run jobs to create grid of smeared spectra from pinched flux files
use File::Copy;
#use strict;
#use warnings;

$influxdir = $ARGV[0];
$outdir = $ARGV[1];

if ($influxdir eq "" || $outdir eq "") {
    print "Usage: ./run_pinched.pl influxdir outdir\n"; 
    exit();
}

foreach $fluxfile (`ls $influxdir`)
{
    print $fluxfile;
    chop($fluxfile);

# Copy flux file to the SNOwGLoBES fluxes directory so SNOwGLoBES knows about it

    $fromfile = $influxdir."/".$fluxfile;
    $tofile = "fluxes/".$fluxfile;

    copy($fromfile,$tofile);
    
    $fluxfile =~ s/\.dat//;

    $command = "./supernova.pl $fluxfile argon ar17kt";

    print $command,"\n";

# Run SNOwGLoBES for this flux

    system $command;

# Clean up the copy of the flux

    unlink $tofile;
    
# Now sum the relevant SNOwGLoBES output files and put in an output file

# Initialize
    $numbins = 200;
    for ($j=0;$j<$numbins;$j++) {
	$energy[$j]=0.;
	$totevents[$j]=0.;
    }

# Add up all the smeared files for this flux
    foreach $sgfile (`ls out/$fluxfile*_smeared.dat`) {

	chop($sgfile);

	print "Opening ",$sgfile,"\n";
	open(SMEARFILE,$sgfile);

	$i=0;
	while(<SMEARFILE>) {

	    $_=~s/\s*//;
	    $_=~s/\s+/ /g;

	    @stuff = split(/\ /,$_);

	    $energy[$i]=$stuff[0];
	    $totevents[$i]+=$stuff[1];
#	    print $i," ",$energy[$i]," ",$totevents[$i],"\n";
	    $i++;
	    
	}

	close(SMEARFILE);
	
    }  # End of loop over smeared files for this flux

# Now put the sum in the output file
    $sumfile = ">".$outdir."/".$fluxfile."_smeared_sum.dat";
    print "Output file: ",$sumfile,"\n";
    open(SUMFILE,$sumfile);
 
   
    for ($j=0;$j<$numbins;$j++) {

	$tot = sprintf("%e",$totevents[$j]);
	print SUMFILE $energy[$j]," ",$tot,"\n";
    }


    close(SUMFILE);



# Clean up all the output

    foreach $sgfile (`ls out/$fluxfile*`) {
	print "Removing Snowglobes file: ",$sgfile;

	chop($sgfile);
	unlink $sgfile or warn "could not unlink!";
	
    }




} # End of loop over flux files
