#!/usr/bin/perl
# Script to compile event rates from Globes output
# K. Scholberg Aug 2010
# arguments:  flux name, channel file name, experiment configuration name
# e.g.  ./make_event_table.pl livermore argon ar17

$fluxname = $ARGV[0];
$channame = $ARGV[1];
$expt_config = $ARGV[2];
$nosmeared = $ARGV[3];

# Now go through the channels 

$chanfilename = "channels/channels_".$channame.".dat";

open(CHANFILE,$chanfilename);

$tot_events=0;
$tot_es=0;
$tot_nc=0;
while(<CHANFILE>) {

# Grab the channel name

    $_=~s/\s*//;
    $_=~s/\s+/ /g;

    @stuff = split(/\ /,$_);

    $chan_name = $stuff[0];

    if ($nosmeared == 1) {
	$outfile = "out/".$fluxname."_".$chan_name."_".$expt_config."_events.dat";
    } else {
	$outfile = "out/".$fluxname."_".$chan_name."_".$expt_config."_events_smeared.dat";

    }

# Target weighting factor
#    $target_fact = $stuff[4];

# Target factor only needed for unweighted
    $target_fact = 1;


#    print $outfile,"\n";
    open (OUTFILE,$outfile);

    while (<OUTFILE>) {

	if (/Total/) {
	    $_=~s/\s*//;
	    $_=~s/\s+/ /g;

	    @stuff2 = split(/\ /,$_);
	    
	    $numevents = $stuff2[1]*$target_fact;
	    print $chan_name," ",$numevents,"\n";
	    $tot_events += $numevents;
	    if ($chan_name =~ /\_e/) {
		$tot_es += $numevents;
	    }
	    if ($chan_name =~ /nc/) {
		$tot_nc += $numevents;
	    }


	}
    }

    close(OUTFILE);
#    print $output_line;

}

close(CHANFILE);

print "Total ES: ",sprintf("%12.3f \n",$tot_es);
print "Total NC: ",$tot_nc,"\n";
print "Total events: ",$tot_events,"\n";


