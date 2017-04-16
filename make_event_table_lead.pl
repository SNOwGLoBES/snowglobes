#!/usr/bin/perl
# Script to run compile event rates from Globes output
# K. Scholberg Aug 2010
# arguments:  flux name, channel file name, experiment configuration name
# e.g.  ./make_event_table_lead.pl gvkm lead halo2

$fluxname = $ARGV[0];
$channame = $ARGV[1];
$expt_config = $ARGV[2];
#$nosmeared = $ARGV[3];
#Don't use smeared rates for lead
$nosmeared = 1;

# Now go through the channels 

$chanfilename = "channels/channels_".$channame.".dat";

open(CHANFILE,$chanfilename);

$tot_events=0;
$tot_es=0;
$tot_nc=0;
$tot_1n=0;
$tot_2n=0;
$tot_nc_nu_1n=0;
$tot_nc_nu_2n=0;
$tot_nc_nubar_1n=0;
$tot_nc_nubar_2n=0;

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
	    if  ($chan_name =~ /nc_nue_Pb208_1n/|| $chan_name =~ /nc_numu_Pb208_1n/||
                     $chan_name =~ /nc_nutau_Pb208_1n/) {
		$tot_nc_nu_1n += $numevents;
	    }

	    if  ($chan_name =~ /nc_nue_Pb208_2n/|| $chan_name =~ /nc_numu_Pb208_2n/ ||
                     $chan_name =~ /nc_nutau_Pb208_2n/) {
		$tot_nc_nu_2n += $numevents;
	    }

	    if  ($chan_name =~ /nc_nuebar_Pb208_1n/|| $chan_name =~ /nc_numubar_Pb208_1n/ || $chan_name =~ /nc_nutaubar_Pb208_1n/) {
		$tot_nc_nubar_1n += $numevents;
	    }
	    if  ($chan_name =~ /nc_nuebar_Pb208_2n/|| $chan_name =~ /nc_numubar_Pb208_2n/ || $chan_name =~ /nc_nutaubar_Pb208_2n/) {
		$tot_nc_nubar_2n += $numevents;
	    }
	    if ($chan_name =~ /1n/) {
		$tot_1n += $numevents;
	    }
	    if ($chan_name =~ /2n/) {
		$tot_2n += $numevents;
	    }


	}
    }

    close(OUTFILE);
#    print $output_line;

}

close(CHANFILE);

print "Total ES: ",$tot_es,"\n";
print "Total NC: ",$tot_nc,"\n";
print "Total 1n: ",$tot_1n,"\n";
print "Total 2n: ",$tot_2n,"\n";
print "Total NC nu 1n: ",$tot_nc_nu_1n,"\n";
print "Total NC nu 2n: ",$tot_nc_nu_2n,"\n";
print "Total NC nubar 1n: ",$tot_nc_nubar_1n,"\n";
print "Total NC nubar 2n: ",$tot_nc_nubar_2n,"\n";
print "Total neutron events: ",$tot_1n+$tot_2n,"\n";

print "Total events: ",$tot_events,"\n";


