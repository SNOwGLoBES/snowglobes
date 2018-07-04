#!/usr/bin/perl
# Script to run supernova rates programs
# K. Scholberg July 2010
# arguments:  flux name, channel file name, experiment configuration name, noweight
# e.g.  ./supernova.pl livermore water wc100kt30prct 0
# If noweight is set, then the output spectra will not have channel weighting factors applied  (default is that weighting factors are applied)


$fluxname = $ARGV[0];
$channame = $ARGV[1];
$expt_config = $ARGV[2];
$noweight = $ARGV[3];


$exename = "bin/supernova";
unless (-f $exename)  {
    print $exename," executable not found.  Please compile and install, with SNOWGLOBES variable set.", "\n";
    exit;
}


$chanfilename = "channels/channels_".$channame.".dat";

# Create the globes file

$globesfilename = "supernova.glb";

open(GLOBESFILE,">$globesfilename");


open(PREAMBLE,"glb/preamble.glb");
while(<PREAMBLE>) {
    print GLOBESFILE $_;
}
close(PREAMBLE);

# Put in the flux info

$fluxfilename = "fluxes/".$fluxname.".dat";
unless(-f $fluxfilename) {
    print "Flux file name ",$fluxfilename," not found\n";
    exit;
}


open(FLUX,"glb/flux.glb");

while(<FLUX>) {
# Replace the flux file name with the input argument
    if (/flux_file/) {
	$_ = "        \@flux_file=  \"".$fluxfilename."\"\n";
    }
    print GLOBESFILE $_;
}
close(FLUX);

# Now go through the channels and put in the relevant lines

# First, smearing for each channel

unless (-f $chanfilename) {
    print "Channel file name ",$chanfilename," not found\n";
    exit;
}
open(CHANFILE,$chanfilename);

while(<CHANFILE>) {

# Grab the channel name

    $_=~s/\s*//;
    $_=~s/\s+/ /g;

    @stuff = split(/\ /,$_);

    $chan_name = $stuff[0];

    $output_line = "include \"smear/smear_".$chan_name."_".$expt_config.".dat\"\n";


#    print $output_line;
    print GLOBESFILE $output_line;
}

close(CHANFILE);

#  Detector info

$detfilename = "detector_configurations.dat";

unless (-f $detfilename) {
    print "Detector file name ",$detfilename," not found\n";
    exit;
}

open(DETFILENAME,$detfilename);

while(<DETFILENAME>) {
    chop($_);
    
# Skip comments
    if (/\#/) {next;}
    $_=~s/\s*//;
    $_=~s/\s+/ /g;

    @stuff = split(/\ /,$_);

    $detname = $stuff[0];
    $masses{$detname} = $stuff[1];
    $normfact{$detname} = $stuff[2];

    if ($detname eq "" || $masses{$detname} eq "" || $normfact{$detname} eq ""){ next;}
}
close(DETFILENAME);

# Mass and target normalization by species

#%masses = ("sk1",22.5,"sk2",22.5,"100kt15pc",100, "100kt30pc",100,"ar17",17,"scint50kt",50,"halo",1.0);
#%normfact = ("sk1",2/18,"sk2",2/18,"100kt15pc",2/18, "100kt30pc",2/18,"ar17kt",1/40,"scint50kt",2/14,"halo",1/208);

if ($masses{$expt_config} == 0){
    print "Error: please enter a valid experiment configuration\n";
    exit;

}

$target_mass = sprintf("%13.6f",$masses{$expt_config}*$normfact{$expt_config});  # This is ktons of free protons

print "Experiment config: ",$expt_config,"  Mass: ",$masses{$expt_config}," kton \n";
#print " Target mass: ",$normfact{$expt_config}," ",$target_mass,"\n";

# Add the background smearing here, for the given detector configuration
#  (Not yet implemented: for multiple background channels, 
# read them from a file labeled by detector configuration)
$do_bg = 0;
$bg_chan_name = "bg_chan";

$bg_filename = "backgrounds/".$bg_chan_name."_".$expt_config.".dat";
if (-e $bg_filename) {
    $do_bg = 1;
    print "Using background file ",$bg_filename,"\n";
} else {
    print "No background file for this configuration\n";
}

if ($do_bg == 1) {
    $output_line = "include \"smear/smear_".$bg_chan_name."_".$expt_config.".dat\"\n";
    print GLOBESFILE $output_line;
}

open(DETECTOR,"glb/detector.glb");
while(<DETECTOR>) {
# Replace the flux file name with the input argument
    if (/mass/) {
	$_ = "\$target_mass=  ".$target_mass."\n";
    }
    print GLOBESFILE $_;
}
close(DETECTOR);


print GLOBESFILE "\n /******** Cross-sections *********/\n \n";

# Now the cross-sections.  Note that some of these are repeated 
# even thougn it is not necessary (xscns for several flavors can be in the 
# same file).
#  This is just to make a consistent loop over channels.

open(CHANFILE,$chanfilename);

while(<CHANFILE>) {

# Grab the channel name

    $_=~s/\s*//;
    $_=~s/\s+/ /g;

    @stuff = split(/\ /,$_);

    $chan_name = $stuff[0];

    $output_line = "cross(\#".$chan_name.")<\n";
    print GLOBESFILE $output_line;

    $output_line = "      \@cross_file= \"xscns/xs_".$chan_name.".dat\"\n";
#    print $output_line;
    print GLOBESFILE $output_line;

    $output_line = "\>\n";
#    print $output_line;
    print GLOBESFILE $output_line;

}

close(CHANFILE);

# Now the fake bg channel cross section, if it exists

if ($do_bg == 1) {
    $output_line = "cross(\#".$bg_chan_name.")<\n";
    print GLOBESFILE $output_line;

    $output_line = "      \@cross_file= \"xscns/xs_zero.dat\"\n";
#    print $output_line;
    print GLOBESFILE $output_line;
    
    $output_line = "\>\n";
#    print $output_line;
    print GLOBESFILE $output_line;
}

print GLOBESFILE "\n \/******** Channels *********\/\n \n";

# Now, the channel definitions

unless(-f $chanfilename) {
    print "Channel file name ",$chanfilename," not found\n";
    exit;
}

open(CHANFILE,$chanfilename);

while(<CHANFILE>) {

# Grab the channel name

    $_=~s/\s*//;
    $_=~s/\s+/ /g;

    @stuff = split(/\ /,$_);

    $chan_name = $stuff[0];
    $cpstate = $stuff[2];
    $inflav = $stuff[3];

    $output_line = "channel(\#".$chan_name."_signal)<\n";
    print GLOBESFILE $output_line;

    $output_line = "      \@channel= \#supernova_flux:  ".$cpstate.":    ".$inflav.":     ".$inflav.":    \#".$chan_name.":    \#".$chan_name."_smear\n"; 
#    print $output_line;
    print GLOBESFILE $output_line;

# Get the post-smearing efficiencies by channel

    $eff_file = "effic/effic_".$chan_name."_".$expt_config.".dat";
#    print $eff_file,"\n";
    open(EFF_FILE,$eff_file);
    while(<EFF_FILE>) {
	$output_line = "       \@post_smearing_efficiencies = ".$_;
	print GLOBESFILE $output_line;
    }
    close(EFF_FILE);

# Try crazy reformatting
# or else get mysterious (but apparently harmless?) error from GLoBES
# Should try to track it down...  -> error goes away with development GLoBES version

    $output_line = "\n\>\n\n";
    print GLOBESFILE $output_line;

}

close(CHANFILE);

if ($do_bg == 1) {
# Now make a fake channel for the background (just one possible bg file for now)

# This is dummy info
    $cpstate = "-"; $inflav = "e";

    $output_line = "channel(\#".$bg_chan_name."_signal)<\n";
    print GLOBESFILE $output_line;

    $output_line = "      \@channel= \#supernova_flux:  ".$cpstate.":    ".$inflav.":     ".$inflav.":    \#".$bg_chan_name.":    \#".$bg_chan_name."_smear\n"; 
#    print $output_line;
    print GLOBESFILE $output_line;

# Get the pre-smearing backgrounds by channel

    $bg_file = "backgrounds/".$bg_chan_name."_".$expt_config.".dat";
    print $bg_file,"\n";
    open(BG_FILE,$bg_file);
    while(<BG_FILE>) {
	$output_line = "       \@pre_smearing_background = ".$_;
	print GLOBESFILE $output_line;
    }
    close (BG_FILE);

# Try crazy reformatting
# or else get mysterious (but apparently harmless?) error from GLoBES
# Should try to track it down...  -> error goes away with development GLoBES version

$output_line = "\n\>\n\n";
print GLOBESFILE $output_line;
}

# End-matter

open(POSTAMBLE,"glb/postamble.glb");

while(<POSTAMBLE>) {
    print GLOBESFILE $_;
}


close(POSTAMBLE);

close(GLOBESFILE);


# Now run the executable

$comstring = $exename." ".$fluxname." ".$chanfilename." ".$expt_config;

system($comstring);

# Now apply channel-weighting factors-- this is the default

unless($noweight==1) {
 print "Applying channel weighting factors to output\n";


&apply_weights("");
&apply_weights("_smeared");

}



sub apply_weights
  {

open(CHANFILE,$chanfilename);

while(<CHANFILE>) {

# Grab the channel name

    $_=~s/\s*//;
    $_=~s/\s+/ /g;

    @stuff = split(/\ /,$_);

    $chan_name = $stuff[0];
    $cpstate = $stuff[2];
    $inflav = $stuff[3];
    $num_target_factor = $stuff[4];

# Open the unweighted output file to read and the weighted file to write

    $unweightfilename = "out/".$fluxname."_".$chan_name."_".$expt_config."_events".$_[0]."_unweighted.dat";

    $weightedfilename = "> out/".$fluxname."_".$chan_name."_".$expt_config."_events".$_[0].".dat";

    open (UNWEIGHT,$unweightfilename);
    open (WEIGHTED,$weightedfilename);
    while(<UNWEIGHT>) {
      $_=~s/\s*//;
      $_=~s/\s+/ /g;
     
      if (/---/) {
	print WEIGHTED $_,"\n";
      } else {
	@stuff2 = split(/\ /,$_);
	$enbin = $stuff2[0];
	$evrate = $stuff2[1];
	if ($enbin ne "") {
	  print WEIGHTED $enbin," ",$evrate*$num_target_factor,"\n";
	}
      }


    }

    close(UNWEIGHT);
    close(WEIGHTED);


  }

 close(CHANFILE);




  }

