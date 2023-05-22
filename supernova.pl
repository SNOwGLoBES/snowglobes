#!/usr/bin/perl
# Script to run supernova rates programs
# K. Scholberg July 2010
# arguments:  run mode, flux name, channel file name, experiment configuration name, noweight
# e.g.  ./supernova.pl 0 livermore water wc100kt30prct 0
# If noweight is set, then the output spectra will not have channel weighting factors applied  (default is that weighting factors are applied)

##
$run_choice = $ARGV[0];
$fluxname    = $ARGV[1];
$channame    = $ARGV[2];
$expt_config = $ARGV[3];
$noweight    = $ARGV[4];

$exename = "bin/supernova";
unless ( -f $exename ) {
    print $exename,
" executable not found.  Please compile and install, with SNOWGLOBES variable set.",
      "\n";
    exit;
}

$chanfilename = "channels/channels_" . $channame . ".dat";
$globesfilename = "supernova.glb";

# open( GLOBESFILE, ">$globesfilename" );
$energywindowset = "standard";

# This block of code will look for essential files
print
"\nBefore running SNOWGLOBES, it will check see if you have the corresponding files for each channel\n\n* Smearing matrix\n* Cross Section \n* Energy Effic\n\n";

# Open channel file
open( CHANFILE, $chanfilename );
while (<CHANFILE>) {

    # body...

    #split each line into an array
    if (/%/) { next; }

    $_ =~ s/\s*//;
    $_ =~ s/\s+/ /g;

    # if it has the SN_ID skip it
    @stuff = split( /\ /, $_ );

    if ( scalar(@stuff) == '12' ) {
        shift(@stuff);
    }

    $chan_name = $stuff[0];

    print( "\nChecking for: " . $chan_name . "\n" );
    $file_found = 1;

    # Formated names of each file type
    $xsc_file = "xscns/xs_" . $chan_name . ".dat";

    $smear_file = "smear/". $expt_config . "/smear_" . $chan_name . "_" . $expt_config . ".dat";

    $eff_file = "effic/". $expt_config . "/effic_" . $chan_name . "_" . $expt_config . ".dat";

    #File names

    # if not found post message
    unless ( -f $xsc_file ) {
        print "WARNING: XSC file: ", $xsc_file, " not found\n";
        $file_found = 0;
    }
    unless ( -f $smear_file ) {
        print "WARNING: Smear file: ", $smear_file, " not found\n";
    }
    unless ( -f $eff_file ) {
        print "WARNING: Effic file: ", $eff_file, " not found\n";
        $file_found = 0;
    }
    if ($file_found==1){
      print("All good !\n")
    }

}    # End of while loop

# Close chanfile
close(CHANFILE);
# Done looking for essential files

# #####################################
# ##### START of old supernova.pl #####
# #####################################


if ( $run_choice == '0' ) {

    # Time benchmark
    use Time::HiRes qw( time );
    my $start = time();

    $ewins{'bins'}            = "200";
    $ewins{'emin'}            = "0.0005";    #GeV
    $ewins{'emax'}            = "0.100";     #GeV
    $ewins{'sampling_points'} = "200";
    $ewins{'sampling_min'}    = "0.0005";    #GeV
    $ewins{'sampling_max'}    = "0.100";     #GeV

    $energywindowset = "standard";

    open( CHANFILE, $chanfilename );
    $firstchanline = <CHANFILE>;
    $firstchanline =~ s/^s+//;               # remove leading whitespace
    if ( index( $firstchanline, "%" ) != -1 ) {
        @binningarray = split /\s+/, $firstchanline;
        $ewins{'bins'}            = @binningarray[1];
        $ewins{'emin'}            = @binningarray[2];
        $ewins{'emax'}            = @binningarray[3];
        $ewins{'sampling_points'} = @binningarray[4];
        $ewins{'sampling_min'}    = @binningarray[5];
        $ewins{'sampling_max'}    = @binningarray[6];
    }
    close(CHANFILE);

    # these energy window/binning values are now ready for the preamble

# now that we've read in the channel file, we should know what energy windows we have
    if ( $ewins{'emax'} > 0.199 and $ewins{'emax'} < 0.201 ) {
        $energywindowset = "_he";
        print "using high energy window!\n";
    }

    # Create the globes file

    $globesfilename = "supernova.glb";

    open( GLOBESFILE, ">$globesfilename" );

# Here we add the globes preamble, with any modifications needed for rebinning taken care of
    open( PREAMBLE, "glb/preamble.glb" );
    $ready_to_modify = 0;
    while (<PREAMBLE>) {
        if ( index( $_, "Energy window" ) != -1 ) { $ready_to_modify = 1; }
        $ourline = $_;
        if ($ready_to_modify) {
            keys %ewins
              ; # reset the internal iterator so a prior each() doesn't affect the loop
            while ( my ( $k, $v ) = each %ewins ) {
                if ( index( $ourline, $k ) != -1 ) {
                    @vallist = split /\s+/, $ourline;
                    $ourline =~ s/$vallist[2]/$v/;
                    last;    # we have done the replacement, no need to continue
                }
            }    # end while loop over ewin replacements
        }    # any necessary modification is done
        print GLOBESFILE $ourline;
    }
    close(PREAMBLE);

    # Put in the flux info

    $fluxfilename = "fluxes/" . $fluxname . ".dat";
    unless ( -f $fluxfilename ) {
        print "Flux file name ", $fluxfilename, " not found\n";
        exit;
    }

    open( FLUX, "glb/flux.glb" );

    while (<FLUX>) {

        # Replace the flux file name with the input argument
        if (/flux_file/) {
            $_ = "        \@flux_file=  \"" . $fluxfilename . "\"\n";
        }
        print GLOBESFILE $_;
    }
    close(FLUX);

    # Now go through the channels and put in the relevant lines

    # First, smearing for each channel

    unless ( -f $chanfilename ) {
        print "Channel file name ", $chanfilename, " not found\n";
        exit;
    }
    open( CHANFILE, $chanfilename );

    while (<CHANFILE>) {

        # skip energy window info line
        if (/%/) { next; }

        # Grab the channel name
        $_ =~ s/\s*//;
        $_ =~ s/\s+/ /g;

        @stuff = split( /\ /, $_ );

        $chan_name = $stuff[0];

        $output_line =
            "include \"smear/". $expt_config . "/smear_"
          . $chan_name . "_"
          . $expt_config
          . ".dat\"\n";

        print GLOBESFILE $output_line;
    }

    close(CHANFILE);

    #  Detector info

    $detfilename = "detector_configurations.dat";

    unless ( -f $detfilename ) {
        print "Detector file name ", $detfilename, " not found\n";
        exit;
    }

    open( DETFILENAME, $detfilename );

    while (<DETFILENAME>) {
        chop($_);

        # Skip comments
        if (/\#/) { next; }
        $_ =~ s/\s*//;
        $_ =~ s/\s+/ /g;

        @stuff = split( /\ /, $_ );

        $detname            = $stuff[0];
        $masses{$detname}   = $stuff[1];
        $normfact{$detname} = $stuff[2];

        if (   $detname eq ""
            || $masses{$detname} eq ""
            || $normfact{$detname} eq "" )
        {
            next;
        }
    }
    close(DETFILENAME);

    # Mass and target normalization by species

#%masses = ("sk1",22.5,"sk2",22.5,"100kt15pc",100, "100kt30pc",100,"ar17",17,"scint50kt",50,"halo",1.0);
#%normfact = ("sk1",2/18,"sk2",2/18,"100kt15pc",2/18, "100kt30pc",2/18,"ar17kt",1/40,"scint50kt",2/14,"halo",1/208);

    if ( $masses{$expt_config} == 0 ) {
        print "Error: please enter a valid experiment configuration\n";
        exit;

    }

    $target_mass =
      sprintf( "%13.6f", $masses{$expt_config} * $normfact{$expt_config} )
      ;    # This is ktons of free protons

    print "Experiment config: ", $expt_config, "  Mass: ",
      $masses{$expt_config}, " kton \n";

    #print " Target mass: ",$normfact{$expt_config}," ",$target_mass,"\n";

    # Add the background smearing here, for the given detector configuration
    #  (Not yet implemented: for multiple background channels,
    # read them from a file labeled by detector configuration)
    $do_bg        = 0;
    $bg_chan_name = "bg_chan";

    $bg_filename = "backgrounds/" . $bg_chan_name . "_" . $expt_config . ".dat";
    if ( -e $bg_filename ) {
        $do_bg = 1;
        print "Using background file ", $bg_filename, "\n";
    }
    else {
        print "No background file for this configuration\n";
    }

    if ( $do_bg == 1 ) {
        $output_line =
            "include \"smear/". $expt_config . "/smear_"
          . $bg_chan_name . "_"
          . $expt_config
          . ".dat\"\n";
        print GLOBESFILE $output_line;
    }

    open( DETECTOR, "glb/detector.glb" );
    while (<DETECTOR>) {

        # Replace the flux file name with the input argument
        if (/mass/) {
            $_ = "\$target_mass=  " . $target_mass . "\n";
        }
        print GLOBESFILE $_;
    }
    close(DETECTOR);

    print GLOBESFILE "\n /******** Cross-sections *********/\n \n";

    # Now the cross-sections.  Note that some of these are repeated
    # even thougn it is not necessary (xscns for several flavors can be in the
    # same file).
    #  This is just to make a consistent loop over channels.

    open( CHANFILE, $chanfilename );

    while (<CHANFILE>) {

        # skip energy window info line
        if (/%/) { next; }

        # Grab the channel name

        $_ =~ s/\s*//;
        $_ =~ s/\s+/ /g;

        @stuff = split( /\ /, $_ );

        $chan_name = $stuff[0];

        $output_line = "cross(\#" . $chan_name . ")<\n";
        print GLOBESFILE $output_line;

        $output_line =
          "      \@cross_file= \"xscns/xs_" . $chan_name . ".dat\"\n";

        #    print $output_line;
        print GLOBESFILE $output_line;

        $output_line = "\>\n";

        #    print $output_line;
        print GLOBESFILE $output_line;

    }

    close(CHANFILE);

    # Now the fake bg channel cross section, if it exists

    if ( $do_bg == 1 ) {
        $output_line = "cross(\#" . $bg_chan_name . ")<\n";
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

    unless ( -f $chanfilename ) {
        print "Channel file name ", $chanfilename, " not found\n";
        exit;
    }

    open( CHANFILE, $chanfilename );

    while (<CHANFILE>) {

        # skip energy window info line
        if (/%/) { next; }

        # Grab the channel name

        $_ =~ s/\s*//;
        $_ =~ s/\s+/ /g;

        @stuff = split( /\ /, $_ );

        $chan_name = $stuff[0];
        $cpstate   = $stuff[2];
        $inflav    = $stuff[3];

        $output_line = "channel(\#" . $chan_name . "_signal)<\n";
        print GLOBESFILE $output_line;

        $output_line =
            "      \@channel= \#supernova_flux:  "
          . $cpstate . ":    "
          . $inflav
          . ":     "
          . $inflav
          . ":    \#"
          . $chan_name
          . ":    \#"
          . $chan_name
          . "_smear\n";

        #    print $output_line;
        print GLOBESFILE $output_line;

        # Get the post-smearing efficiencies by channel

        $eff_file =
            "effic/"
          . $expt_config
          . "/effic_"
          . $chan_name_post . "_"
          . $expt_config . ".dat";

        #    print $eff_file,"\n";
        open( EFF_FILE, $eff_file );
        while (<EFF_FILE>) {
            $output_line = "       \@post_smearing_efficiencies = " . $_;
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

    if ( $do_bg == 1 ) {

# Now make a fake channel for the background (just one possible bg file for now)

        # This is dummy info
        $cpstate = "-";
        $inflav  = "e";

        $output_line = "channel(\#" . $bg_chan_name . "_signal)<\n";
        print GLOBESFILE $output_line;

        $output_line =
            "      \@channel= \#supernova_flux:  "
          . $cpstate . ":    "
          . $inflav
          . ":     "
          . $inflav
          . ":    \#"
          . $bg_chan_name
          . ":    \#"
          . $bg_chan_name
          . "_smear\n";

        #    print $output_line;
        print GLOBESFILE $output_line;

        # Get the pre-smearing backgrounds by channel

        $bg_file = "backgrounds/" . $bg_chan_name . "_" . $expt_config . ".dat";
        print $bg_file, "\n";
        open( BG_FILE, $bg_file );
        while (<BG_FILE>) {
            $output_line = "       \@pre_smearing_background = " . $_;
            print GLOBESFILE $output_line;
        }
        close(BG_FILE);

# Try crazy reformatting
# or else get mysterious (but apparently harmless?) error from GLoBES
# Should try to track it down...  -> error goes away with development GLoBES version

        $output_line = "\n\>\n\n";
        print GLOBESFILE $output_line;
    }

    # End-matter
    if ( $energywindowset eq "standard" ) {
        open( POSTAMBLE, "glb/postamble.glb" );
    }
    elsif ( $energywindowset eq "_he" ) {
        open( POSTAMBLE, "glb/postamble_he.glb" );
    }

    while (<POSTAMBLE>) {
        print GLOBESFILE $_;
    }

    close(POSTAMBLE);

    close(GLOBESFILE);

    # Now run the executable

    $comstring =
      $exename . " " . $fluxname . " " . $chanfilename . " " . $expt_config;

    system($comstring);

    # Now apply channel-weighting factors-- this is the default

    unless ( $noweight == 1 ) {
        print "Applying channel weighting factors to output\n";

        &apply_weights("");
        &apply_weights("_smeared");

    }

    sub apply_weights {

        open( CHANFILE, $chanfilename );

        while (<CHANFILE>) {

            # skip energy window info line
            if (/%/) { next; }

            # Grab the channel name

            $_ =~ s/\s*//;
            $_ =~ s/\s+/ /g;

            @stuff = split( /\ /, $_ );

            $chan_name         = $stuff[0];
            $cpstate           = $stuff[2];
            $inflav            = $stuff[3];
            $num_target_factor = $stuff[4];

        # Open the unweighted output file to read and the weighted file to write

            $unweightfilename = "out/"
              . $fluxname . "_"
              . $chan_name . "_"
              . $expt_config
              . "_events"
              . $_[0]
              . "_unweighted.dat";

            $weightedfilename =
                "> out/"
              . $fluxname . "_"
              . $chan_name . "_"
              . $expt_config
              . "_events"
              . $_[0] . ".dat";

            open( UNWEIGHT, $unweightfilename );
            open( WEIGHTED, $weightedfilename );
            while (<UNWEIGHT>) {
                $_ =~ s/\s*//;
                $_ =~ s/\s+/ /g;

                if (/---/) {
                    print WEIGHTED $_, "\n";
                }
                else {
                    @stuff2 = split( /\ /, $_ );
                    $enbin  = $stuff2[0];
                    $evrate = $stuff2[1];
                    if ( $enbin ne "" ) {
                        print WEIGHTED $enbin, " ",
                          $evrate * $num_target_factor, "\n";
                    }
                }

            }

            close(UNWEIGHT);
            close(WEIGHTED);

        }

        close(CHANFILE);

    }
    print("All done !!\n");

    # Time benchmark
    my $end = time();
    printf( "Execution Time: %0.02f s\n", $end - $start );
}    #END of old supernova.pl

# ###################################
# #### START of new supernova.pl ####
# ###################################
# 08/10/20 Now supernova.pl will treat each interaction in the channel file as its own simulation
# This while loop will iterate through each interaction.
# the rest of the CHANFILE loops were deleted.
#  S. Torres-Lara
if ( $run_choice == '1' ) {

    # Time benchmark
    use Time::HiRes qw( time );
    my $start = time();
    $my_chan_num = 0;

    #Start of main while loop
    open( CHANFILE, $chanfilename );
    while (<CHANFILE>) {

        print "\n";

        "===================================================================";
        print "\n\n";
        @arr = split /\s+/, $_;
        shift(@arr);

        # skip energy window info line
        if (/%/) { next; }

        # Grab the channel name

        $_ =~ s/\s*//;
        $_ =~ s/\s+/ /g;

        @stuff = split( /\ /, $_ );

        # print @stuff,"\n";
        $general_id = $stuff[0];
        shift(@stuff);

        # print @stuff,"\n";

        $chan_name         = $stuff[0];
        $cpstate           = $stuff[2];
        $inflav            = $stuff[3];
        $num_target_factor = $stuff[4];

      # $chan_POSTAMBLE = sprintf( "postamble_%s.glb", $chan_name );
      # $path_to_POSTAMBLE = sprintf( "glb/%s/%s", $channame, $chan_POSTAMBLE );
        $GLOBAL_POSTAMBLE = "glb/postamble_SN.glb";

        ## If binning parameters are found use them
        if ( scalar(@arr) == '11' ) {
            print "| Using custom binning for: ", $arr[0], " |\n\n";
            print "_________________________________________________\n\n";
            # Flux parameters
            $ewins{'sampling_points'} = $arr[5];
            $ewins{'sampling_min'}    = $arr[6];     #GeV
            $ewins{'sampling_max'}    = $arr[7];    #GeV
            # Detected Energy parameters
            $ewins{'bins'}            = $arr[8];
            $ewins{'emin'}            = $arr[9];     #GeV
            $ewins{'emax'}            = $arr[10];     #GeV

            $energywindowset = "standard";

        }

        # If NO binning parameters are found use the standard set
        else {

            print "| Using standard binning for:", $arr[0], " |\n";
            print "_________________________________________________\n\n";

            $ewins{'bins'}            = "200";
            $ewins{'emin'}            = "0.0005";    #GeV
            $ewins{'emax'}            = "0.100";     #GeV
            $ewins{'sampling_points'} = "200";
            $ewins{'sampling_min'}    = "0.0005";    #GeV
            $ewins{'sampling_max'}    = "0.100";     #GeV

            $energywindowset = "standard";
        }

        if ( $ewins{'emax'} > 0.199 and $ewins{'emax'} < 0.201 ) {
            $energywindowset = "_he";
            print "using high energy window!\n";
        }

        open( GLOBESFILE, ">$globesfilename" );

# Here we add the globes preamble, with any modifications needed for rebinning taken care of
        open( PREAMBLE, "glb/preamble.glb" );
        $ready_to_modify = 0;
        while (<PREAMBLE>) {
            if ( index( $_, "Energy window" ) != -1 ) { $ready_to_modify = 1; }
            $ourline = $_;
            if ($ready_to_modify) {
                keys %ewins
                  ; # reset the internal iterator so a prior each() doesn't affect the loop
                while ( my ( $k, $v ) = each %ewins ) {
                    if ( index( $ourline, $k ) != -1 ) {
                        @vallist = split /\s+/, $ourline;
                        $ourline =~ s/$vallist[2]/$v/;
                        last
                          ;  # we have done the replacement, no need to continue
                    }
                }    # end while loop over ewin replacements
            }    # any necessary modification is done
            print GLOBESFILE $ourline;

        }
        close(PREAMBLE);

        $fluxfilename = "fluxes/" . $fluxname . ".dat";
        unless ( -f $fluxfilename ) {
            print "Flux file name ", $fluxfilename, " not found\n";
            exit;
        }

        open( FLUX, "glb/flux.glb" );

        while (<FLUX>) {

            # Replace the flux file name with the input argument
            if (/flux_file/) {
                $_ = "        \@flux_file=  \"" . $fluxfilename . "\"\n";
            }
            print GLOBESFILE $_;
        }
        close(FLUX);

        # Now go through the channels and put in the relevant lines

        # First, smearing for each channel

        unless ( -f $chanfilename ) {
            print "Channel file name ", $chanfilename, " not found\n";
            exit;
        }

        $output_line =
            "include \"smear/". $expt_config . "/smear_"
          . $chan_name . "_"
          . $expt_config
          . ".dat\"\n";

        print GLOBESFILE $output_line;

        #  Detector info

        $detfilename = "detector_configurations.dat";

        unless ( -f $detfilename ) {
            print "Detector file name ", $detfilename, " not found\n";
            exit;
        }

        open( DETFILENAME, $detfilename );

        while (<DETFILENAME>) {
            chop($_);

            # Skip comments
            if (/\#/) { next; }
            $_ =~ s/\s*//;
            $_ =~ s/\s+/ /g;

            @stuff_det = split( /\ /, $_ );

            $detname            = $stuff_det[0];
            $masses{$detname}   = $stuff_det[1];
            $normfact{$detname} = $stuff_det[2];

            if (   $detname eq ""
                || $masses{$detname} eq ""
                || $normfact{$detname} eq "" )
            {
                next;
            }
        }
        close(DETFILENAME);

        # Mass and target normalization by species

#%masses = ("sk1",22.5,"sk2",22.5,"100kt15pc",100, "100kt30pc",100,"ar17",17,"scint50kt",50,"halo",1.0);
#%normfact = ("sk1",2/18,"sk2",2/18,"100kt15pc",2/18, "100kt30pc",2/18,"ar17kt",1/40,"scint50kt",2/14,"halo",1/208);

        if ( $masses{$expt_config} == 0 ) {
            print "Error: please enter a valid experiment configuration\n";
            exit;

        }

        $target_mass =
          sprintf( "%13.6f", $masses{$expt_config} * $normfact{$expt_config} )
          ;    # This is ktons of free protons

        print "Experiment config: ", $expt_config, "  Mass: ",
          $masses{$expt_config},
          " kton \n";

        #print " Target mass: ",$normfact{$expt_config}," ",$target_mass,"\n";

        # Add the background smearing here, for the given detector configuration
        #  (Not yet implemented: for multiple background channels,
        # read them from a file labeled by detector configuration)
        $do_bg        = 0;
        $bg_chan_name = "bg_chan";

        $bg_filename =
          "backgrounds/" . $bg_chan_name . "_" . $expt_config . ".dat";
        if ( -e $bg_filename ) {
            $do_bg = 1;
            print "Using background file ", $bg_filename, "\n";
        }
        else {
            print "No background file for this configuration\n";
        }

        if ( $do_bg == 1 ) {
            $output_line =
                "include \"smear/". $expt_config . "/smear_"
              . $bg_chan_name . "_"
              . $expt_config
              . ".dat\"\n";
            print GLOBESFILE $output_line;
        }

        open( DETECTOR, "glb/detector.glb" );
        while (<DETECTOR>) {

            # Replace the flux file name with the input argument
            if (/mass/) {
                $_ = "\$target_mass=  " . $target_mass . "\n";
            }
            print GLOBESFILE $_;
        }
        close(DETECTOR);

        print GLOBESFILE "\n /******** Cross-sections *********/\n \n";

      # Now the cross-sections.  Note that some of these are repeated
      # even thougn it is not necessary (xscns for several flavors can be in the
      # same file).
      #  This is just to make a consistent loop over channels.

        $output_line = "cross(\#" . $general_id . ")<\n";
        print GLOBESFILE $output_line;

        $output_line =
          "      \@cross_file= \"xscns/xs_" . $chan_name . ".dat\"\n";

        #    print $output_line;
        print GLOBESFILE $output_line;

        $output_line = "\>\n";

        #    print $output_line;
        print GLOBESFILE $output_line;

        # Now the fake bg channel cross section, if it exists

        if ( $do_bg == 1 ) {
            $output_line = "cross(\#" . $bg_chan_name . ")<\n";
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

        unless ( -f $chanfilename ) {
            print "Channel file name ", $chanfilename, " not found\n";
            exit;
        }

        $chan_name_post = $stuff[0];
        $cpstate_post   = $stuff[2];
        $inflav_post    = $stuff[3];



        $output_line = "channel(\#" . $general_id . "_signal)<\n";
        print GLOBESFILE $output_line;

        $output_line =
            "      \@channel= \#supernova_flux:  "
          . $cpstate_post . ":    "
          . $inflav_post
          . ":     "
          . $inflav_post
          . ":    \#"
          . $general_id
          . ":    \#"
          . $chan_name
          . "_smear\n";

        print GLOBESFILE $output_line;

        # Get the post-smearing efficiencies by channel

        $eff_file =
            "effic/"
          . $expt_config
          . "/effic_"
          . $chan_name_post . "_"
          . $expt_config . ".dat";

        print $eff_file, "\n";
        open( EFF_FILE, $eff_file );
        while (<EFF_FILE>) {
            $output_line = "       \@post_smearing_efficiencies = " . $_;
            print GLOBESFILE $output_line;
        }
        close(EFF_FILE);

# Try crazy reformatting
# or else get mysterious (but apparently harmless?) error from GLoBES
# Should try to track it down...  -> error goes away with development GLoBES version

        $output_line = "\n\>\n\n";
        print GLOBESFILE $output_line;

        if ( $do_bg == 1 ) {

# Now make a fake channel for the background (just one possible bg file for now)

            # This is dummy info
            $cpstate = "-";
            $inflav  = "e";

            $output_line = "channel(\#" . $bg_chan_name . "_signal)<\n";
            print GLOBESFILE $output_line;

            $output_line =
                "      \@channel= \#supernova_flux:  "
              . $cpstate . ":    "
              . $inflav
              . ":     "
              . $inflav
              . ":    \#"
              . $bg_chan_name
              . ":    \#"
              . $bg_chan_name
              . "_smear\n";

            #    print $output_line;
            print GLOBESFILE $output_line;

            # Get the pre-smearing backgrounds by channel

            $bg_file =
              "backgrounds/" . $bg_chan_name . "_" . $expt_config . ".dat";
            print $bg_file, "\n";
            open( BG_FILE, $bg_file );
            while (<BG_FILE>) {
                $output_line = "       \@pre_smearing_background = " . $_;
                print GLOBESFILE $output_line;
            }
            close(BG_FILE);

# Try crazy reformatting
# or else get mysterious (but apparently harmless?) error from GLoBES
# Should try to track it down...  -> error goes away with development GLoBES version

            $output_line = "\n\>\n\n";
            print GLOBESFILE $output_line;
        }

        # End-matter
        if ( $energywindowset eq "standard" ) {
            open( POSTAMBLE, $GLOBAL_POSTAMBLE );
        }
        elsif ( $energywindowset eq "_he" ) {
            open( POSTAMBLE, $GLOBAL_POSTAMBLE );
        }

        while (<POSTAMBLE>) {
            print GLOBESFILE $_;
        }

        close(POSTAMBLE);

        close(GLOBESFILE);

        # Now run the executable

        $comstring =
            $exename . " "
          . $fluxname . " "
          . $chanfilename . " "
          . $expt_config . " "
          . $my_chan_num;

        system($comstring);

        # Now apply channel-weighting factors-- this is the default

        unless ( $noweight == 1 ) {
            print "Applying channel weighting factors to output\n";

            &apply_weights("");
            &apply_weights("_smeared");

        }

        sub apply_weights {

        # Open the unweighted output file to read and the weighted file to write
        # print  $_[0],"\n";
            $unweightfilename = "out/"
              . $fluxname . "_"
              . $chan_name . "_"
              . $expt_config
              . "_events"
              . $_[0]
              . "_unweighted.dat";

            $weightedfilename =
                "> out/"
              . $fluxname . "_"
              . $chan_name . "_"
              . $expt_config
              . "_events"
              . $_[0] . ".dat";

            open( UNWEIGHT, $unweightfilename );
            open( WEIGHTED, $weightedfilename );
            while (<UNWEIGHT>) {
                $_ =~ s/\s*//;
                $_ =~ s/\s+/ /g;

                if (/---/) {
                    print WEIGHTED $_, "\n";
                }
                else {
                    @stuff2 = split( /\ /, $_ );
                    $enbin  = $stuff2[0];
                    $evrate = $stuff2[1];

                    if ( $enbin ne "" ) {
                        print WEIGHTED $enbin, " ",
                          $evrate * $num_target_factor,
                          "\n";

                    }
                }

            }

            close(UNWEIGHT);
            close(WEIGHTED);

        }

        print "\n";
        $my_chan_num++;
    }
    close(CHANFILE);

    print("All done !!\n");

    # Time benchmark
    my $end = time();
    printf( "Execution Time: %0.02f s\n", $end - $start );
    exit;

}    #END of custom binning
#
#
#
#
