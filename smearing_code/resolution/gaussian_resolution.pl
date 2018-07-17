#!/usr/bin/perl
# K. Scholberg July 2013
# Create resolution matrix e.g. ./gaussian_resolution.pl 0.18

$resolution = $ARGV[0];

$matsize = 200;
$emin = 0.5; # In MeV
$emax = 100;
$binsize = ($emax-$emin)/$matsize;

open(OUTFILE,">resolution_$resolution.ssv");

print "Bin size ",$binsize,"\n";

# First fill the matrix, by column

for ($icol=0;$icol<$matsize;$icol++) {

  # Find the energy corresponding to the column

    $encol = $binsize*($icol+1);
    
    # Get the sigma for the column energy

# 2/9/14: factor of 2
#    $sigma = 2.*(.11/sqrt($encol)+0.02)*$encol;
# 2/9/14: factor of 2
    $sigma = $encol*$resolution;

    print $encol," ",$sigma,"\n";

    $sumrow = 0;
    for ($irow=0;$irow<$matsize;$irow++) {

# Find the energy corresponding to the row entry
      $enrow = $binsize*($irow+1);

# Gaussian with given sigma
      $gaussval = exp(-($enrow-$encol)**2/$sigma**2); 
      $sumrow += $gaussval;
      $resmat[$irow][$icol] = $gaussval;


    }  # End of loop over row entries for this column

# Normalize the column 

    for ($irow=0;$irow<$matsize;$irow++) {
      $resmat[$irow][$icol] /= $sumrow; 

    }

}  # End of loop over columns

# Now output the matrix, by row


for ($irow=0;$irow<$matsize;$irow++) {
  for ($icol=0;$icol<$matsize;$icol++) {

    $printval = sprintf("%8.4f",$resmat[$irow][$icol]);
    print OUTFILE $printval," ";
    
  }
  print OUTFILE "\n";

}
#print OUTFILE "\n";

close(OUTFILE);
