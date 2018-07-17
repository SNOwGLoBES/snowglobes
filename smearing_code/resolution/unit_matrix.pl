#!/usr/bin/perl
# K. Scholberg July 2013
# Create resolution matrix for LAr

$resname = "unit_matrix";

$matsize = 200;
$emin = 0.5; # In MeV
$emax = 100;
$binsize = ($emax-$emin)/$matsize;

open(OUTFILE,">resolution_$resname.ssv");


# First fill the matrix, by column

for ($icol=0;$icol<$matsize;$icol++) {

  # Find the energy corresponding to the column

    $encol = $binsize*($icol+1);
    
    print $encol," ",$sigma,"\n";

    $sumrow = 0;
    for ($irow=0;$irow<$matsize;$irow++) {

# Find the energy corresponding to the row entry
      $enrow = $binsize*($irow+1);

# Gaussian with given sigma
      if ($irow==$icol) {
	  $resmat[$irow][$icol] = 1.;
      } else {
	  $resmat[$irow][$icol] = 0.;
      }
      $sumrow += $resmat[$irow][$icol];

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
