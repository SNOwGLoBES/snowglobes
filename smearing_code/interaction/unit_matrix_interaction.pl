#!/usr/bin/perl
# K. Scholberg July 2013


$matsize = 200;
$emin = 0.5; # In MeV
$emax = 100;
$binsize = ($emax-$emin)/$matsize;

open(OUTFILE,">interaction_unit_matrix.ssv");


# First fill the matrix

for ($icol=0;$icol<$matsize;$icol++) {
  for ($irow=0;$irow<$matsize;$irow++) {

# Find the energy corresponding to the row and column    
    $enrow = $binsize*($irow+1);
    $encol = $binsize*($icol+1);


    $edep = $encol;
    
    if (abs($enrow-$edep)<$binsize/2.) {
      $intmat[$irow][$icol] = 1.;
      
    } else {
      $intmat[$irow][$icol] = 0.;

    }

  }
}

# Now output the matrix, by row

for ($irow=0;$irow<$matsize;$irow++) {
  for ($icol=0;$icol<$matsize;$icol++) {


#    print $irow," ",$icol," ",$enrow," ",$encol," ",$intmat[$irow][$icol],"\n";

    print OUTFILE $intmat[$irow][$icol]," ";
    
  }
  print OUTFILE "\n";

}
#print OUTFILE "\n";

close(OUTFILE);
