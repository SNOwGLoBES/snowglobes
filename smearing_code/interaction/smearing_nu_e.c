/* 
Originally from an old version of Rachel's smearing code

 Kate Scholberg Feb 2011
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef char * string;
int main(void) 
{

  string outfilename[6];
  outfilename[0]="interaction_nue_e.ssv";
  outfilename[1]="interaction_nuebar_e.ssv";
  outfilename[2]="interaction_numu_e.ssv";
  outfilename[3]="interaction_numubar_e.ssv";
  outfilename[4]="interaction_nutau_e.ssv";
  outfilename[5]="interaction_nutaubar_e.ssv";


 /* Specify energy bin parameters. */

  /* Energy in MeV */
 int number_bins=200;                  // number of energy bins
 double E_min=0.5;                     // first energy 
 double E_max=100.0;                   // max energy 

 double bin_size=(E_max-E_min)/number_bins;


 printf("bin_size %g\n",bin_size);
 /* Create y-distribution matrix. */
 
 double D[number_bins][number_bins];   // y-"D"istribution matrix
 int i, j;                             // i=column, j=row
 double E_e, E_nu;

 double norm, tmax;
 double me = 0.511;

 double s2tw = 0.231;

 int iflav;

 double gL[6];
 double gR[6];

 // For nue
 gL[0] = 2*s2tw+1;
 gR[0] = 2*s2tw;
 // For nuebar
 gL[1] = 2*s2tw;
 gR[1] = 2*s2tw+1;
 // For numu
 gL[2] = 2*s2tw-1;
 gR[2] = 2*s2tw;
 // For numubar
 gL[3] = 2*s2tw;
 gR[3] = 2*s2tw-1;
 // For nutau
 gL[4] = 2*s2tw-1;
 gR[4] = 2*s2tw;
 // For nutaubar
 gL[5] = 2*s2tw;
 gR[5] = 2*s2tw-1;


 for (iflav=0;iflav<6;iflav++) {
   E_nu=E_min+bin_size/2;

   for(j=0;j<number_bins;j++)            // step through columns
     {

       // Energy for this column j
       tmax = 2* pow(E_nu,2)/(me+2*E_nu);
     
       E_e = E_min+bin_size/2;

       for(i=0;i<number_bins;i++)        // step down a column, get image of the nu energy
	 {

	   norm = 1./(pow(gL[iflav],2)*tmax+pow(gR[iflav],2)*tmax-pow(gR[iflav],2)*pow(tmax,2)/E_nu
		      - gL[iflav]*gR[iflav]*me*pow(tmax,2)/(2*pow(E_nu,2))
		      + pow(gR[iflav],2)*pow(tmax,3)/(3*pow(E_nu,2)));


	 /* 	 if(i>=j)   */
	   if (E_e<=tmax)
	     D[i][j] = norm*(pow(gL[iflav],2)+pow(gR[iflav],2)*pow((1-E_e/E_nu),2)
			     -gR[iflav]*gL[iflav]*me*E_e/pow(E_nu,2))*bin_size;
	   else  D[i][j] = 0.0;

	   printf("%d %d %d %g %g %g %g\n",iflav,i,j,E_nu,E_e, tmax, D[i][j]);

	   E_e += bin_size;

	 }

       E_nu += bin_size;    

     }


   FILE* output=fopen(outfilename[iflav], "w");


   for(i=0;i<number_bins;i++)
     {
     
       for(j=0;j<number_bins;j++)
	 {

	   fprintf(output, " %g", D[i][j]);    
	    
	 }
     fprintf(output,"\n");    
     }

   fclose(output);                                                           // close output file


 }

 exit(0);

} 






