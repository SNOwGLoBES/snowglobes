/* GLoBES -- General LOng Baseline Experiment Simulator
 * (C) 2002 - 2007,  The GLoBES Team
 *
 * GLoBES is mainly intended for academic purposes. Proper
 * credit must be given if you use GLoBES or parts of it. Please
 * read the section 'Credit' in the README file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdio.h>
#include <string.h>
#include "myio.h"

static char* THEFILE;

void InitOutput(char* filename, char* headline)
{
 THEFILE=filename;
 if(strlen(THEFILE)==0) printf(headline);
 else 
 {
   FILE* f=fopen(THEFILE, "w");
   if (!f)
   {
     printf("File cannot be opened!\n");
     THEFILE[0]=0;
   }
   else {
    fprintf(f,headline);
    fclose(f);
   }
 }
}

FILE *InitDiffSpecOutput(char* filename)
{
   /*    FILE* f=fopen(THEFILE, "w"); */
  FILE* f=fopen(filename, "w"); 
  if (!f)
    {
      printf("File %s cannot be opened!\n",filename);
    }
  else {
    return f;
  }
 
}


void AddToOutput(double n1,double n2,double n3)
{
 if(strlen(THEFILE)==0) printf("%g %g %f \n",n1,n2,n3);
 else 
 {
   FILE* f=fopen(THEFILE, "a");
   if (!f)
   {
     printf("File cannot be opened!\n");
     THEFILE[0]=0;
   }
   else
   {
    fprintf(f,"%g %g %f \n",n1,n2,n3);
    fclose(f);
   }
 }
}

void AddToOutput2(double n1,double n2)
{
 if(strlen(THEFILE)==0) printf("%g %g \n",n1,n2);
 else 
 {
   FILE* f=fopen(THEFILE, "a");
   if (!f)
   {
     printf("File cannot be opened!\n");
     THEFILE[0]=0;
   }
   else
   {
    fprintf(f,"%g %g \n",n1,n2);
    fclose(f);
   }
 }
}
