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

/*
 * Example: Correlation between s22th13 and deltacp
 * Compile with ``make example1''
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include <globes/globes.h>   /* GLoBES library */
#include "myio.h"             /* my input-output routines */

/* If filename given, write to file; for empty filename write to screen */
char MYFILE[]="supernova.dat";

int ret;

char flux_file_name[256];
char channel_file_name[256];
char expt_config_name[256];

static int maxchans = 32;

char outfile[256];
char outfile_smeared[256];
char bgfile[256];




struct stat buf;

int main(int argc, char *argv[])
{
        /* Initialize libglobes */
        glbInit(argv[0]);

        if (argc<2) {
                printf("Arguments required: flux filename, channels filename, detector configuration name\n");
                exit(0);
        }

        strncpy(flux_file_name,argv[1],256);

        /* Read the channels to process from the file */
        strncpy(channel_file_name,argv[2],256);
        printf("Channels from %s\n",channel_file_name);

        strncpy(expt_config_name,argv[3],256);

        FILE* fp_chans = fopen(channel_file_name,"r");

        if (fp_chans == NULL) {
                printf ("Cannot open file \n");
                exit(0);
        }

        int chan_num[maxchans];
        char chan_name[maxchans][256];
        char chbuf[1000];

        int num_chans=0;
        char gen_id[256];
        char cp;
        char flav;
        int num_target_factor=1;
        int skipline = 0;
        char* eofcheck; /* need to check return of fgets for whether line was read successfully */
        // Set propper name formating for custom binning
        if(argc==5) {
                while (!feof(fp_chans)) {
                        skipline = 0;
                        eofcheck = fgets(chbuf,1000,fp_chans);
                        for (int ichar = 0; ichar != 5; ichar++) {
                                if (chbuf[ichar] == '%') {
                                        skipline = 1;
                                }
                        }
                        if (skipline) {continue;}
                        if (eofcheck != NULL) {
                                ret = sscanf(chbuf,"%s %s %i %s %s %i",gen_id, chan_name[num_chans],&chan_num[num_chans],&cp, &flav, &num_target_factor);
                                num_chans++;
                        }
                }

                fclose(fp_chans);
        }
          // Set propper name formating for pre-set binning
        if (argc==4) {
                while (!feof(fp_chans)) {
                        skipline = 0;
                        eofcheck = fgets(chbuf,1000,fp_chans);
                        for (int ichar = 0; ichar != 5; ichar++) {
                                if (chbuf[ichar] == '%') {
                                        skipline = 1;
                                }
                        }
                        if (skipline) {continue;}
                        if (eofcheck != NULL) {
                                ret = sscanf(chbuf,"%s %i %s %s %i",chan_name[num_chans],&chan_num[num_chans],&cp, &flav, &num_target_factor); /*      printf("%s %i\n",chan_name[num_chans],chan_num[num_chans]); */
                                num_chans++;
                        }
                }

                fclose(fp_chans);

        }


        // printf("Number of channels found: %i\n",num_chans);


        /* Initialize experiment NFstandard.glb */
        glbInitExperiment("supernova.glb",&glb_experiment_list[0],&glb_num_of_exps);


        /* Intitialize output */
        /*  InitOutput(MYFILE,"Format: Log(10,s22th13)   deltacp   chi^2 \n"); */

        /* Define standard oscillation parameters */
        /*  double theta12 = asin(sqrt(0.8))/2;
           double theta13 = asin(sqrt(0.001))/2;
           double theta23 = M_PI/4;
           double deltacp = M_PI/2;
           double sdm = 7e-5;
           double ldm = 2e-3; */

        /* Zero oscillation parameters */
        double theta12 = 0.;
        double theta13 = 0.;
        double theta23 = 0.;
        double deltacp = 0.;
        double sdm = 0.;
        double ldm = 0;


        /* Initialize parameter vector(s) */
        glb_params true_values = glbAllocParams();
        glb_params test_values = glbAllocParams();

        glbDefineParams(true_values,theta12,theta13,theta23,deltacp,sdm,ldm);
        glbSetDensityParams(true_values,1.0,GLB_ALL);
        glbDefineParams(test_values,theta12,theta13,theta23,deltacp,sdm,ldm);
        glbSetDensityParams(test_values,1.0,GLB_ALL);

        /* The simulated data are computed */
        glbSetOscillationParameters(true_values);
        glbSetRates();

        int ifile;
        // ## make outfiles for custom binning run
        if(argc ==5) {


                for (ifile=0; ifile<num_chans; ifile++) {

                        if(atoi(argv[4]) == ifile) {

                                sprintf(outfile,"out/%s_%s_%s_events_unweighted.dat",flux_file_name,chan_name[ifile],expt_config_name);
                                printf("%i %s\n",ifile, outfile);


                                FILE* f_out = InitDiffSpecOutput(outfile);

                                ret = glbShowChannelRates(f_out,0,chan_num[ifile],GLB_PRE,GLB_WO_EFF,GLB_WO_BG); // Failing here ?
                                // printf("%d\n", ret );

                                fclose(f_out);

                                sprintf(outfile_smeared,"out/%s_%s_%s_events_smeared_unweighted.dat",flux_file_name,chan_name[ifile],expt_config_name);
                                printf("%i %s\n",ifile,outfile_smeared);


                                FILE* f_out_smeared=InitDiffSpecOutput(outfile_smeared);

                                ret = glbShowChannelRates(f_out_smeared,0,chan_num[ifile],GLB_POST,GLB_W_EFF,GLB_W_BG);


                                fclose(f_out_smeared);
                                break;
                        }
                        else{continue;}

                }
        }   // ## end of make outfiles for custom binning run

        // ##  make outfiles for static binning run
        if(argc==4) {
                for (ifile=0; ifile<num_chans; ifile++) {

                        sprintf(outfile,"out/%s_%s_%s_events_unweighted.dat",flux_file_name,chan_name[ifile],expt_config_name);
                        printf("%i %s\n",ifile, outfile);


                        FILE* f_out = InitDiffSpecOutput(outfile);

                        ret = glbShowChannelRates(f_out,0,chan_num[ifile],GLB_PRE,GLB_WO_EFF,GLB_WO_BG);

                        fclose(f_out);

                        sprintf(outfile_smeared,"out/%s_%s_%s_events_smeared_unweighted.dat",flux_file_name,chan_name[ifile],expt_config_name);
                        printf("%i %s\n",ifile,outfile_smeared);


                        FILE* f_out_smeared=InitDiffSpecOutput(outfile_smeared);

                        ret = glbShowChannelRates(f_out_smeared,0,chan_num[ifile],GLB_POST,GLB_W_EFF,GLB_W_BG);

                        fclose(f_out_smeared);


                }


        }

        /* Now, if necessary, handle the fake background channel which should be tacked onto the end, if present */


        /* Decide to handle background according to whether background file exists */
        sprintf(bgfile,"backgrounds/bg_chan_%s.dat",expt_config_name);

        /*    if (stat(bgfile,&buf) == 0) {
           printf("Background file exists: %s\n",bgfile);
           } else {
           printf("Background file does not exist: %s\n",bgfile);

           } */

        if (stat(bgfile,&buf) == 0)

        {

                sprintf(outfile,"out/%s_bg_chan_%s_events_unweighted.dat",flux_file_name,expt_config_name);
                printf("%i %s\n",ifile, outfile);


                FILE* f_out = InitDiffSpecOutput(outfile);

                ret = glbShowChannelRates(f_out,0,num_chans,GLB_PRE,GLB_WO_EFF,GLB_W_BG);

                fclose(f_out);

                sprintf(outfile_smeared,"out/%s_bg_chan_%s_events_smeared_unweighted.dat",flux_file_name,expt_config_name);
                printf("%i %s\n",ifile,outfile_smeared);



                FILE* f_out_smeared=InitDiffSpecOutput(outfile_smeared);

                ret = glbShowChannelRates(f_out_smeared,0,num_chans,GLB_POST,GLB_W_EFF,GLB_W_BG);

                fclose(f_out_smeared);
        } else {
                printf("No background file\n");
        }



        /* Destroy parameter vector(s) */
        glbFreeParams(true_values);
        glbFreeParams(test_values);

        exit(0);
}
