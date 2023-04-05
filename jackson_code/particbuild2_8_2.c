//program to generate new files of particles at a single timestep after integration using mercury
//In version 2.4:   introduced switch between faster (but breaks for >~1000 timesteps) and slower (but
//                  more robust) methods automatically depending on timesteps > or < 1000.
//In version 2.5:   corrected problem of number of outputs not being correctly determined if the particle
//                  TS000001 does not survive until end of simulation.  Tidied program slightly using functions.
//In version 2.6:   changed variables that hold numbers of particles to longs to deal with large simulations.
//In version 2.6.1: changed output format to 'g' type for all floating point numbers for higher precision.
//In version 2.6.2: changed output to include argument of pericentre, longitude of ascending node and mean 
//                  and true anomalies
//In version 2.6.3: added check for infinite semi-major axis (NB: uses C99 INFINITY macro for output).
//                  Also changed back to 'f' output format with a fixed width output equal to standard Mercury input.
//In version 2.7:   made major changes to allow for arbitrary number of output fields rendering many previous
//                  tweaks obsolete.
//In version 2.7.1: Removed a few more bits rendered unnecessary by version 2.7 changes.
//In version 2.8:   Changed old 'fast' and 'slow' methods such that now re-packing is done in batches of 1000 timesteps.
//In version 2.8.1: Changed output timestep file format to 6 digits.
//In version 2.8.2: Adjusted batch reading to ensure information from particles that have been lost is not retained.
//Written by Alan Jackson

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>

//function to determine number of particles
long pcount()
{
    long i, npart;
    char partfile[11], partfileno[6];
    FILE *fin;
    
    i=1;
    while(1==1){
        strcpy(partfile, "TS");
        sprintf(partfileno, "%06ld", i);
        strcat(partfile, partfileno);
        strcat(partfile, ".aei");
        
        fin=fopen(partfile, "r");
        if(fin==NULL){
            break;
        }
        fclose(fin);
        i++;
    }
    npart=i-1;
    return(npart);
}

//function to determine maximum number of outputs in files
//requires end time of simulation and number of particles
long outcount(double tmax, long npart){
    double time;
    char partfile[11], partfileno[6], line[500];
    long nout, k, i;
    FILE *fin;
    
    nout=0;
    i=1;
    while(i<=npart && time<(tmax-1.0)){
        strcpy(partfile, "TS");
        sprintf(partfileno, "%06ld", i);
        strcat(partfile, partfileno);
        strcat(partfile, ".aei");
        
        fin=fopen(partfile, "r");
        fgets(line, 500, fin); 
        fgets(line, 500, fin); 
        fgets(line, 500, fin); 
        fgets(line, 500, fin); 
        
        k=0;
        while(!feof(fin)){
            fscanf(fin, "%lf", &time);
            fgets(line, 500, fin);
            k++;
        }
        fclose(fin);
        if (k>nout){
            nout=k;
        }
        i++;
    }
    return(nout);
}


int main(int argc, char *argv[]){
    FILE *fin;
    
    char partfile[11], partfileno[6], line[500], head[500], junk[20];
    double tmax, time;
    int nswitch, nbatch, nrem;
    long npart, i, j, k, n, m, nfiles;
    
    if (argc<2){
        printf("End time of simulation is required.  Please try again.\n");
        exit(1);
    }
    
    tmax=atof(argv[1]);  //require end time of simulation as an input
    
    printf("End time of simulation: %f\n",tmax);
    
    nswitch=1000;   //nswitch=number of timesteps (and thus output files) present at which to
    //switch between fast and slow methods
    
    //determine no. of particles
    npart=pcount();
    printf("Simulation contains %ld particles\n", npart);
    
    //determine number of outputs in files
    
    nfiles=outcount(tmax, npart);
    
    int Fr[nfiles];
    
    printf("Simulation has %ld outputs\n", nfiles);
    printf("Time samples will be stored in files partic00001.dat to partic%06ld.dat\n", nfiles);
    
    if(nfiles<nswitch){                  //For less than nswitch timesteps use fast method
        //Fast method opens all output files initially and then 
        FILE *fout[nfiles];                //opens, reads and outputs each particle file in turn
        char outfile[16], outfileno[6];
        
        for(i=0; i<nfiles; i++){           //Create output files
            strcpy(outfile, "partic");
            sprintf(outfileno, "%06ld", i+1);
            strcat(outfile, outfileno);
            strcat(outfile, ".dat");
            fout[i]=fopen(outfile, "w");
            Fr[i]=0;                         //set run checks to 0;
        }
        
        printf("Reading TS000001.aei");
        
        for(j=1; j<=npart; j++){          //primary loop for fast method
            k=0;
            
            strcpy(partfile, "TS");
            sprintf(partfileno, "%06ld", j);
            strcat(partfile, partfileno);
            strcat(partfile, ".aei");
            
            fin=fopen(partfile, "r");
            
            printf("\b\b\b\b\b\b\b\b\b\b\b\b%s", partfile);
            
            fgets(line, 500, fin);
            fgets(line, 500, fin); 
            fgets(line, 500, fin);
            
            fscanf(fin, "%s", junk);
            fscanf(fin, "%s", junk);
            fgets(head, 500, fin); 
            
            while(!feof(fin)){                   //reading through individual particles
                fscanf(fin, "%lf", &time);
                fgets(line, 500, fin);
                
                //outputting
                if(Fr[k]==0){     //For first run through print headers to output files
                    fprintf(fout[k],"Time= %f\n", time);
                    fprintf(fout[k], "id   %s", head);
                    for(i=0; i<strlen(head); i++) fputc('-', fout[k]);
                    fputc('\n', fout[k]);
                    Fr[k]=1;
                }
                fprintf(fout[k],"%06ld %s", j, line);
                k++;
            }
            fclose(fin);
        } 
        for(i=0; i<nfiles; i++){
            fclose(fout[i]);    //close output files
        }
    }
    else{                       //For more than nswitch output files, do batches of nswitch output files at a time
        
        int nbatch;
        nbatch = nfiles/nswitch;
        nrem = nfiles % nswitch;
        if(nrem == 0){
            printf("Re-packaging will be run in %d batches\n", nbatch);
        }
        else{
            printf("Re-packaging will be run in %d batches\n", nbatch+1);
        }
        
        char outfile[nfiles][17], outfileno[6];
        
        for(i=0; i<nfiles; i++){                 //Create names of output files
            strcpy(outfile[i], "partic");
            sprintf(outfileno, "%06ld", i+1);
            strcat(outfile[i], outfileno);
            strcat(outfile[i], ".dat");
            Fr[i]=0;                             //set run checks to zero
        }
        
        printf("Batch01 TS000001.aei");
        
        for(i=0; i<nbatch; i++){
            FILE *fout[nswitch];
            for(j=0; j<nswitch; j++){
                if(Fr[j+i*nswitch]==0){    //On first run through create output file
                    fout[j]=fopen(outfile[j+i*nswitch], "w");
                }
                else{        //for subsequent cycles just append to output file
                    fout[j]=fopen(outfile[j+i*nswitch], "a");
                }
            }
            for(j=1; j<=npart; j++){
                k=0;
                
                strcpy(partfile, "TS");
                sprintf(partfileno, "%06ld", j);
                strcat(partfile, partfileno);
                strcat(partfile, ".aei");
                
                fin=fopen(partfile, "r");
                
                printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%02ld %s", i, partfile);
                
                fgets(line, 500, fin);
                fgets(line, 500, fin); 
                fgets(line, 500, fin);
                
                fscanf(fin, "%s", junk);
                fscanf(fin, "%s", junk);
                fgets(head, 500, fin);
                
                if(i>0){
                    for(n=0; n<(i*nswitch); n++){    //Read through particles and dump lines prior to current batch
                        if(!feof(fin)){
                            fscanf(fin, "%lf", &time);
                            fgets(line, 500, fin);
                        }
                    }
                }
                for(n=0; n<nswitch; n++){        //Read through particles and output lines in current batch
                    if(!feof(fin)){
                        fscanf(fin, "%lf", &time);
                        fgets(line, 500, fin);
                    
                        if(Fr[k+i*nswitch]==0){   //On first run through only, print header
                            fprintf(fout[k],"Time= %lf\n", time);
                            fprintf(fout[k], "id        %s", head);
                            for(m=0; m<strlen(head); m++) fputc('-', fout[k]);
                            fputc('\n', fout[k]);
                            Fr[k+i*nswitch]=1;
                        }
                    
                        fprintf(fout[k],"%06ld %s", j, line);
                    }
                    k++;
                }
                fclose(fin);
            }
            for(j=0; j<nswitch; j++){
                fclose(fout[j]);
            }
        }
        
        //Final remainder batch
        
        if(nrem > 0){
            FILE *fpout[nrem];
            for(j=0; j<nrem; j++){
                if(Fr[j+nbatch*nswitch]==0){    //On first run through create output file
                    fpout[j]=fopen(outfile[j+nbatch*nswitch], "w");
                }
                else{        //for subsequent cycles just append to output file
                    fpout[j]=fopen(outfile[j+nbatch*nswitch], "a");
                }
            }
            for(j=1; j<=npart; j++){
                k=0;
                              
                strcpy(partfile, "TS");
                sprintf(partfileno, "%06ld", j);
                strcat(partfile, partfileno);
                strcat(partfile, ".aei");
                
                fin=fopen(partfile, "r");
                
                printf("\b\b\b\b\b\b\b\b\b\b\b\b%s", partfile);
                
                fgets(line, 500, fin);
                fgets(line, 500, fin); 
                fgets(line, 500, fin);
                
                fscanf(fin, "%s", junk);
                fscanf(fin, "%s", junk);
                fgets(head, 500, fin);
                
                for(n=0; n<(nbatch*nswitch); n++){    //Read through particles and dump lines prior to current batch
                    if(!feof(fin)){
                        fscanf(fin, "%lf", &time);
                        fgets(line, 500, fin);
                    }
                }
                
                for(n=0; n<nrem; n++){        //Read through particles and output lines in current batch
                    if(!feof(fin)){
                        fscanf(fin, "%lf", &time);
                        fgets(line, 500, fin);                        
                    
                        if(Fr[k+nbatch*nswitch]==0){   //On first run through only, print header
                            fprintf(fpout[k],"Time= %f\n", time);
                            fprintf(fpout[k], "id        %s", head);
                            for(m=0; m<strlen(head); m++) fputc('-', fpout[k]);
                            fputc('\n', fpout[k]);
                            Fr[k+i*nswitch]=1;
                        }
                   
                        fprintf(fpout[k],"%06ld %s", j, line);
                    }
                    k++;
                }
                fclose(fin);
            }
            for(j=0; j<nrem; j++){
                fclose(fpout[j]);
            }
        }
    }
    
    printf("\nSampling complete\n");
    return(0);
}
