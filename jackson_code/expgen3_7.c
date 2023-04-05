//Generates a spherically symmetric explosion of particles with a truncated gaussian, or uniform, 
//kick velocity distribution and generates orbital elements of particles suitable for plotting in
//IDL or other and produces a set of particles for input into a Mercury simulation.
//It is completely dimensionless with the kick velocity measured in terms of the orbital velocity
//and the (new) semi-major axis measured in terms of the old semi-major axis.
//Set up to work with hard parallelisation via init_sim bash script, can also be used independently however.
//
//created by Alan Jackson
//
//In version 3.2 switched from variables describing velocity distribution being hard-coded to being drawn
//from an input file.
//In version 3.3 added calculation of true anomaly and argument of pericentre and correction to dv
//for loss of KE to escape from parent (only needed for orbits).
//In version 3.4 removed assumption of initial eccentricity=0.  Reference plane is still plane of orbit,
//particle cartesian positions are still centred on x-axis, though orbits are given with x-axis in direction
//of old orbit pericentre.
//In version 3.4.1 made minor adjustments to inputs to allow zero planet mass and zero start radius (only when
//planet mass is zero)
//In version 3.4.2 changed output formats to 'g' type to provide higher precision.
//In version 3.4.3 fixed bug introduced in velocity calculation in 3.4.1
//In version 3.5 introduced option for uniform distributions
//In version 3.6 introduced option to define distribution after escape.
//In version 3.7 introduced command line option to set mass of debris particles and adjusted command line processing
//note, version 3.7 will only work with init_sim5 or greater.

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<time.h>
#include<string.h>
#include<ctype.h>
#include<unistd.h>

#define AU 1.496E11     //1 AU in m
#define PI 3.14159265
#define ds 86400        //1 day in seconds
#define MSOL 1.98892E30 //1 solar mass in kg
#define G 6.673E-11     //gravitational constant
#define ME 5.9736E24    //1 Earth mass in kg

int main(int argc, char** argv)
{
  int nstart, nproc, np, i, j, k, dcheck=0, bacheck=0;
  //np=no. of particles
  //nproc=processor number, need to make sure all particles running on each processor have unique IDs
  //nstart=position at which to start numbering particles

  double debmass=0.0;
  //debmass=mass of each debris particle (kg)

  char *optval = NULL;

  while ((i = getopt (argc, argv, "p:n:m:")) != -1){
    switch (i)
      {
      case 'p':
	optval=optarg;
	nproc=atoi(optval);
	break;
      case 'n':
	optval=optarg;
	np=atoi(optval);
	break;
      case 'm':
	optval=optarg;
	debmass=atof(optval);
	break;
      case '?':
	if (optopt == 'p' || optopt == 'n' || optopt == 'm'){
	  printf("Option -%c requires an argument\n", optopt);
	  exit(1);
	}
	else{
	  printf("Unrecognised option\n");
	  exit(1);
	}
      }
  }

  nstart=np*(nproc-1);
  debmass=debmass/MSOL;

  double sinthetap, theta, stheta, ctheta, phi, sphi, cphi, dv, dvc, aprim, hprimh2, ein, ecc, inc; 
  //sinthetap=sin(Pi/2-theta) with range -1 to 1, to get correct distribution of theta, 
  //theta=angle to z axis(latitude), phi=angle to x axis in xy plane(longitude)
  //dv=magnitude of kick velocity
  //dvc=magnitude of kick after escape from parent
  //aprim=new semi-major axis, hprimh2=(h'/h)^2
  //ein=(initial) eccentricity, ecc=(new) eccentricity, inc=(new) inclination

  double fin, sfin, cfin, fp, sfp, cfp, argperi, lascn, meanv=0, sdv=0, vcut, sphif, cphif;
  //fin=initial true anomaly
  //fp=(new) true anomaly, argperi=(new) argument of pericentre, laccn=(new)longitude of ascending node
  //meanv=mean kick velocity, sdv=std. dev. of gaussian, vcut=cut-off velocity, all in units of Earth orb. vel.
  //sphif=sin(phi-f), cphif=cos(phi-f)

  double x, y, z, vx, vy, vz, xpl=0, stmass=0, vorb, vorbm, radst=0, radpl=0, pmass=0, rad, radm;
  //x, y, z, vx, vy, vz are variables for cartesian positions and velocities of particles.
  //xpl=semi-major axis of parent planet (AU), position of parent is set as being on the x axis at time of explosion
  //stmass=mass of central star (solar masses), rad=distance of particles from parent at start of integration (km)
  //radst=distance of particles from parent at start of integration (parent radii)
  //radpl=radius of parent (km), vorb=orbital velocity of parent (km/s)
  //pmass=mass of parent (ME)
  //vorbm, radm, conversions of vorb, rad into Mercury units (AU/day, AU)

  char car, junk[200], var1[200], unit[20], var2[20], *varfind, infile[]="exppar3_6.in", dist[20], baesc[20];
  //car=single character storage variable, junk=junk line storage, var1,var2=variable name storage
  //varfind=substring location, unit=variable to store units of vel dist, infile=input file name

  fpos_t posit;
  //file position indicator

  FILE *fout, *fmer, *finput;
  fout=fopen("exp3elem.data", "w");
  fmer=fopen("small.in", "w");
  finput=fopen(infile, "r");

  if(fout==NULL){
    printf("Cannot open exp3elem.data\n");
    exit(1);
  }
  if(fmer==NULL){
    printf("Cannot open small.in\n");
    exit(1);
  }
  if(finput==NULL){
    printf("Cannot open %s\n", infile);
    exit(1);
  }

  strcpy(unit,"rel"); //default units of velocity distribution to relative

  //Read in initialisation parameters
  while(!feof(finput)){
      fgetpos(finput, &posit);
      car=fgetc(finput);
      if (car == ')'){
	  fgets(junk,200, finput);
      }
      else if (feof(finput)){
	break;
      }
      else
	{
	  fsetpos(finput, &posit);
	  j=0;
	  do{
	    car=fgetc(finput);
	    car=tolower(car);
	    var1[j]=car;
	    j++;
	  }while(car!='=');
	  
	  if((varfind = strstr(var1, "type"))!=NULL) fgets(dist,20,finput);
	  else if((varfind = strstr(var1, "units"))!=NULL) fgets(unit,20,finput);
	  else if((varfind = strstr(var1, "defined"))!=NULL) fgets(baesc,20,finput);
	  else if((varfind = strstr(var1, "mean"))!=NULL){
	    fgets(var2,20,finput);
	    meanv=atof(var2);
	  }
	  else if((varfind = strstr(var1, "standard"))!=NULL){
	    fgets(var2,20,finput);
	    sdv=atof(var2);
	  }
	  else if((varfind = strstr(var1, "trunc"))!=NULL){
	    fgets(var2,20,finput);
	    vcut=atof(var2);
	  }
	  else if((varfind = strstr(var1, "start"))!=NULL){
	    fgets(var2,20,finput);
	    radst=atof(var2);
	  }
	  else if((varfind = strstr(var1, "semi"))!=NULL){
	    fgets(var2,20,finput);
	    xpl=atof(var2);
	  }
	  else if((varfind = strstr(var1, "ecc"))!=NULL){
	    fgets(var2,20,finput);
	    ein=atof(var2);
	  }
	  else if((varfind = strstr(var1, "true"))!=NULL){
	    fgets(var2,20,finput);
	    fin=atof(var2)*PI/180.0;
	  }
	  else if((varfind = strstr(var1, "star"))!=NULL){
	    fgets(var2,20,finput);
	    stmass=atof(var2);
	  }
	  else if((varfind = strstr(var1, "parent radius"))!=NULL){
	    fgets(var2,20,finput);
	    radpl=atof(var2);
	  }
	  else if((varfind = strstr(var1, "parent mass"))!=NULL){
	    fgets(var2,20,finput);
	    pmass=atof(var2);
	  }
	  else{
	    printf("Invalid parameter entry detected, please check %s\n", infile);
	    exit(1);
	  }
	  memset(var1, '\0', strlen(var1));
	  memset(var1, '\0', strlen(var2));
	}
  }
  fclose(finput);

  j=0;
  while (j<=strlen(baesc)){
    baesc[j]=tolower(baesc[j]);
    j++;
  }

  j=0;
  while (j>=strlen(dist)){
    dist[j]=tolower(dist[j]);
    j++;
  }

  if ((varfind = strstr(dist, "uni"))!=NULL){
    printf("Using uniform distribution\n");
    dcheck=1;
  }
  else if (sdv==0){
    printf("Standard deviation of velocity distribution input as zero, switching to delta function velocity distribution\n");
  }
  if ((varfind = strstr(baesc, "aft"))!=NULL){
    printf("Velocity distribution defined after escape\n");
    bacheck=1;
  }
  if (radst==0){
    printf("Start radius set to zero\n");
  }
  if (pmass==0){
    printf("Mass of parent set to zero\n");
  }
  else{
    if (radst==0){
      printf("Zero start radius not valid with non-zero planet mass, please check %s\n", infile);
      exit(1);
    }
    if (xpl==0){
      printf("Semi-major axis of parent not set correctly, please check %s\n", infile);
      exit(1);
    }
    if (radpl==0){
      printf("Radius of parent not set correctly, please check %s\n", infile);
      exit(1);
    }
  }

  //If velocity distribution parameters were entered in absolute terms convert to relative (non-dimensional)
  j=0;
  while (j<=strlen(unit)){
    unit[j]=tolower(unit[j]);
    j++;
  }

  stmass=stmass*MSOL;

  vorb=sqrt(G*stmass/(xpl*AU));
  //printf("orbital velocity=%f m/s\n",vorb);

  if ((varfind = strstr(unit, "abs"))!=NULL){
    meanv=meanv*1000/vorb;
    sdv=sdv*1000/vorb;
    vcut=vcut*1000/vorb;
  }

  rad=radst*radpl*1000.0;     //Initial distance from parent at which particles are started, in metres
 
  //convert parent parameters into Mercury units
  vorbm=vorb/AU*ds;
  radm=rad/AU;


  const gsl_rng_type * T;    //initialize random number generator 
  gsl_rng * r;               //T=(uniform) random number generator to use, here default is used
                             //r=pointer to instance of random number generator
  gsl_rng_env_setup();

  T=gsl_rng_mt19937;         //Use Mersenne Twister
  r=gsl_rng_alloc(T);

  while(clock() < CLOCKS_PER_SEC){}

  gsl_rng_set(r, time(NULL));      //seed random number generator with system time, above waiting loop ensures
                                   //system time has advanced by at least 1 second between runs

  //write header to output files
  fprintf(fout, "id      dv      dvc      theta      phi        a       ecc      inc    argperi    lascn   tranom\n");
  fprintf(fout, "------------------------------------------------------------------------------------------------\n");

  fprintf(fmer, ")O+_06 Small-body initial data  (WARNING: Do not delete this line!!)\n");
  fprintf(fmer, ") Lines beginning with `)' are ignored.\n");
  fprintf(fmer, ")---------------------------------------------------------------------\n");
  fprintf(fmer, " style (Cartesian, Asteroidal, Cometary) = Cartesian\n");
  fprintf(fmer, ")---------------------------------------------------------------------\n");


  dv=0.0;

  //generate directions of kick velocities and new orbital elements and output to files
  for(i=1; i<=np; i++){
    if(dcheck == 1){
      dv=gsl_rng_uniform(r)*(meanv-vcut)+vcut;
    }
    else if(sdv>0.0){
      do{
	dv=gsl_ran_gaussian(r, sdv)+meanv;
      }while(dv<vcut);
    }
    else dv=meanv;

    //If velocity distribution defined before escape correct dv for loss of kinetic energy
    //in escaping parent gravity for calculating orbits
    if((pmass > 0.0) && (bacheck == 0)){
      dvc=sqrt(dv*dv-(2*G*pmass*ME/rad)/(vorb*vorb));
    }
    else dvc=dv;

    //If velocity distribution defined after escape correct dv for kinetic energy required
    //for escape for calculating positions and velocities
    if ((bacheck == 1) && (pmass > 0.0)) dv=sqrt(dv*dv+(2*G*pmass*ME/rad)/(vorb*vorb));

    phi=gsl_rng_uniform(r)*2*PI;
    sinthetap=(gsl_rng_uniform(r)*2)-1;

    theta=PI/2-asin(sinthetap);

    stheta=sin(theta);
    ctheta=cos(theta);
    sphi=sin(phi);
    cphi=cos(phi);

    z=radm*ctheta;
    x=radm*stheta*cphi+xpl*(1.0-ein*ein)/(1.0+ein*cos(fin));
    y=radm*stheta*sphi;

    vz=dv*vorbm*ctheta;
    vx=(dv*stheta*cphi-ein*sin(2*PI-fin))*vorbm/sqrt(1.0-ein*ein);
    vy=(dv*stheta*sphi+ 1.0+ein*cos(2*PI-fin))*vorbm/sqrt(1.0-ein*ein);

    if(ein==0.0){ //If ein=0 use simpler (faster) circular form of equations
      //printf("%f\n", ein);

      aprim=1.0/(1.0-dvc*(dvc+2*stheta*sphi));
      
      hprimh2=1+2*dvc*stheta*sphi+dvc*dvc*(ctheta*ctheta+stheta*stheta*sphi*sphi);
      
      ecc=sqrt(1.0-hprimh2/aprim);
      
      inc=acos((1 + dvc*stheta*sphi)/sqrt(hprimh2));
      
      sfp=(1.0/ecc)*sqrt(hprimh2)*stheta*cphi*dvc;
      
      cfp=(1.0/ecc)*(hprimh2-1.0);
      
      fp=atan2(sfp, cfp); //true anomaly should be (0,2*PI), but atan2 gives (-PI, PI), so correct
      if(fp < 0.0){
	  fp=fp+2*PI;
	}
      
      aprim=aprim*xpl;   //Re-dimensionalise a in terms of parent semi-major axis
      
      if(ctheta>=0.0){
	  lascn=0.0;
	  argperi=PI-fp;
	  if(argperi < 0.0){
	    argperi=argperi+2*PI;
	  }
      }
      else{
	lascn=PI;
	argperi=2*PI-fp;
      }
    }
    else{
      sphif=sin(phi-fin);
      cphif=cos(phi-fin);
      sfin=sin(fin);
      cfin=cos(fin);

      aprim=1.0/(1.0 - dvc*(dvc + (2.0/sqrt(1.0-ein*ein))*stheta*(sphif + ein*sphi)));

      hprimh2=1 + 2*sqrt(1.0-ein*ein)*dvc*stheta*sphif/(1.0+ein*cfin) + (1-ein*ein)*dvc*dvc*(ctheta*ctheta + stheta*stheta*sphif*sphif)/pow((1.0 + ein*cfin),2);

      ecc=sqrt(1.0-(1-ein*ein)*hprimh2/aprim);

      inc=acos((1 + sqrt(1-ein*ein)*dvc*stheta*sphif/(1+ein*cfin))/sqrt(hprimh2));

      sfp=(1.0/ecc)*sqrt(hprimh2)*(ein*sfin + sqrt(1-ein*ein)*stheta*cphif*dvc);

      cfp=(1.0/ecc)*(hprimh2*(1+ein*cfin)-1);

      fp=atan2(sfp, cfp); //true anomaly should be (0,2*PI), but atan2 gives (-PI, PI), so correct
      if(fp < 0.0){
	  fp=fp+2*PI;
	}
      
      aprim=aprim*xpl;   //Re-dimensionalise a in terms of parent semi-major axis
      
      if(ctheta>=0.0){
	lascn=fin;
	argperi=PI-fp;
	if(argperi < 0.0) argperi=argperi+2*PI;
      }
      else{
	lascn=fin+PI;
	argperi=2*PI-fp;
      }
    }
    fprintf(fout, "%06d %g %g %g %g %g %g %g %g %g %g\n", i+nstart, dv, dvc, theta, phi, aprim, ecc, inc, argperi, lascn, fp);

    fprintf(fmer, "TS%06d ep=0.0 m=%g d=2.5\n",i+nstart, debmass);
    fprintf(fmer,"%g %g %g %g %g %g 0 0 0\n", x, y, z, vx, vy, vz);
  }

  fclose(fout);
  fclose(fmer);

  return(0);
}
