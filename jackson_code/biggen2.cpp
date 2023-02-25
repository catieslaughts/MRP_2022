//program to generate big.in for Mercury simulations of giant impact debris
//assumes all planetary orbits are co-planar
//Originating planet is started on x-axis, others placed at random position on their orbit
//and all orbit orientations randomised.
//Essentially the purpose of this program is to enable orbits of planets to be input
//but to use Cartesian format for Mercury big.in to allow accurate placement of particles

//This program should only be run once for each simulation, not each segment
//of the simulation to ensure that the initial conditions are constant across
//all segments.

#include<cstdio>
#include<cmath>
#include<cstdlib>
#include<gsl/gsl_rng.h>
#include<cstring>
#include<cctype>
#include<iostream>
#include<fstream>
#include"kepler2.h"


//#define RE 4.25796E-5  //radius of the Earth in AU
#define AUm 1.496E11   //1 AU in m
#define ds 86400       //1 day in s
//#define PI 3.14159265
#define G 6.67E-11     //gravitational constant
#define Msol 1.988E30  //Solar mass in kg
#define ME 5.9736E24   //Earth mass in kg

using namespace std;

struct planet //structure to hold planet information
{
  //input
  double a;
  double e;
  double mass;
  double dens;
  string name;
  //calculated
  double tranom;
  double x;
  double y;
  double vx;
  double vy;
};

void plposvel (planet *plan, double stmass, gsl_rng * r, int op=0)
{
  double peri, manom, Eanom, vorb, rad;  
  int count=0;

  peri=gsl_rng_uniform(r)*2*PI;

  if(op==1){
    plan->tranom = 2*PI-peri;
    cout << "True anomaly of originating planet " << plan->tranom*180.0/PI <<" (deg)" << endl;
    rad=plan->a*(1.0-plan->e*plan->e)/(1.0+plan->e*cos(plan->tranom));
    plan->y = 0.0;
    plan->x = rad;
  }
  else{
    do{
      manom=gsl_rng_uniform(r)*2*PI;
      Eanom=kepler(plan->e, manom, 1.0e-14);
      if(Eanom==100)count++;
    }while(count <=10 && Eanom==100);
    if(Eanom==100 && count>10){
      cerr << "Cannot calculate position of planet " << plan->name
           << " please check input file" << endl;
      exit(1);
    }
    plan->tranom = 2.0*atan2((sqrt((1.0+plan->e)/(1.0-plan->e))*tan(Eanom/2.0)),1.0);
    rad=plan->a*(1.0-plan->e*cos(Eanom));
    plan->x = rad*cos(peri+plan->tranom);
    plan->y = rad*sin(peri+plan->tranom);
  }

  vorb=sqrt(G*stmass*Msol/(plan->a*AUm*(1.0-plan->e*plan->e)));
  vorb=vorb/AUm*ds;          //convert to AU/day
  plan->vx=-vorb*(sin(peri+plan->tranom)+plan->e*sin(peri));
  plan->vy=vorb*(cos(peri+plan->tranom)+plan->e*cos(peri));
}

int main(int argc, char *argv[])
{
  int nplan, cpl=0, i, op;
  double stmass;
  string junk, var1, var2;
  char car;
  size_t varfind;

  if(argc!=2){
    cerr << "Total number of planets in system required as command line argument" << endl;
    exit(1);
  }

  nplan=atoi(argv[1]);

  planet planets[nplan];
  for(i=0; i<=nplan; i++) planets[i].mass=0;

  ifstream fin("bigpar.in");
  ofstream fout("big.in");

  if(!fin.is_open()){
    cerr << "Cannot open initialisation file bigpar.in" << endl;
    exit(1);
  }
  if(!fout.is_open()){
    cerr << "Cannot open output file big.in" << endl;
    exit(1);
  }

  while(!fin.eof()){
    car=fin.peek();
    if (car == ')'){
      getline(fin, junk);
    }
    else if (fin.eof()){
      break;
    }
    else{
      do{
	car=fin.get();
	car=tolower(car);
	var1.push_back(car);
      }while(car!='=');

      if((varfind = var1.find("stellar"))!=string::npos){
	getline(fin, var2);
	stmass=atof(var2.c_str());
      }
      else if((varfind = var1.find("total"))!=string::npos){
	getline(fin, var2);
	//nplan=atoi(var2);
	//planet planets[nplan];
      }
      else if(nplan==0){
	cerr << "Number of planets must be declared before planets can be initialised" << endl;
	exit(1);
      }
      else if((varfind = var1.find("name"))!=string::npos){
	getline(fin, var2);
	if(planets[cpl].mass != 0.0) cpl++;
	planets[cpl].name=var2;
      }
      else if((varfind = var1.find("density"))!=string::npos){
	getline(fin, var2);
	planets[cpl].dens=atof(var2.c_str());
      }
      else if((varfind = var1.find("mass"))!=string::npos){
	getline(fin, var2);
	planets[cpl].mass=atof(var2.c_str());
      }
      else if((varfind = var1.find("semi"))!=string::npos){
	getline(fin, var2);
	planets[cpl].a=atof(var2.c_str());
      }
      else if((varfind = var1.find("ecc"))!=string::npos){
	getline(fin, var2);
	planets[cpl].e=atof(var2.c_str());
      }
      else{
	cerr << "Invalid parameter entry detected, please check bigpar.in" << endl;
	exit(1);
      }
      var1.clear();
      var2.clear();
    }
  }
  fin.close();

  const gsl_rng_type * T;    //initialize random number generator 
  gsl_rng * r;               //T=(uniform) random number generator to use, here default is used
                             //r=pointer to instance of random number generator
  gsl_rng_env_setup();

  T=gsl_rng_default;
  r=gsl_rng_alloc(T);

  gsl_rng_set(r, time(NULL));      //seed random number generator with system time

  //output standard Mercury header
  fout << ")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n"
       << ") Lines beginning with `)' are ignored.\n"
       << ")---------------------------------------------------------------------\n"
       << " style (Cartesian, Asteroidal, Cometary) = Cartesian\n"
       << " epoch (in days) = 0.0" << endl;

  //Originating planet
  plposvel(&planets[0], stmass, r, op=1);

  fout << ")True anomaly of originating planet: " << planets[0].name << " = " << planets[0].tranom*180.0/PI << " (deg)\n"
       << ")---------------------------------------------------------------------" << endl;

  fout << planets[0].name << " m=" << planets[0].mass*ME/Msol << " r=3.d0 d=" << planets[0].dens << "\n"
       << planets[0].x << " 0.0 0.0\n"
       << planets[0].vx << " " << planets[0].vy << " 0.0\n"
       << "0.0 0.0 0.0" << endl;

  for(i=1; i<=cpl; i++){
    plposvel(&planets[i], stmass, r, op=0);
    fout << planets[i].name << " m=" << planets[i].mass*ME/Msol << " r=3.d0 d=" << planets[i].dens << "\n"
	 << planets[i].x << " "<< planets[i].y << " 0.0\n"
	 << planets[i].vx << " " << planets[i].vy << " 0.0\n"
	 << "0.0 0.0 0.0" << endl;
  }
  
  fout.close();

  return(0);
}
