//function to take in eccentricity and mean anomaly and solve the Kepler equation using Danby's (1998)
//iterative method, returning the eccentric anomaly.
//If there is not convergence after 10 iterations then the returned eccentric anomaly is set to 100
//as an error indicator.

#define PI 3.14159265358979323846264338327950288419716939937510

double kepler(double ecc, double man, double acp=1.0e-14)
{
  //ecc=eccentricity, man=mean anomaly, acp=accuracy parameter (defaults to 1e-14)
  int sign, nc;
  //sign=stores sign of a number, nc=counts number of iterations
  float k=0.85;    //factor defining initial value

  double E0, sinm, sinE, esinE, ecosE, E, F, Fp, Fpp, Fppp, g;
  //E0=initial estimate of eccentric anomaly, sinm=sin(mean anomaly)
  //E=eccentric anomaly

  man=fmod(man,(2*PI));       //ensure that mean anomaly is in the range 0 to 2pi
  sinm=sin(man);

  if (sinm>=0.0)
    {
      sign=1;
    }
  else
    {
      sign=-1;
    }

  E0=man+sign*k*ecc;       //set initial estimate of eccentric anomaly
  E=E0;                    //set eccentric anomaly to initial guess
  F=1.0;
  nc=0;

  //iteration loop, stop iteration when solution is found to within accuracy parameter
  while(fabs(F)>acp)
    {
      esinE=ecc*sin(E);
      ecosE=ecc*cos(E);
      F=E-esinE-man;
      Fp=1-ecosE;
      Fpp=esinE;
      Fppp=ecosE;
      g=F/(Fp-0.5*(F*Fpp)/Fp);
      E=E-F/(Fp-0.5*g*Fpp+g*g*Fppp/6);
      nc++;
      if(nc>=10)
	{
	  E=100.0;
	  break;
	}
    }
  return E;
}

