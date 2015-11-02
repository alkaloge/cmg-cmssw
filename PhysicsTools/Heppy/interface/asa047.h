#ifndef __ASA047__
#define __ASA047__

namespace Topness{ class gfuncMET;};
void nelmin ( Topness::gfuncMET& fn, int n, double start[], double xmin[],
              double *ynewlo, double reqmin, double step[], int konvge, int kcount,
              int *icount, int *numres, int *ifault );
/*
void nelmin ( double fn ( double x[] ), //< function to minimize
	      int n,                    //< number of variables
              double start[],           //< start point
              double xmin[],            //> minimum point
              double *ynewlo,           //> function value at minimum
              double reqmin,            //< terminating limit for the variance of function values
              double step[],            //< size and shape of initial simplex
              int konvge,               //< the convergence check is carried out every KONVGE iterations
              int kcount,               //< the maximum number of function evaluations.
              int *icount,              //> the number of function evaluations used
              int *numres,              //> the number of restarts
              int *ifault               //> error indicator: 0, no errors;
                                        //                   1, REQMIN, N, or KONVGE has an illegal 
                                        //                   2, iteration terminated because KCOUNT was exceeded
            );
example:
double rosenbrock ( double x[2] );
int main(){
 ....
  int i;
  int icount;
  int ifault;
  int kcount;
  int konvge;
  int n;
  int numres;
  double reqmin;
  double *start;
  double *step;
  double *xmin;
  double ynewlo;
  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Apply NELMIN to ROSENBROCK function.\n";

  start[0] = -1.2;
  start[1] =  1.0;

  reqmin = 1.0E-08;

  step[0] = 1.0;
  step[1] = 1.0;

  konvge = 10;
  kcount = 500;

  ynewlo = rosenbrock ( start );

  nelmin ( rosenbrock, n, start, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault );

  delete [] start;
  delete [] step;
  delete [] xmin;

}
*/
#endif
