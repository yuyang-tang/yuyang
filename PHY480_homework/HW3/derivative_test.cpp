//  file: derivative_test.cpp
// 
//  Program to study the error in differentiation rules
//
//  Programmer:  Yuyang Tang   yuyangta@msu.edu
//
//  Original code:session_04  derivative_test.cpp
//   
//
//
//
//*****************************************************************
// include files
#include <iostream>		
#include <iomanip>		
#include <fstream>		
#include <cmath>
using namespace std;		
#include <gsl/gsl_math.h>
#include <gsl/gsl_diff.h>
 

typedef struct
{
double alpha;
double beta;
}
osu_parameters;


double test_function (double x, void *params_ptr);
double test_function_derivative (double x, void *params_ptr);

double forward_diff (double x, double h,
		     double (*f) (double x, void *params_ptr),
		     void *params_ptr);
double central_diff (double x, double h,
		     double (*f) (double x, void *params_ptr),
		     void *params_ptr);
double extrap_diff (double x, double h,
		    double (*f) (double x, void *params_ptr),
		    void *params_ptr);
double extrap_diff2 (double x, double h,
		    double (*f) (double x, void *params_ptr),
		    void *params_ptr);

//************************** main program ***************************
int
main (void)
{
  void *params_ptr;		// void pointer passed to functions 
  osu_parameters alpha_beta;
  alpha_beta.alpha=2.;
  alpha_beta.beta=3./2.;
  const double hmin = 1.e-10;	// minimum mesh size 
  double x = 2.;		// find the derivative at x 
  double diff_cd, diff_fd;	// central, forward difference 
  double diff_extrap;		// extrapolated derivative using central_diff
  double diff_extrap2;          // extrapolated derivative
  double diff_gsl_cd;		// gsl adaptive central derivative 
  gsl_function My_F;		// gsl_function type 
  double abserr;                // absolute error

  ofstream out ("derivative_test.dat");	// open the output file 

  params_ptr = &alpha_beta;		// double to pass to function 

  // exact answer for test 
  double answer = test_function_derivative (x, params_ptr);	

  My_F.function = &test_function;	// set up the gsl function 
  My_F.params = params_ptr;
  gsl_diff_central (&My_F, x, &diff_gsl_cd, &abserr);	// gsl calculation

  cout << "gsl_diff_central(" << x << ") = " << scientific
    << setprecision (16) << diff_gsl_cd << " +/- "
    << setprecision (6) << abserr << endl;
  cout << " actual relative error: " << setprecision (8)
    << fabs((diff_gsl_cd - answer)/answer) << endl;

  double h = 0.5;		// initialize mesh spacing 
  out<<"log10(h)  "<<"diff_fd  "<<"diff_cd  "<<"diff_extrap  "<< "diff_extrap2  "<<endl;
  while (h >= hmin)
  {
    diff_fd = forward_diff (x, h, &test_function, params_ptr);
    diff_cd = central_diff (x, h, &test_function, params_ptr);
    diff_extrap = extrap_diff (x, h, &test_function, params_ptr);
    diff_extrap2 = extrap_diff2 (x, h, &test_function, params_ptr);
    // print relative errors to output file 
    out << scientific << setprecision (8)
      << log10 (h) << "   "
      << log10 (fabs ((diff_fd - answer) / answer)) << "   "
      << log10 (fabs ((diff_cd - answer) / answer)) << "   "
      << log10 (fabs ((diff_extrap - answer) / answer)) <<"   "
      <<log10 (fabs ((diff_extrap2 - answer)/ answer))<< endl;

    h /= 2.;		// reduce mesh by 2 
  }

  out.close ();         // close the output stream
  return (0);		// successful completion 
}

//************************** funct ***************************
double
test_function (double x, void *params_ptr)
{
  osu_parameters *passed_ptr;
  passed_ptr = (osu_parameters *) params_ptr;
  double alpha=((osu_parameters *) params_ptr)->alpha;
  double beta=passed_ptr->beta;
  
  return (alpha * pow(x,beta));
}

//************************** funct_deriv *********************
double
test_function_derivative (double x, void *params_ptr)
{
  osu_parameters *passed_ptr;
  passed_ptr = (osu_parameters *) params_ptr;
  double alpha=((osu_parameters *) params_ptr)->alpha;
  double beta=passed_ptr->beta;

  return (alpha * beta * pow(x,beta-1.));
}


//************************** forward_diff *********************
double
forward_diff (double x, double h,
	      double (*f) (double x, void *params_ptr), void *params_ptr)
{
  return ( f(x + h, params_ptr) - f(x, params_ptr) ) / h;
}

//************************** central_diff *********************
double
central_diff (double x, double h,
	      double (*f) (double x, void *params_ptr), void *params_ptr)
{
  return ( f(x + h/2., params_ptr) - f(x - h/2., params_ptr) ) / h;
}

//************************** extrap_diff *********************
double
extrap_diff (double x, double h,
	     double (*f) (double x, void *params_ptr), void *params_ptr)
{
 // return ( 8.*(f(x + h/4., params_ptr) - f(x - h/4., params_ptr))
 //	  - (f(x + h/2., params_ptr) - f(x - h/2., params_ptr)) ) 
 //	  / (3.*h);
   return (4.*central_diff(x,h/2.,f,params_ptr)-central_diff(x,h,f,params_ptr))/3.;
}
//************************** extrap_diff2 *********************
double
extrap_diff2 (double x, double h,
	     double (*f) (double x, void *params_ptr), void *params_ptr)
{
   return ( 8.*(f(x + h/4., params_ptr) - f(x - h/4., params_ptr))
 	  - (f(x + h/2., params_ptr) - f(x - h/2., params_ptr)) ) 
 	  / (3.*h);
}
