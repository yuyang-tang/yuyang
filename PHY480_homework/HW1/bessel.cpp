//  file: bessel.cpp 
//  added comment//
//  Spherical Bessel functions via up and down recursion and relative difference between the results for each x     
//                                                                     
//
//  Programmer:  Yuyang Tang yuyangta@msu.edu
//
//  Revision history:
//      07-Feb-2020 original version
//
//  homework3 (b):In the plot without logscale,when x is close to 0,the relative difference is about 1,which means that ans_up and ans_down are opposite.When x is from 0 to 40,the relative difference is close to 0,which means that ans_up and ans_down are the same.when x is from 40 to 100 the funtion is like "sin" function,which means that ans_up changes from -ans_down to ans_down.
//************************************************************************

// include files
#include <iostream>		
#include <iomanip>		
#include <fstream>		
#include <cmath> 

using namespace std;		 
#include <gsl/gsl_sf_bessel.h>  //use the GSL library

double down_recursion (double x, int n, int m);	// downward algorithm 
double up_recursion (double x, int n);	        // upward algorithm 

// global constants  
const double xmax = 100.0;	// max of x  
const double xmin = 0.1;	// min of x >0  
const double step = 0.1;	// delta x  
const int order = 10;		// order of Bessel function 
const int start = 50;		// used for downward algorithm 

//********************************************************************
int
main ()
{
  double ans_down, ans_up, j1;

  // open an output file stream
  ofstream my_out ("bessel.dat");

  my_out << "# Spherical Bessel functions via up and down recursion" 
         << endl;

  // step through different x values
  for (double x = xmin; x <= xmax; x += step)
    {
      ans_down = down_recursion (x, order, start);
      ans_up = up_recursion (x, order);
      j1=gsl_sf_bessel_j1(x);      //use GSL function to calculate j1
      my_out << fixed << setprecision (6) << setw (8) << x << " "
	<< scientific << setprecision (6)
	<< setw (13) << ans_down << " "
	<< setw (13) << ans_up << " "
	<< setw (13) << abs(ans_up -ans_down)/(abs(ans_up)+abs(ans_down))<< " "
	<< setw (13) << j1 
        << endl;                   //output the data in a file
    }
  cout << "data stored in bessel.dat." << endl;

  // close the output file
  my_out.close ();
  return (0);
}


//------------------------end of main program----------------------- 

// function using downward recursion  
double
down_recursion (double x, int n, int m)
{
  double j[start + 2];		// array to store Bessel functions 
  j[m + 1] = j[m] = 1.;		// start with "something" (choose 1 here) 
  for (int k = m; k > 0; k--)
    {
      j[k - 1] = ((2.* double(k) + 1.) / x) * j[k] - j[k + 1];  // recur. rel.
    }
  double scale = (sin (x) / x) / j[0];	// scale the result 
  return (j[n] * scale);
}


//------------------------------------------------------------------ 

// function using upward recursion  
double
up_recursion (double x, int n)
{
  double term_three = 0.;
  double term_one = (sin (x)) / x;	// start with lowest order 
  double term_two = (sin (x) - x * cos (x)) / (x * x);	// next order
  for (int k = 1; k < n; k += 1)	// loop for order of function     
    { // recurrence relation
      term_three = ((2.*double(k) + 1.) / x) * term_two - term_one;	       
      term_one = term_two;
      term_two = term_three;
    }
  return (term_three);
}
