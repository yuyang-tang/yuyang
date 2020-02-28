//  file: integration_test.cpp               
//                                                                     
//  Programmer:  Yuyang Tang    yuyangta@msu.edu 
// 
//  Revision history: 
//      02/26/2020 original version
//
//homework 2:
//Simpson's rule:using f(x)=10**(a*log10(x)+b) to fit the data,then find a=-3.78138,b=1.86655.It looks like a line in logscale.
//Milne rule:using f(x)=10**(a*log10(x)+b) to fit the data,then find a=-5.15001,b=3.38446.It looks like a line in logscale.
//homework 3:
//In the plot we see that by the increase of x(max_intervals),the error will be smaller.In notes,it says that e_total=N**(-beta)+sqrt(N)e_m,so when N get larger,e_total will get smaller.They actually get the same conclusion.so we can get the max N using for Milne rule(as it possible),it will get the best answer.
//************************************************************************
// include files 
#include <iostream> 
#include <iomanip> 
#include <fstream> 
#include <cmath> 
using namespace std;
#include "integration_routines.h"// prototypes for integration routines
#include <gsl/gsl_integration.h>
double my_integrand (double x);
double My_integrand (double x, void *params);
int main()
{
  gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc(1000);
  //set up the integration variables
  const int max_intervals = 501;
  const double a = 10;              //upper limit of integration
  const double b = 1;               //lower limit of integration
  double abs_error= 1.0e-8;         //to avoid round-off problems
  double rel_error= 1.0e-8;         //the result will usually be much better
  double result = 0;                //set up initial value of result
  double error;                     //the estimated error from the integration
  double alpha = 1.0;               //parameter in integrand
  gsl_function My_function;
  void *params_ptr=&alpha;
  My_function.function = &My_integrand;
  My_function.params = params_ptr;
  //open the output file stream
  ofstream integration_out("integration.dat");
  integration_out << "#N   Simpsons           Milne          GSL        error_1         error_2"<<endl;
  integration_out << "-------------------------------------"<<endl;
  //Milne rule should have at least 5 parts and will get more parts by adding 4i
  for (int i= 5; i<=max_intervals; i +=4)
  {
    integration_out << setw(4) << i;
    result = Simpsons_rule (i, a, b, &my_integrand);
    double m=result;                   //get the value of Simpsons integration
    integration_out << " " << setprecision(12) << result;
    result = Milne_rule (i, a, b, &my_integrand);
    double n=result;                   //get the value of Milne integration
    integration_out << " " << setprecision(12) << result;
    gsl_integration_qags (&My_function, b, a, abs_error, rel_error, 1000, work_ptr, &result, &error);
    integration_out << " " << setprecision(12) << result;
    integration_out << " " << scientific << result-m << " " << result-n << endl;
    //get the error of two ways
  }
  cout << "data stored in integration.dat" << endl;
  integration_out.close();           //close the file
  return 0;
}
//****************************************************************************

//the function to be integrated
double my_integrand (double x)
{
  return  (log(2.*x)/sqrt(x));//(log (sqrt(6.*pow(x,3)))-.3*pow(x,4)+5.*exp (x));//+sin(.2*pow(x,4)+0.1)));
}
double My_integrand (double x, void *params)
{
  double alpha = *(double * ) params;
  return  (log(2.*alpha*x)/sqrt(x));//(log (alpha *sqrt(6.*pow(x,3)))-.3*pow(x,4)+5.*exp (x));//+sin(.2*pow(x,4)+0.1)));
}
