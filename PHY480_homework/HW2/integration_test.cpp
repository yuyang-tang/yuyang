//  file: integration_test.cpp               
//                                                                     
//  Programmer:  Yuyang Tang    yuyangta@msu.edu 
// 
//  Revision history: 
//      02/26/2020 original version
//      02/28/2020 add some notes and finish homework 4 and 5
//
//
//homework b:
//Simpson's rule:using f(x)=10**(a*log10(x)+b) to fit the data,then find a=-3.78138,b=1.86655.It looks like a line in logscale.
//Milne rule:using f(x)=10**(a*log10(x)+b) to fit the data,then find a=-5.15001,b=3.38446.It looks like a line in logscale.
//
//homework c:
//In the plot we see that by the increase of x(max_intervals),the error will be smaller.In notes,it says that e_total=N**(-beta)+sqrt(N)e_m,so when N get larger,e_total will get smaller.They actually get the same conclusion.so we can get the max N using for Milne rule(as it possible),it will get the best answer.
//
//homework d:
//I used the three methods to calculate the integration of 1./((1.+x)*sqrt(x)) from 0 to 2.It turns out that with Simpson's rule and Milne rule the result is inf,but with GSL function it gets the same answer with the handout(1.9106).I think it is because that using Simpson's rule or Milne rule,it will calculate sum of (f(x)*h),while h is a limited value.But when x get very close to 0,f(x) will be approach to inf,so f(x)*h will be approach to inf,it can't get the correct answer.However,using the function of GSL will not get into this troble.
//
//homework e:
//For Simpson's rule,I add some codes to calculate the integration of 2N,then get the plot of (A(N)-A(2N))/A(2N).From the note,I set f(x)=a*x**b to fit the plot,and get a=0.138764,b=-0.956696.I get the result that beta=-0.956656,so the approxiation error should be alpha/N^(0.956696) (We can not know the value of alpha because we don't know the A_exact.)
//
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
  const double a = 10;//2;              //upper limit of integration
  const double b = 1;//0;               //lower limit of integration
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
  integration_out << "#N   Simpsons           Milne          GSL        error_1         error_2         (A(N)-A(2N))/A(2N)"<<endl;
  integration_out << "-------------------------------------"<<endl;
  //Milne rule should have at least 5 parts and will get more parts by adding 4i
  for (int i= 5; i<=max_intervals; i +=4)
  {
    integration_out << setw(4) << i;
    result = Simpsons_rule (i, a, b, &my_integrand);
    double result_2 = Simpsons_rule (2.*i, a, b, &my_integrand);//get the value when N->2N
    double m=result;                   //get the value of Simpsons integration
    integration_out << " " << setprecision(12) << result;
    result = Milne_rule (i, a, b, &my_integrand);
    double n=result;                   //get the value of Milne integration
    integration_out << " " << setprecision(12) << result;
    gsl_integration_qags (&My_function, b, a, abs_error, rel_error, 1000, work_ptr, &result, &error);
    integration_out << " " << setprecision(12) << result;
    integration_out << " " << scientific << result-m << " " << result-n ;//<< endl;
    integration_out << " " << scientific <<  " " << (m - result_2)/result_2 << endl;
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
  return  (log(2.*x)/sqrt(x));//(1./((1.+x)*sqrt(x)));
}
double My_integrand (double x, void *params)
{
  double alpha = *(double * ) params;
  return  (alpha*log(2.*x)/sqrt(x));//(1./((1.+x)*sqrt(x)*alpha));   
}
