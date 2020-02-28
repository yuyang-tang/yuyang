//  file: integration_routines.cpp 
// 
//  Program to do integration of an integration of an integrand function 
// 
//  Programmer:  Yuyang Tang  yuyangta@msu.edu 
// 
//  Revision history: 
//      02/22/20  original version  
//      02/26/20  adding notes
// 
//***************************************************************** 
// include files 
#include <cmath>
#include <gsl/gsl_math.h>
//#include <integration_routines.h>   //integration routine prototypes

//using Simpson's rule to do the integration
double Simpsons_rule(int num_pts,double a,double b,double (*integrand)(double x)) //a is the x_max,b is the x_min               
{
   double interval = ((a-b)/double(num_pts - 1.));  // called h in notes
   double sum=  0.;                       // initialize integration sum to zero 
   for (int n=2; n<num_pts; n+=2)                // loop for odd points     
   {     
      double x = b + interval * double(n-1.);     
      sum += (4./3.)*interval * integrand(x);   
   }   
   for (int n=3; n<num_pts; n+=2)                // loop for even points     
   {     
      double x = b + interval * double(n-1.);     
      sum += (2./3.)*interval * integrand(x);   
   }      // add in the end and start point contributions      
   sum +=  (interval/3.) * (integrand(b) + integrand(a));      
   return (sum);
}
//using Milne rule to do the integration
double Milne_rule(int num_pts,double a,double b,double (*integrand)(double x))                  //using Milne integration rule to calculate
{
   double interval = ((a-b)/double(num_pts - 1.));  // called h in notes
   double sum=  0.;                       // initialize integration sum to zero
   for (int n=2; n<num_pts;n+=4)                //loop for 4i+2 points
   {
     double x = b + interval * double(n-1.);
     sum+=(64./45.)*interval*integrand(x); 
   }
   for (int n=3; n<num_pts;n+=4)                //loop for 4i+3 points
   {
     double x = b + interval * double(n-1.);
     sum+=(24./45.)*interval*integrand(x); 
   }
   for (int n=4; n<num_pts;n+=4)                //loop for 4i+4 points
   {
     double x = b + interval * double(n-1.);
     sum+=(64./45.)*interval*integrand(x); 
   }
   for (int n=5; n<num_pts;n+=4)                //loop for 4i+5 points
   {
     double x = b + interval * double(n-1.);
     sum+=(28./45.)*interval*integrand(x); 
   }   //add in the end and start points contributions
   sum += (14./45.) *interval * (integrand(b) + integrand(a)); 
   return sum;
}


