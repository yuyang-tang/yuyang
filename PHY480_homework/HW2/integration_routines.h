//  file: integration_routines.h
// 
//  Header file for integration_routines.cpp
//
//
//  Programmer:  Yuyang Tang        yuyangta@msu.edu
//
//  Revision History:
//    02/26/2020       orignal version
//
//
//************************************************************************

//  begin: function prototypes 
 
extern double Simpsons_rule ( int num_pts, double a, double b, 
                       double (*integrand) (double x) );    // Simpson's rule 
extern double Milne_rule ( int num_pts, double a, double b, 
                       double (*integrand) (double x) );    // Milne rule 

//  end: function prototypes 
