//file: round_off_errors.cpp
//
//this is a test program for the round-off errors in single-precision
//
//Programmer: Yuyang Tang yuyangta@msu.edu
//
//Revision history:
//          07-Feb-2020    original version
//          10-Feb-2020    update codes
//
//                    didn't work out
//***************************************************************************
//including files
#include <iostream> 
#include <iomanip> 
#include <fstream> 
#include <cmath> 
using namespace std;
int main()
{
  float z=1.,sum=0;            //set a original value of single-precision numbers
  ofstream errors_out ("errors.dat");    //open the output file stream  
  for (z;z<=1000;z++)
  {
    sum=sum+z;    //calculate the sum of z
  }
  float sqrt_0=sqrt(sum);      //calculate a value of the sqrt as original value
  errors_out<<"sqrt_0="<<setw(16)<<sqrt_0<<endl;
  for(int i=1;i<=1000;i++)
  {
    float sqrt_x=sqrt(sum);
    errors_out<<i<<"    "<<setw(16)<<abs(sqrt_x-sqrt_0)/(sqrt_x+sqrt_0)<<endl;  //output the relative difference in the file
  }
  errors_out.close();
  cout<<"data stored in 'errors.dat'"<<endl;
  return 0;
}
