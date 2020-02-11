//  file: summing_up_down.cpp 
//                
//  Programmer:  Yuyang Tang yuyangta@msu.edu
// 
//  Revision history: 
//  04-Feb-2020 original version
//  07-Feb-2020 finish the unfinished code
//*******************************************************************
//homework 2.(c):It has different regions.It seems that there is always a better calculationIn a exact region,and in logscale,it will decrease like linear function,so the power is 10.
//include files
#include <iostream> 
#include <iomanip> 
#include <fstream> 
#include <cmath> 
using namespace std;
int 
main()
{
  ofstream  fplot ("summing_up_down.dat");   //open the plot file stream
  fplot<<"N  "<<"difference"<<endl;
  float sum_up=0 ,sum_down=0;    //set variables equal to 0
  int N;
  float difference;
  cout<< "input the N you want to sum up" <<endl;
  cin>>N;                      //get N
  for (int m=1;m<=N;m++)       //for different N summing up
  {
    for (int n=1;n<=N;n++)     //for different 1/n summing up
    { 
      sum_up=sum_up+1./n;          //sum up 1/n by up order
      sum_down=sum_down+1./(N-n+1);    //sum up 1/n by down order
    }
    difference =2.*abs(sum_up-sum_down)/(sum_up+sum_down); //calculate the difference
    fplot<<m<<"  "<<difference<<endl;      //output the data into file
  }
  fplot.close();
  cout<<"data stored in summing_up_down.dat"<<endl;
  return 0; 
}
