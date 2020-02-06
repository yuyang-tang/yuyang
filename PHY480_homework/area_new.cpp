//  file: area_new.cpp
//
//  This program calculates the area of a circle, given the radius.
//
//  Programmer: Yuyang Tang   yuyangta@msu.edu
//
//  Revision history:
//      04-Feb-2020 original version
//
//  Notes:  
//   * compile with:  "g++ -o area_new.x area_new.cpp"
//**********************************************************************
//include files
#include<iostream>       // this has the cout, cin definitions
#include<cmath>          // this has many  math functions
#include<fstream>        // this can output to a file
using namespace std;     // if omitted, then need std::cout, std::cin
const double pi=3.1415926535897932385 ;//predefined value of pi
double Area(double radius);           //set a fuction to caculate the area
int
main ()
{
  ofstream areaout ("area_new.dat");  //output the value of area into a file
  double radius;    // every variable is declared as int or double or ...

  cout << "Enter the radius of a circle: " << endl;	// ask for radius
  cin >> radius;
  if(radius<0)    //add checks of the input
  {
   areaout<<"radius should be a positive value"<<endl;
   cout<<"radius should be a positive value"<<endl;
  }
  else
  {
   double area = Area(radius);  	// caculate area by a function
   areaout << "radius = " << scientific << radius << ",  area = " <<  scientific    << area << endl;
   cout << "radius = " << scientific << radius << ",  area = " <<  scientific <<    area << endl;
  }

  return 0;			// "0" for successful completion
}
double Area(double radius)
{
  double area=pi*radius*radius;
  return area;
}
