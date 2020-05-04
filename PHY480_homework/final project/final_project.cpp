//program: final_project.cpp
//This program is to solve Jacobian iteration
//
//programer: Yuyang Tang    yuyangta@msu.edu
//
//
//One attept is as follow:
//--------Jacobian iteration--------input the number of unknowns x and the //number of equations y: 3 3
//input the augmented matrix:10 -1 -2 7.2
//-1 10 -2 8.3
//-1 -1 5 4.2
//input the accuracy e and the max number of iterations 0.01 100
//10.000000 -1.000000 -2.000000 7.200000
//-1.000000 10.000000 -2.000000 8.300000
//-1.000000 -1.000000 5.000000 4.200000
//result after d interation0.720000 0.830000 0.840000
//result after d interation0.971000 1.070000 1.150000
//result after d interation1.057000 1.157100 1.248200
//result after d interation1.085350 1.185340 1.282820
//result after d interation1.095098 1.195099 1.294138
//*************************************************


#include<iostream>
#include<cmath>
#include<fstream>
 using namespace std;
 
double A[10][10];
double re[10];
void swapA(int i,int j,int n)   //swap the i line and j line
 {

     for(int x = 0;x<=n;x++) 
     {
         double temp = A[i][x];
         A[i][x] = A[j][x];
         A[j][x] = temp;
     }
 }
void getResult(int n,double e,int N) 
 {
     for(int i = 0;i<n;i++)
     {
         re[i] = 0.0;
     } 
     for(int i = 0;i<n;i++)   //to judge that if the diagonal elements are equal to 0.If they are,then swap the two column
     {
        if(fabs(A[i][i]-0)<=1e-9)
        {
              int j;
              for(j = 0;j<n;j++)
	      {
                 if(fabs(A[j][i]-0)>1e-9)
                 {
                    swapA(i,j,n);
                    break;
                 }
              }
             if(j>=n)
             {
                printf("it's invalid matrix!");
             }
             i = 0;//again from 0 
        }
     }
     for(int i = 0;i<n;i++) 
     {
         for(int j = 0;j<n+1;j++)
	 {
             printf("%lf ",A[i][j]);
         } 
         cout<<endl;
     }
     //start interating from the bottom
     int k = 0;
     double x[10];
     for(int i = 0;i<n;i++){
         x[i] = 0.0;
     }
     while(k<=N){
         k++;
         if(k>N) {
             printf("failed!");
             exit(0);
         }
         for(int i = 0;i<n;i++){
             re[i] = A[i][n];
             for(int j = 0;j<n;j++){
                 if(j!=i){
                     re[i] = re[i] - A[i][j]*x[j];        
                 }
             }
             re[i] = re[i] / A[i][i];
         }
	 //when the max x is smaller than unit matrix I.
         double maxXerror = 0.0; 
         for(int i = 0;i<n;i++){
             if(fabs(x[i]-re[i]) >maxXerror){
                 maxXerror = fabs(x[i] - re[i]);
             }
         }
         if(maxXerror < e){
             return;
         }
         printf("result after d interation",k); 
         for(int i = 0;i<n;i++) {
             printf("%lf ",re[i]);
         }
         cout<<endl;
	 //if not,then go on
         for(int i = 0;i<n;i++){
             x[i] = re[i];
         }
     }
 }
 
 int main() {
     printf("--------Jacobian iteration--------");
     int x,y;
     cout<<"input the number of unknowns x and the number of equations y: ";
     cin>>x>>y;
     if(x!=y) {
         cout<<"invalid number!"<<endl;
         return 0;
     }
     //
     printf("input the augmented matrix:");
     for(int i = 0;i<x;i++){
         for(int j = 0;j<x+1;j++){
             cin>>A[i][j];
         }
     }
     //double re[10];
     cout<<"input the accuracy e and the max number of iterations" ;
     double e;
     int N;
     cin>>e>>N;
     getResult(x,e,N);
     for(int i = 0;i<x;i++){
         cout<<re[i]<<" ";
     }
 }
