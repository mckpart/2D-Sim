#include <iostream>
#include <vector> 

using namespace std;

void calc_xy_dens(double x, double y){
   int ind_1 = 0;
   int ind_2 = 0;
   double boxLength = 2.5; 
   double cell_L = .25; 
   int val = boxLength/cell_L + 1; 
   vector<vector<double>> F(val,vector<double>(val,0)); 

   ind_1 = (x + boxLength/2)/cell_L; // cell_L = delta x = delta y 
   ind_2 = (y + boxLength/2)/cell_L;

   if(x > (ind_1 + 0.5) * cell_L - boxLength/2){
      ++ind_1;
   }
   if(y > (ind_2 + 0.5) * cell_L - boxLength/2){
      ++ind_2;
   }
   std::cout << "ind1 " << ind_1 << " ind2 " << ind_2 << std::endl;
   F[ind_1][ind_2] = F[ind_1][ind_2] + 1; 

   for(int k = 0; k < val; ++k){
      for(int n = 0; n < val; ++n){
         std::cout << F[k][n] << " ";
      }
      std::cout << "\n"; 
   }
 }

void attemptCorr(){
   double L = 10; 
   double cell_L = 1; 
   int val = L/cell_L + 1; 
   int x = 0; 
   int y = 0; 
   
   int ind1 = 0; 
   int ind2 = 0; 
   vector<vector<double>> M = {{2,1},{-1,3},{4,2},{-3,2}}; 
   vector<vector<double>> F(val,vector<double>(val,0)); 
   
   for(int k = 0; k < 4; k++){
      x = M[k][0]; 
      y = M[k][1];
//      std::cout << x << " " << y << std::endl
      ind1 = x + L/2; // delta x = delta y = 1  
      ind2 = y + L/2; 
      std::cout << "ind1: " << ind1
	        << " ind2: " << ind2 << "\n"; 
      F[ind1][ind2] = F[ind1][ind2] + 1; 
   }
   for(int k = 0; k < val; k++){
      for(int n = 0; n < val; n++){
         std::cout << F[k][n] << " "; 
      }
      std::cout <<"\n"; 
   }
}


int main(){
   int val = 3; 

//   std::vector<std::vector<double>> matrix(val,std::vector<double>(val,0)); 
//   for(int k = 0; k < val; ++k){
//      for(int n = 0; n < val; ++n){
//         matrix[k][n] = k + 1; 
//      }
//   }
//   for(int k = 0; k < val; ++k){
//      for(int n = 0; n < val; ++n){
//         std::cout << matrix[n][k] << std::endl;
//      }
//   }
//attemptCorr(); 
double x = .376;
double y = .13; 

calc_xy_dens(x,y); 
}
