#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <cctype>
#include <functional> // std::plus // isinf()

#include <stdlib.h>

#include "../helper.h"

using namespace std;

const int schemaNO = 0;

const int edim = 2;
const int ndim = 2;

const vector<int> elem =
{
  2, 0,
  0, 1,
  1, 2
};

const int esize = elem.size() / edim;

vector<double> node =
{// [cm]
  0, 0,
500, 0,
900, 692.82032
};

const int nsize = node.size() / ndim;

const vector<double> elastic = {20e+6, 11e-6, 5};// E, alpha, theta
const vector<double> area = 
{
  10,
  20, 
  20
};

vector<double> load =
{ // [N]
  0.0,
  0.0,
  0.0,
  0.0,
  0.0,
  -100.0e+3
};

vector<double> bound =
{ // [cm]
  0.0,
  0.0,
  INFINITY,
  0.0,
  INFINITY,
  INFINITY
};

void getLRG(vector<double> a, vector<double> b, double &L, vector<double> &r, vector<double> &g)
{

}

void applyKinematic(int N, vector<double> u, vector<double> &f, vector<double> &stiffness)
{

  printf("\nKinematic boundary condition is applyed\n\n");
}

void blockStiffness(double kx, vector<double> gradG, vector<double> &block)
{

}

int kinematicTest()
{ // Л.Сегерлинд. Применение метода конечных элементов стр. 110-112
  vector<double> u = {
      150,
      INFINITY,
      INFINITY,
      INFINITY,
      40,
  };
  vector<double> f = {
      500,
      2000,
      1000,
      2000,
      900,
  };
  vector<double> stiffness = {
      55, -46,   4,   0,   0,
     -46, 140, -46,   0,   0,
       4, -46, 110, -46,   4,
       0,   0, -46, 142, -46,
       0,   0,   4, -46,  65,
  };
  const int N = u.size();
  printStiffnessF(N, stiffness, f);
  applyKinematic(N, u, f, stiffness);
  printStiffnessF(N, stiffness, f);

  double sum =
      stiffness[0 * N + 0] - 55.0 + f[0] - 8250 +
      stiffness[1 * N + 1] - 140.0 + stiffness[1 * N + 2] + 46.0 + f[1] - 8900 +
      stiffness[2 * N + 1] + 46.0 + stiffness[2 * N + 2] - 110.0 + stiffness[2 * N + 3] + 46.0 + f[2] - 240 +
      stiffness[3 * N + 2] + 46.0 + stiffness[3 * N + 3] - 142.0 + f[3] - 3840 +
      stiffness[4 * N + 4] - 65.0 + f[4] - 2600;

  return sum == 0.0 ? 0 : 1;
}

// Преобразовать вектор-матрицу жёсткости в формат CSR https://docs.nvidia.com/cuda/cusparse/index.html#csr-format
void transformMatrixToCsr(int M, int N, int &nz, vector<double> &stiffness, vector<double> &val, vector<int> &I, vector<int> &J)
{

}

void transformMatrixTest()
{
  const int M = 4;
  const int N = 5;

  vector<double> stiffness = {
      1.0, 4.0, 0.0, 0.0, 0.0,
      0.0, 2.0, 3.0, 0.0, 0.0,
      5.0, 0.0, 0.0, 7.0, 8.0,
      0.0, 0.0, 9.0, 0.0, 6.0,
  };

  vector<double> val;
  vector<int> I;
  vector<int> J;
  int nz;
  transformMatrixToCsr(M, N, nz, stiffness, val, I, J);
  printSCR(M, N, nz, val, I, J);
}

int main(int argc, char **argv)
{
  double E = elastic[0];
  double al = elastic[1];
  double th = elastic[2];

  int N = nsize * ndim;

  vector<double> globalStiffness(N * N);
  for (int i = 0; i < N * N; i++)
  {
    globalStiffness[i] = 0;
  }

  printStiffness(N, globalStiffness);

  vector<int> P{0, 0}, Q{0, 0};
  vector<double> r{0, 0};
  vector<double> g{0, 0, 0, 0};
  const vector<double> gradV{-1, 1};  
  const vector<double> gradG{1, -1, -1, 1};

  vector<double> a{0, 0};
  vector<double> b{0, 0};   

  vector<double> elemStiffness(edim * edim);

  for (int ie = 0; ie < esize; ie++)
  {
    const double A = area[ie];

    // P[0] = 
    // P[1] = 

    // Q[0] = 
    // Q[1] = 

    double L = 1;
    // a[0] = 
    // a[1] = 
    // b[0] = 
    // b[1] =   

    getLRG(a, b, L, r, g);

    const double kx = E * A / L;

    blockStiffness(kx, gradG, elemStiffness);
    printf("ie=%d\n", ie);    
    printf("\nL=%lf; ri=(%lf, %lf)\n", L, r[0], r[1]); 
    printf("\ngij=\n"); 
    printStiffness3(ndim, g);
    printf("Kpq=\n"); 
    printStiffness(edim, elemStiffness);


    vector<double> elemStiffness2(edim *  edim * ndim * ndim);
    for (int p = 0; p < edim; p++)
    {
      // add termal force
      // here
      // here

      for (int q = 0; q < edim; q++)
      {
        for (int i = 0; i < ndim; i++)
        {
          for (int j = 0; j < ndim; j++)
          {
            // elemStiffness2[...] += elemStiffness[...] * g[...];
            // globalStiffness[...] += elemStiffness[...] * g[...];
          }
        }
      }
    }
    printf("K^ij_pq=\n");
    printStiffness(edim * ndim, elemStiffness2);   
  }

  printStiffnessF(N, globalStiffness, load);
  applyKinematic(N, bound, load, globalStiffness);
  printStiffnessF(N, globalStiffness, load);

  int status = kinematicTest();

  if (status != 0)
  {
    printf("\nKinematic error\n");
    int exit(status);
  }
  else
  {
    printf("\nKinematic test ok\n");
  }

  // transformMatrixTest();

  vector<double> val;
  vector<int> I;
  vector<int> J;
  int nz;

  // transformMatrixToCsr(N, N, nz, globalStiffness, val, I, J);
  // printSCR(N, N, nz, val, I, J);

  vector<double> result(N);
  for (int i = 0; i < N; i++)
  {
    result[i] = 0.0;
  }

  // status = cg_host(N, nz, I, J, val, load, result, 1);

  resultPrint(nsize, ndim, node, result);

  vector<double> aforces(esize);
  for (int ie = 0; ie < esize; ie++)
  {
    const double A = area[ie];

    // Q[0] =
    // Q[1] =

    double L = 1;
    // a[0] = 
    // a[1] = 
    // b[0] = 
    // b[1] =    

    getLRG(a, b, L, r, g);

    vector<double> dN{-1/L, 1/L};  
    // add axial force calculation here 
  }

  axialForcesPrint(aforces);

  std::cin.get();
  int exit(status);
}