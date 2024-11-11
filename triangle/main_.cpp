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

const int edim = 3;
const int ndim = 2;

const vector<int> elem =
{
  0,  1,  3,
  1,  4,  3,
  0,  3,  2,
};
const int esize = elem.size() / edim;// sizeof(elem) / sizeof(int) / edim;

vector<double>  node = // [cm]
  {
    -4.0,   0.0,
    -2.0,   0.0,
    -5.0,   1.7,
    -3.0,   1.7,
    -1.0,   1.7,
  };
const int nsize = node.size() / ndim;

const vector<double> elastic = {20.6, 0.3}; // { [MN/cm2], [1]}

const double thickness = 1.0; // [cm]

vector<double> load = // [MN]
  {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
   -3.0e-3,
    0.0,
   -3.0e-3,
    0.0,
    0.0,
  };

vector<double> bound = { // [cm]
    0.0,
    0.0,
    0.0,
    0.0,
    INFINITY,
    INFINITY,
    INFINITY,
    INFINITY,
    INFINITY,
    INFINITY,
};
// bound[ndim * 2 + 0]

void gradientA(vector<int> P, vector<double> &result) 
{ // ai = eijm xj ym; a0 = x1 y2 - x2 y1; a1 = x2 y0 - x0 y2; a2 = x0 y1 - x1 y0;
  result[0] = node[P[1] * ndim + 0] * node[P[2] * ndim + 1] - node[P[2] * ndim + 0] * node[P[1] * ndim + 1];
//   result[1] = ;
//   result[2] = ;
}

void gradientB(vector<int> P, vector<double> &result) 
{ // bj = eijm 1i ym; b0 = y1 - y2; b1 = y2 - y0; b2 = y0 - y1;
//   result[0] = ;
//   result[1] = ;
//   result[2] = ;
}

void gradientC(vector<int> P, vector<double> &result)
{ // cm = eijm 1i xj; c0 = x2 - x1; c1 = x0 - x2; c2 = x1 - x0;

//   result[0] = ;
//   result[1] = ;
//   result[2] = ;
}

void getDet(vector<double> a, double &det)
{
//   det = ;
}

void blockStressStiffnessPS(double det, double bp, double bq, double cp, double cq, vector<double> &block)
{

  double v = elastic[1];
  double E = elastic[0];

//   block[0 * ndim + 0] = ;
//   block[0 * ndim + 1] = ;

//   block[1 * ndim + 0] = ;
//   block[1 * ndim + 1] = ;
}

void blockStressStiffnessPE(double det, double bp, double bq, double cp, double cq, vector<double> &block)
{
  double v = elastic[1];
  double E = elastic[0];

//   block[0 * ndim + 0] = ;
//   block[0 * ndim + 1] = ;

//   block[1 * ndim + 0] = ;
//   block[1 * ndim + 1] = ;
}

void applyKinematic(int N, vector<double> u, vector<double> f, vector<double> &stiffness)
{
  for (int p = 0; p < N; p++)
  {
    for (int q = 0; q < N; q++)
    {
      //
    }
  }
  printf("\nKinematic boundary condition is applyed\n\n");
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

int sumStiffnessTest(int N, vector<double> stiffness)
{
  const double tol = 1e-9;
  double sum = 0;

// 

  printf("\n\nSum stiffness is %5.9e\n\n", sum);

  return abs(sum) < tol ? 0 : 1;
}

// Преобразовать вектор-матрицу жёсткости в формат CSR https://docs.nvidia.com/cuda/cusparse/index.html#csr-format
void transformMatrixToCsr(int M, int N, int &nz, vector<double> &stiffness, vector<double> &_val, vector<int> &_I, vector<int> &_J)
{

  vector<double> val;
  vector<int> I(N + 1);
  vector<int> J;

  nz = 0;
  int icounter = 0;

  for (int i = 0; i < M; i++)
  {
    // I[icounter] = ;
    icounter++;
    for (int j = 0; j < N; j++)
    {

      if (stiffness[i * N + j] == 0) continue;

    //   val.push_back(//);
    //   J.push_back(//);
      nz++;
    }
  }
  // I[icounter] = 

  _I.resize(M + 1);
  _J.resize(nz);
  _val.resize(nz);

  for (int i = 0; i < nz; i++)
  {
    if (i < M + 1)
    {
      _I[i] = I[i];
    }
    _val[i] = val[i];
    _J[i] = J[i];
  }
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
  int N = nsize * ndim;

  vector<double> globalStiffness(N * N);
  for (int i = 0; i < N * N; i++)
  {
    globalStiffness[i] = 0;
  }

  printStiffness(N, globalStiffness);

  vector<double> a{0,0,0}, b{0,0,0}, c{0,0,0};
  vector<int> P{0,0,0}, Q{0,0,0};         
  vector<double> block(ndim * ndim);

  for (int ie = 0; ie < esize; ie++)
  {
    // P[0] = ;
    // P[1] = ;
    // P[2] = ;

    // Q[0] = ;
    // Q[1] = ;
    // Q[2] = ;

    // double det = 0;
    // gradientA(P, a);
    // gradientB(P, b);
    // gradientC(P, c);
    // getDet(a, det);

    // printABCDet(ie, edim, det, a, b, c);

    // double elemStiffness[ndim * edim * ndim * edim];

    for (int p = 0; p < edim; p++)
    {

      for (int q = 0; q < edim; q++)
      {
        // blockStressStiffnessPE(det, b[p], b[q], c[p], c[q], block);

        // printf("\nPQ[%d%d] pq[%d%d]\n", P[p], Q[q], p, q);
        // printStiffness(ndim, block);

        for (int i = 0; i < ndim; i++)
        {
          for (int j = 0; j < ndim; j++)
          {
            // globalStiffness[...] += block[...];
          }
        }
      }
    }

    // printStiffness(ndim * edim, elemStiffness);
  }

  int status = sumStiffnessTest(N, globalStiffness);

  if (status != 0)
  {
    printf("\nSum stiffness error\n");
    int exit(status);
  }
  else
  {
    printf("\nSum stiffness ok\n");
  }

//   printStiffnessF(N, globalStiffness, load);
//   applyKinematic(N, bound, load, globalStiffness);
//   printStiffnessF(N, globalStiffness, load);

//   status = kinematicTest();

  if (status != 0)
  {
    printf("\nKinematic error\n");
    int exit(status);
  }
  else
  {
    printf("\nKinematic test ok\n");
  }

//   transformMatrixTest();

  vector<double> val;
  vector<int> I;
  vector<int> J;
  int nz;

//   transformMatrixToCsr(N, N, nz, globalStiffness, val, I, J);
//   printSCR(N, N, nz, val, I, J);

  vector<double> result(N);
  for (int i = 0; i < N; i++)
  {
    result[i] = 0.0;
  }

  // status = cgTest();

  // status = cg_host(N, nz, I, J, val, load, result, 1);
  // resultPrint(nsize, ndim, node, result);

  // Cochi equation
  // Constutive equation
  vector<double> epsilon(3);
  vector<double> sigma(3);

  double v = elastic[1];
  double E = elastic[0];
  // double la = 
  // double mu = 
  // double ka = 

  printTensorH();
  for (int ie = 0; ie < esize; ie++)
  {



    for (int p = 0; p < edim; p++)
    {
      // epsilon[0] += 
      // epsilon[1] += 
      // epsilon[2] += 

      // plane strain
      // sigma[0] += 
      // sigma[1] += 
      // sigma[2] += 

      // plane stress
      // sigma[0] += 
      // sigma[1] += 
      // sigma[2] += 
    }
    // printTensor(ie, epsilon, sigma);
  }

  std::cin.get();
  int exit(status);
}