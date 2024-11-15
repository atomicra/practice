Реализация модели статической задачи линейной теории упругости методом конечных элементов в стеке с++, c.
Литературу см. ./docs

Задание по компьютерному практикуму для студентов 6 курса кафедры механики композитов:
- реализовать недостающие процедуры/методы для расчёта статической задачи линейной теории упругости с использованием метода конечных элементов;
- учесть случай плоских напряжений/(плоских деформаций);
- после отладки и реализации программы необходимо выполнить расчёт модели в соответствии с заданием;
- (по согласованию) учесть массовые, поверхностные силы, поле тепловых деформаций;
- результатом выполненной работы является демонстрация рабочего кода проекта и вывод вектора перемещений по всем узлам модели, 
  при условии что вектор перемещений совпадёт с эталонным расчётом;
- (по согласованию) определить деформации и напряжения по элементу;  
- за основу взять курс видео лекций по МКЭ https://youtu.be/_vsjFA28IL8

Список заданий и схемы сеток см. assignment.pdf
 
Проект содержит исходную модель состоящую из
- узлов конечных элементов node, 
- списка элементов elem, 
- вектора узловых сил load, 
- вектора кинематических граничных условий bound, 
- толщину и материальные константы.

Приложения и плагины, которые необходимо установить:
1. https://code.visualstudio.com/
2. https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools
3. https://code.visualstudio.com/docs/languages/cpp (компилятор и дебаггер g++)

логи для написания, отладки и сравнения результатов работы программы представлены ниже:

//--------------------------------------------- случай плоских напряжений ---------------------------------------
matrix 10x10
0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000

elem 0
det 3.400
a               b               c
-3.400          -1.700          -1.000
6.800           1.700           -1.000
0.000           0.000           2.000

PQ[00] pq[00]
matrix 2x2
10.7860     3.6786
3.6786      6.6963


PQ[01] pq[01]
matrix 2x2
-8.4557     -0.2830
0.2830      -0.0383


PQ[03] pq[02]
matrix 2x2
-2.3303     -3.3956
-3.9615     -6.6580


PQ[10] pq[10]
matrix 2x2
-8.4557     0.2830
-0.2830     -0.0383


PQ[11] pq[11]
matrix 2x2
10.7860     -3.6786
-3.6786     6.6963


PQ[13] pq[12]
matrix 2x2
-2.3303     3.3956
3.9615      -6.6580


PQ[30] pq[20]
matrix 2x2
-2.3303     -3.9615
-3.3956     -6.6580


PQ[31] pq[21]
matrix 2x2
-2.3303     3.9615
3.3956      -6.6580


PQ[33] pq[22]
matrix 2x2
4.6606      0.0000
0.0000      13.3161



elem 1
det 3.400
a               b               c
3.400           0.000           -2.000
3.400           1.700           1.000
-3.400          -1.700          1.000

PQ[11] pq[00]
matrix 2x2
4.6606      -0.0000
-0.0000     13.3161


PQ[14] pq[01]
matrix 2x2
-2.3303     -3.9615
-3.3956     -6.6580


PQ[13] pq[02]
matrix 2x2
-2.3303     3.9615
3.3956      -6.6580


PQ[41] pq[10]
matrix 2x2
-2.3303     -3.3956
-3.9615     -6.6580


PQ[44] pq[11]
matrix 2x2
10.7860     3.6786
3.6786      6.6963


PQ[43] pq[12]
matrix 2x2
-8.4557     -0.2830
0.2830      -0.0383


PQ[31] pq[20]
matrix 2x2
-2.3303     3.3956
3.9615      -6.6580


PQ[34] pq[21]
matrix 2x2
-8.4557     0.2830
-0.2830     -0.0383


PQ[33] pq[22]
matrix 2x2
10.7860     -3.6786
-3.6786     6.6963



elem 2
det 3.400
a               b               c
3.400           0.000           -2.000
6.800           1.700           1.000
-6.800          -1.700          1.000

PQ[00] pq[00]
matrix 2x2
4.6606      -0.0000
-0.0000     13.3161


PQ[03] pq[01]
matrix 2x2
-2.3303     -3.9615
-3.3956     -6.6580


PQ[02] pq[02]
matrix 2x2
-2.3303     3.9615
3.3956      -6.6580


PQ[30] pq[10]
matrix 2x2
-2.3303     -3.3956
-3.9615     -6.6580


PQ[33] pq[11]
matrix 2x2
10.7860     3.6786
3.6786      6.6963


PQ[32] pq[12]
matrix 2x2
-8.4557     -0.2830
0.2830      -0.0383


PQ[20] pq[20]
matrix 2x2
-2.3303     3.3956
3.9615      -6.6580


PQ[23] pq[21]
matrix 2x2
-8.4557     0.2830
-0.2830     -0.0383


PQ[22] pq[22]
matrix 2x2
10.7860     -3.6786
-3.6786     6.6963

Sum stiffness is -3.552713679e-15
Sum stiffness ok

matrix 10x10
15.4467     3.6786      -8.4557     -0.2830     -2.3303     3.9615      -4.6606     -7.3571     0           0           |0.0000
3.6786      20.0124     0.2830      -0.0383     3.3956      -6.6580     -7.3571     -13.3161    0           0           |0.0000
-8.4557     0.2830      15.4467     -3.6786     0           0           -4.6606     7.3571      -2.3303     -3.9615     |0.0000
-0.2830     -0.0383     -3.6786     20.0124     0           0           7.3571      -13.3161    -3.3956     -6.6580     |0.0000
-2.3303     3.3956      0           0           10.7860     -3.6786     -8.4557     0.2830      0           0           |0.0000
3.9615      -6.6580     0           0           -3.6786     6.6963      -0.2830     -0.0383     0           0           |-0.0030
-4.6606     -7.3571     -4.6606     7.3571      -8.4557     -0.2830     26.2327     -0.0000     -8.4557     0.2830      |0.0000
-7.3571     -13.3161    7.3571      -13.3161    0.2830      -0.0383     -0.0000     26.7088     -0.2830     -0.0383     |-0.0030
0           0           -2.3303     -3.3956     0           0           -8.4557     -0.2830     10.7860     3.6786      |0.0000
0           0           -3.9615     -6.6580     0           0           0.2830      -0.0383     3.6786      6.6963      |0.0000


Kinematic boundary condition is applied

matrix 10x10
15.4467     0           0           0           0           0           0           0           0           0           |0.0000
0           20.0124     0           0           0           0           0           0           0           0           |0.0000
0           0           15.4467     0           0           0           0           0           0           0           |0.0000
0           0           0           20.0124     0           0           0           0           0           0           |0.0000
0           0           0           0           10.7860     -3.6786     -8.4557     0.2830      0           0           |0.0000
0           0           0           0           -3.6786     6.6963      -0.2830     -0.0383     0           0           |-0.0030
0           0           0           0           -8.4557     -0.2830     26.2327     -0.0000     -8.4557     0.2830      |0.0000
0           0           0           0           0.2830      -0.0383     -0.0000     26.7088     -0.2830     -0.0383     |-0.0030
0           0           0           0           0           0           -8.4557     -0.2830     10.7860     3.6786      |0.0000
0           0           0           0           0           0           0.2830      -0.0383     3.6786      6.6963      |0.0000

matrix 5x5
55.0000     -46.0000    4.0000      0           0           |500.0000
-46.0000    140.0000    -46.0000    0           0           |2000.0000
4.0000      -46.0000    110.0000    -46.0000    4.0000      |1000.0000
0           0           -46.0000    142.0000    -46.0000    |2000.0000
0           0           4.0000      -46.0000    65.0000     |900.0000


Kinematic boundary condition is applied

matrix 5x5
55.0000     0           0           0           0           |8250.0000
0           140.0000    -46.0000    0           0           |8900.0000
0           -46.0000    110.0000    -46.0000    0           |240.0000
0           0           -46.0000    142.0000    0           |3840.0000
0           0           0           0           65.0000     |2600.0000


Kinematic test ok

nz=9
val[9]
I[5]
J[9]
val       I         J
1.000     0         0
4.000     2         1
2.000     4         1
3.000     7         2
5.000     9         0
7.000               3
8.000               4
9.000               2
6.000               4



nz=32
val[32]
I[11]
J[32]
val       I         J
15.447    0         0
20.012    1         1
15.447    2         2
20.012    3         3
10.786    4         4
-3.679    8         5
-8.456    12        6
0.283     18        7
-3.679    24        4
6.696     28        5
-0.283    32        6
-0.038              7
-8.456              4
-0.283              5
26.233              6
-0.000              7
-8.456              8
0.283               9
0.283               4
-0.038              5
-0.000              6
26.709              7
-0.283              8
-0.038              9
-8.456              6
-0.283              7
10.786              8
3.679               9
0.283               6
-0.038              7
3.679               8
6.696               9


> GPU device has 5 Multi-Processors, SM 5.0 compute capabilities

iteration =   1, residual = 2.620849e-03
iteration =   2, residual = 1.363684e-03
iteration =   3, residual = 8.892826e-04
iteration =   4, residual = 8.179503e-04
iteration =   5, residual = 4.207631e-04
iteration =   6, residual = 1.818915e-18
Test Summary:  Error amount = 0.000000


i      x0        x1        u0             u1
0     -4.000     0.000     0.000e+00      0.000e+00
1     -2.000     0.000     0.000e+00      0.000e+00
2     -5.000     1.700    -3.693e-04     -6.595e-04
3     -3.000     1.700    -1.879e-04     -1.112e-04
4     -1.000     1.700    -1.880e-04      1.106e-04

// вывод деформаций и напряжений по элементам

elem           e11            e22            g12            s11            s22            s12
0              0.000e+00     -6.540e-05     -1.106e-04     -4.442e-04     -1.481e-03     -8.759e-04     
1             -9.477e-09     -1.845e-07      3.149e-07     -1.468e-06     -4.242e-06      2.495e-06     
2              9.070e-05     -2.267e-04      1.102e-04      5.138e-04     -4.515e-03      8.734e-04 

// с учётом тепловых перемещений
i      x0        x1        u0             u1
0     -4.000     0.000     0.00000e+00    0.00000e+00
1     -2.000     0.000     0.00000e+00    0.00000e+00
2     -5.000     1.700    -2.62801e-03    5.97292e-04   
3     -3.000     1.700    -1.87941e-04    2.43042e-03
4     -1.000     1.700     2.07072e-03    1.36733e-03

// с учётом тепловых деформаций и напряжений
elem           e11            e22            g12            s11            s22            s12
0              0.000e+00      1.429e-03     -1.105e-04     -2.355e-02     -8.906e-04     -8.759e-04     
1              1.129e-03      1.117e-03      2.221e-05     -1.035e-04     -2.992e-04      1.760e-04
2              1.220e-03      8.902e-04      8.831e-05      4.117e-04     -4.810e-03      6.999e-04

//--------------------------------------------- случай плоских деформаций ---------------------------------------


 0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000     
 0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000     
 0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
 0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000     
 0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000     
 0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000     
 0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000
 0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000     
 0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000     
 0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000      0.0000     



elem 0 det 3.400
a               b               c
-3.400         -1.700         -1.000
 6.800          1.700         -1.000         
 0.000          0.000          2.000

PQ[00] pq[00]
 12.9507     4.9519
 4.9519      7.4454


PQ[01] pq[01]
-10.6204     0.9904     
-0.9904      0.7107


PQ[03] pq[02]
-2.3303     -5.9423
-3.9615     -8.1561


PQ[10] pq[10]
-10.6204    -0.9904     
 0.9904      0.7107


PQ[11] pq[11]
 12.9507    -4.9519
-4.9519      7.4454


PQ[13] pq[12]
-2.3303      5.9423
 3.9615     -8.1561


PQ[30] pq[20]
-2.3303     -3.9615
-5.9423     -8.1561


PQ[31] pq[21]
-2.3303      3.9615
 5.9423     -8.1561


PQ[33] pq[22]
 4.6606      0.0000
 0.0000      16.3122



elem 1 det 3.400
a               b               c
 3.400          0.000         -2.000
 3.400          1.700          1.000
-3.400         -1.700          1.000

PQ[11] pq[00]
 4.6606      -0.0000    
 -0.0000     16.3122


PQ[14] pq[01]
-2.3303     -3.9615
-5.9423     -8.1561


PQ[13] pq[02]
-2.3303      3.9615
 5.9423     -8.1561


PQ[41] pq[10]
-2.3303     -5.9423
-3.9615     -8.1561


PQ[44] pq[11]
 12.9507     4.9519
 4.9519      7.4454


PQ[43] pq[12]
-10.6204     0.9904
-0.9904      0.7107


PQ[31] pq[20]
-2.3303      5.9423
 3.9615     -8.1561


PQ[34] pq[21]
-10.6204    -0.9904     
 0.9904      0.7107


PQ[33] pq[22]
 12.9507    -4.9519
-4.9519      7.4454



elem 2 det 3.400
a               b               c
 3.400          0.000         -2.000
 6.800          1.700          1.000
-6.800         -1.700          1.000

PQ[00] pq[00]
 4.6606      -0.0000
 -0.0000     16.3122


PQ[03] pq[01]
-2.3303     -3.9615
-5.9423     -8.1561     


PQ[02] pq[02]
-2.3303      3.9615
 5.9423     -8.1561


PQ[30] pq[10]
-2.3303     -5.9423
-3.9615     -8.1561


PQ[33] pq[11]
 12.9507     4.9519
 4.9519      7.4454


PQ[32] pq[12]
-10.6204     0.9904
-0.9904      0.7107


PQ[20] pq[20]
-2.3303      5.9423
 3.9615     -8.1561


PQ[23] pq[21]
-10.6204    -0.9904
 0.9904      0.7107


PQ[22] pq[22]
 12.9507    -4.9519
-4.9519      7.4454     



Sum stiffness is -7.105427358e-15


Sum stiffness ok
matrix 10x10
 17.6114     4.9519     -10.6204     0.9904     -2.3303      3.9615     -4.6606     -9.9038      0           0          | 0.000
 4.9519      23.7576    -0.9904      0.7107      5.9423     -8.1561     -9.9038     -16.3122     0           0          | 0.000
-10.6204    -0.9904      17.6114    -4.9519      0           0          -4.6606      9.9038     -2.3303     -3.9615     | 0.000
 0.9904      0.7107     -4.9519      23.7576     0           0           9.9038     -16.3122    -5.9423     -8.1561     | 0.000
-2.3303      5.9423      0           0           12.9507    -4.9519     -10.6204    -0.9904      0           0          | 0.000
 3.9615     -8.1561      0           0          -4.9519      7.4454      0.9904      0.7107      0           0          |-0.003
-4.6606     -9.9038     -4.6606      9.9038     -10.6204     0.9904      30.5621     0          -10.6204    -0.9904     | 0.000
-9.9038     -16.3122     9.9038     -16.3122    -0.9904      0.7107      0           31.2029     0.9904      0.7107     |-0.003
 0           0          -2.3303     -5.9423      0           0          -10.6204     0.9904      12.9507     4.9519     | 0.000         
 0           0          -3.9615     -8.1561      0           0          -0.9904      0.7107      4.9519      7.4454     | 0.000


Kinematic boundary condition is applyed

matrix 10x10
 17.6114     0           0           0           0           0           0           0           0           0          | 0.000
 0           23.7576     0           0           0           0           0           0           0           0          | 0.000
 0           0           17.6114     0           0           0           0           0           0           0          | 0.000
 0           0           0           23.7576     0           0           0           0           0           0          | 0.000
 0           0           0           0           12.9507    -4.9519     -10.6204    -0.9904      0           0          | 0.000
 0           0           0           0          -4.9519      7.4454      0.9904      0.7107      0           0          |-0.003
 0           0           0           0          -10.6204     0.9904      30.5621     0          -10.6204    -0.9904     | 0.000
 0           0           0           0          -0.9904      0.7107      0           31.2029     0.9904      0.7107     |-0.003
 0           0           0           0           0           0          -10.6204     0.9904      12.9507     4.9519     | 0.000
 0           0           0           0           0           0          -0.9904      0.7107      4.9519      7.4454     | 0.000

matrix 5x5
 55.0000    -46.0000     4.0000      0           0          | 500.000
-46.0000     140.0000   -46.0000     0           0          | 2000.000      
 4.0000     -46.0000     110.0000   -46.0000     4.0000     | 1000.000
 0           0          -46.0000     142.0000   -46.0000    | 2000.000      
 0           0           4.0000     -46.0000     65.0000    | 900.000


Kinematic boundary condition is applyed

matrix 5x5
 55.0000     0           0           0           0          | 500.000
 0           140.0000   -46.0000     0           0          | 2000.000
 0          -46.0000     110.0000   -46.0000     0          | 1000.000      
 0           0          -46.0000     142.0000    0          | 2000.000
 0           0           0           0           65.0000    | 900.000       


Kinematic error

nz=9 I[5] J[9] val[9]
I         J         val
0         0          1.000         
2         1          4.000
4         1          2.000
7         2          3.000
9         0          5.000
          3          7.000         
          4          8.000
          2          9.000
          4          6.000



nz=30 I[11] J[30] val[30]
I         J         val       
0         0          17.611
1         1          23.758
2         2          17.611
3         3          23.758
4         4          12.951        
8         5         -4.952
12        6         -10.620
17        7         -0.990
22        4         -4.952
26        5          7.445
30        6          0.990
          7          0.711
          4         -10.620
          5          0.990
          6          30.562
          8         -10.620
          9         -0.990
          4         -0.990         
          5          0.711
          7          31.203
          8          0.990
          9          0.711
          6         -10.620        
          7          0.990
          8          12.951
          9          4.952
          6         -0.990
          7          0.711
          8          4.952
          9          7.445


iteration =   1, residual = 2.678562e-03, r1/r0 = 3.985943e-01
iteration =   2, residual = 1.413111e-03, r1/r0 = 2.783230e-01 
iteration =   3, residual = 1.455495e-03, r1/r0 = 1.060885e+00
iteration =   4, residual = 7.540911e-04, r1/r0 = 2.684271e-01
iteration =   5, residual = 4.762667e-04, r1/r0 = 3.988897e-01 
iteration =   6, residual = 1.308503e-18, r1/r0 = 7.548302e-30

total time for CG without preconditioning is 0.034000

Test Summary:  Error amount = 1.206175e-18


i      x0        x1        u0             u1
0     -4.000     0.000     0.000e+00      0.000e+00
1     -2.000     0.000     0.000e+00      0.000e+00     
2     -5.000     1.700    -3.871e-04     -6.291e-04
3     -3.000     1.700    -1.702e-04     -9.095e-05
4     -1.000     1.700    -1.707e-04      9.958e-05  

// вывод деформаций и напряжений по элементам

elem           e11            e22            g12            s11            s22            s12
0              0.000e+00     -5.350e-05     -1.001e-04     -6.358e-04     -1.484e-03     -7.933e-04     
1             -2.460e-07      2.538e-06     -5.009e-06      2.334e-05      6.747e-05     -3.969e-05     
2              1.084e-04     -2.118e-04      1.051e-04      4.900e-04     -4.584e-03      8.330e-04  

// с учётом тепловых перемещений
i      x0        x1        u0             u1
0     -4.000     0.000     0.00000e+00    0.00000e+00
1     -2.000     0.000     0.00000e+00    0.00000e+00   
2     -5.000     1.700    -3.31927e-03    9.97640e-04
3     -3.000     1.700    -1.70223e-04    3.37890e-03
4     -1.000     1.700     2.76147e-03    1.72628e-03   

// с учётом тепловых деформаций и напряжений
elem           e11            e22            g12            s11            s22            s12
0              0.000e+00      1.987e-03     -1.001e-04     -3.457e-02     -3.078e-03     -7.933e-04     
1              1.465e-03      1.501e-03     -6.418e-05      2.991e-04      8.645e-04     -5.085e-04     
2              1.574e-03      1.287e-03      1.643e-04      7.658e-04     -3.787e-03      1.302e-03