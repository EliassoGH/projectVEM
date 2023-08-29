#ifndef __tetrahedron_jaskowiec_rule_
#define __tetrahedron_jaskowiec_rule_

# include <cmath>
# include <cstdlib>
# include <iomanip>
# include <iostream>

using namespace std;

void comp_next ( int n, int k, int a[], bool &more, int &h, int &t );
double *monomial_value ( int m, int n, int e[], double x[] );
double r8mat_det_4d ( double a[] );
void r8vec_copy ( int n, double a1[], double a2[] );
int rule_order ( int p );
void rule00 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule01 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule02 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule03 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule04 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule05 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule06 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule07 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule08 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule09 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule10 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule11 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule12 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule13 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule14 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule15 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule16 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule17 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule18 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule19 ( int n, double a[], double b[], double c[], double d[], double w[] );
void rule20 ( int n, double a[], double b[], double c[], double d[], double w[] );
void tetrahedron_jaskowiec_rule ( int p, int n, double a[], double b[], 
  double c[], double d[], double w[] );
double tetrahedron_unit_monomial_integral ( int expon[] );
double tetrahedron_unit_volume ( );
double tetrahedron_volume ( double tetra[4*3] );

#endif // __tetrahedron_jaskowiec_rule_