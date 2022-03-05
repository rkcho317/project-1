//Rosa Cho, 888244357
//Stephen Merwin, 887500593

//============================================================
// file: main.cpp
//============================================================
#include <iostream>
#include <cstring>
#include <initializer_list>
#include <cassert>

//MATRIX and VECTOR classes assignment
#include "vector3d_T.h"
#include "matrix3d_T.h"

//QUATERNION class
#include "quaternion_T.h"

#define _USE_MATH_DEFINES 
#include <cmath>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
    #define M_PI_2 3.14159265358979323846
#endif

template <typename T>
void print(T v) {
std::cout << v << std::endl;
}

template <typename T>
void show_vect(T v) { std::cout << v.name() << " is: " << v << "\n"; }

template <typename T>
void show_mat(T m) { std::cout << m.name() << " is: " << m << "\n"; }

void plane_rotation(const std::string& msg, const quatD& plane, const std::initializer_list<double>& li) {
 matrix3dD rotate = matrix3dD("rot_matrix", 3, li);
 assert(plane.rot_matrix() == rotate);
 std::cout << msg << " is: " << plane << plane.rot_matrix() << "\n";
}

std::string yes_or_no(bool condition) { return condition ? "YES" : "no"; }


void test_vectors() {
  print("\n====================  TESTING VECTORS  ========================");
  vector3d<double> u("u", 3, {1,  2,  4});
    std::cout << u.name() << "\n";
    std::cout << u << "\n";
    u.zero();
    show_vect(u);
  vector3dD v("v", 3, {8, 16, 32});
  vector3dD i("i", 3, {1, 0, 0}), j("j", 3, {0, 1, 0}), k("k", 3, {0, 0, 1});
  vector3dD w(3 * i + 4 * j - 2 * k);

  show_vect(u);
  show_vect(v);
  show_vect(i);
  show_vect(j);
  show_vect(k);
  j + k;
  show_vect(w);

  assert(u == u);
  assert(u != v);
  
  assert(u + v == v + u);
  assert(u - v == -(v - u));
  assert(-(-u) == u);
  
  assert(3.0 + u == u + 3.0);
  assert(3.0 * u == u * 3.0);
  assert((u - 3.0) == -(3.0 - u));
  assert((5.0 * u) / 5.0 == u);

  assert(u + vector3dD::zero() == u);

  assert((i.dot(j) == j.dot(k)) == (k.dot(i) == 0));
  
  assert(i.cross(j) == k);
  assert(j.cross(k) == i);
  assert(k.cross(i) == j);
  
  assert(u.cross(v) == -v.cross(u));
  assert(u.cross(v + w) == u.cross(v) + u.cross(w));
  assert((u.cross(v)).dot(u) == 0);

  print(i.angle(j));
  print(M_PI/2);

  assert(i.angle(j) == M_PI_2);
  assert(j.angle(k) == M_PI_2);
  assert(k.angle(i) == M_PI_2);
  
  vector3D uhat = u / u.magnitude(); // unit vector in u direction
  show_vect(u);
  show_vect(uhat);
  print(uhat.magnitude());
  assert(uhat.magnitude() - 1.0 < 1.0e-10);

  print("...test vectors assertions passed");
  print("====================  FINISHED testing vectors  ========================");
}

void test_matrices() {
  print("\n====================  TESTING MATRICES  ========================");
  matrix3dD a("a", 3, {3, 2, 0,   0, 0, 1,   2, -2, 1});
  matrix3dD b("b", 3, {1, 0, 5,   2, 1, 6,   3,  4, 0});
  
  matrix3dD id = matrix3dD::identity(3);
  assert(a * id == a);
  assert(a * b != -b * a);
  assert((a * b).transpose() == b.transpose() * a.transpose());

  matrix3dD acopy(a);    // copy constructor
  matrix3dD a2copy = a;  // copy constructor

  matrix3dD bcopy;
  bcopy = b;        // assignment operator

  matrix3dD ainv = a.inverse();
  matrix3dD binv = b.inverse();

  show_mat(a);
  show_mat(b);
  show_mat(-a);
  show_mat(-b);
  show_mat(a * b);
  printf("|a| = %.2f\n", a.determinant());
  printf("|b| = %.2f\n", b.determinant());
  show_mat(a.transpose());
  show_mat(b.transpose());

  show_mat(a.minors());
  show_mat(b.minors());

  show_mat(a.cofactor());
  show_mat(b.cofactor());
  
  print("A adjugate is: ");
  show_mat(a.adjugate());
  print("B adjugate is: ");
  show_mat(b.adjugate());

  print("A Inverse");
  show_mat(ainv);
  print("B Inverse");
  show_mat(binv);
  print("A * AInv");
  show_mat(a * ainv);
  print("B * Binv");
  show_mat(b * binv);
  show_mat(matrix3dD::identity(3));

  assert(a * ainv == matrix3dD::identity(3));
  assert(a * ainv == ainv * a);
  assert(b * binv == matrix3dD::identity(3));
  assert(b * binv == binv * b);
  assert(a.transpose().transpose() == a);
  assert(a.determinant() == a.transpose().determinant());

  assert(a + b == b + a);
  assert(a - b == -(b - a));
  assert(3.0 + a == a + 3.0);
  assert(3.0 * a == a * 3.0);
  assert((a + 3.0) - 3.0 == a);
  assert((3.0 * a) / 3.0 == a);
  assert(-(-a) == a);

  matrix3dD zerod("zerod", 3, {1, 2, 3,   4, 5, 6,   7, 8, 9});
  assert(zerod.determinant() == 0);

  print("...test matrices assertions passed");
  print("====================  FINISHED testing matrices  ========================");
}

void test_matrices_and_vectors() {
  print("\n====================  TESTING MATRICES and VECTORS  ========================");
  vector3dD p("p", 2, {1, 2});
  matrix3dD m("m", 2, {1, 2,   3, 4});
  show_vect(p);
  show_mat(m);
  assert(p * m == m * p);

  vector3dD q("q", 3, {1, 2, 3});
  matrix3dD n("n", 3, {1, 2, 3,   4, 5, 6,   7, 8, 9});
  show_vect(q);
  show_mat(n);
  assert(q * n == n * q);
  print("...test_matrices_and_vectors assertions passed");
  print("====================  FINISHED testing matrices and vectors  ========================");
}

void test_quarternion(){
    print("\n====================  TESTING QUARTERNIONS  ========================");
     quatD a(1, 2, 3, 4), b(4, 0, 0, 7), c(0, 1, 1, 0), d(0, 0, 1, 0);
     quatD e(0, 0, 0, 1), f(0, 0, 0, 0), g(1, 0, 0, 0), h(3, 0, 0, 0);

    std::cout << "a = " << a << ")\nb = " << b << ")\nc = " << c << ")\nd = " << d
           << ")\ne = " << e << ")\nf = " << f << ")\ng = " << g << ")\nh = " <<  h << "\n";

    std::cout << "c + d = " <<  c + d << "\nc + d + e = " << c + d + e << "\n";
    std::cout << "5 * h = " << 5 * h << "\nh * 5 = " << h * 5 << "\nh / 3.0 = " << h / 3.0 << "\n\n";

    std::cout << "h.magnitude() is " << h.magnitude() << "\nh.unit() is " << h.unit();
    std::cout << "g.unit() is " << g.unit() << "\na.unit() is " << a.unit() << ")\n\n";

    std::cout << "a.vector() is " << a.vector() << "\na.scalar() is " << a.scalar() << "\n";
    std::cout << "a.conjugate() is " << a.conjugate() << "\na.inverse() is " << a.inverse()
              << "\na * a.inverse() is " << a * a.inverse() << "\n\n";

    std::cout << "c == d is " << yes_or_no(c == d) << "\nc != d is " << yes_or_no(c != d);
    std::cout << "\ne == e is " << yes_or_no(e == e) << "\ne != e is " << yes_or_no(e != e) << "\n";

    std::cout << "\n\nquat.ij is: " << quat::ij() << "\nquat.jk is: " << quatD::jk()
              << "\nquat.ki is: " << quatD::ki() << "\n";
    assert(quatD::ij() == quatD::k());
    assert(quatD::jk() == quatD::i());
    assert(quatD::ki() == quatD::j());

    std::cout << "\nquat.ji is: " << quatD::ji() << "\nquat.kj is: " << quatD::kj()
              << "\nquat.ik is: " << quatD::ik() << "\nquat.ijk is: " << quatD::ijk() << "\n";
    assert(quatD::ji() == -quatD::k());
    assert(quatD::kj() == -quatD::i());
    assert(quatD::ik() == -quatD::j());

    std::cout << "\nquat.ii is: " << quatD::ii() << "\nquat.jj is: " << quatD::jj()
              << "\nquat.kk is: " << quatD::kk() << "\n";
    assert(quatD::ii()  == -1);
    assert(quatD::jj()  == -1);
    assert(quatD::kk()  == -1);
    assert(quatD::ijk() == -1);

    std::cout << "\nangle (deg) between c and d is: " << c.angle(d) << "\n";
    quatD c_minus_d = c - d;
    std::cout << "c_minus_d is: " << c_minus_d;
    matrix3dD rot_matrix = c_minus_d.rot_matrix();
    std::cout << "rot_matrix of c_minus_d is: " << c_minus_d.rot_matrix() << "\n";

  double rad2_2 = sqrt(2)/2.0;
  std::cout << "// -------------- LEVEL FLIGHT -------------------')\n";
    plane_rotation("levelFlight(E)", quatD(1),                     {  1,  0,  0,   0,  1,  0,   0,  0,  1 });
    plane_rotation("levelFlight(N)", quatD(rad2_2, 0, rad2_2,  0), {  0,  0,  1,   0,  1,  0,  -1,  0,  0 });
    plane_rotation("levelFlight(W)", quatD(0,      0,  1,      0), { -1,  0,  0,   0,  1,  0,   0,  0, -1 });
    plane_rotation("levelFlight(S)", quatD(rad2_2, 0, -rad2_2, 0), {  0,  0, -1,   0,  1,  0,   1,  0,  0} );
    std::cout << "LEVEL FLIGHT assertions passed ..............................................\n";
    std::cout << "// --------- end LEVEL FLIGHT ------------------------)\n";

    std::cout << "// -------------- STRAIGHT UP -------------------')\n";
    plane_rotation("straightUp(E)", quatD(rad2_2, 0, 0, rad2_2),   {  0, -1,  0,   1,  0,  0,   0,  0,  1 } );
    plane_rotation("straightUp(N)", quatD(0.5, 0.5, 0.5, 0.5),     {  0,  0,  1,   1,  0,  0,   0,  1,  0 } );
    plane_rotation("straightUp(W)", quatD(0, rad2_2, rad2_2, 0),   {  0,  1,  0,   1,  0,  0,   0,  0, -1 } );
    plane_rotation("straightUp(S)", quatD(0.5, -0.5, -0.5, 0.5),   {  0,  0, -1,   1,  0,  0,   0, -1,  0 } );
    std::cout << "STRAIGHT UP assertions passed..............................................\n";
    std::cout << "// --------- end STRAIGHT UP ------------------------)\n\n";


    std::cout << "// -------------- STRAIGHT DOWN ------------------')\n";
    plane_rotation("straightDown(E)", quatD(rad2_2, 0, 0, -rad2_2), {  0,  1,  0,  -1,  0,  0,   0,  0,  1 } );
    plane_rotation("straightDown(E)", quatD(0.5, -0.5, 0.5, -0.5),  {  0,  0,  1,  -1,  0,  0,   0, -1,  0 });
    plane_rotation("straightDown(E)", quatD(0, -rad2_2, rad2_2, 0), {  0, -1,  0,  -1,  0,  0,   0,  0, -1 } );
    plane_rotation("straightDown(E)", quatD(0.5, 0.5, -0.5, -0.5),  {  0,  0, -1,  -1,  0,  0,   0,  1,  0 });
    std::cout << "STRAIGHT DOWN assertions passed..............................................\n";
    std::cout << "// --------- end STRAIGHT DOWN ----------------------)\n\n";


    std::cout << "\n\n -------- BANK/ROLL ----------------\n";
    std::cout << "\nBanking/Rolling 90 degrees left...\n";
    plane_rotation("plane_E_bankLeft90", quatD(rad2_2, rad2_2, 0, 0),  {  1,  0,  0,   0,  0, -1,   0,  1,  0 } );
    plane_rotation("plane_N_bankLeft90", quatD(0.5, 0.5, 0.5, -0.5),   {  0,  1,  0,   0,  0, -1,  -1,  0,  0 } );
    plane_rotation("plane_W_bankLeft90", quatD(0, 0, rad2_2, -rad2_2), { -1,  0,  0,   0,  0, -1,   0, -1,  0 }  );
    plane_rotation("plane_W_bankLeft90", quatD(0.5, 0.5, -0.5, 0.5),   {  0, -1,  0,   0,  0, -1,   1,  0,  0 } );
    std::cout << "ROLL 90 deg left assertions passed..............................................\n";

    std::cout << "\n\nBanking/Rolling 180 degrees...\n";
    plane_rotation("plane_E_bankLeft180", quatD(0, 1, 0, 0),            {  1,  0,  0,   0, -1,  0,   0,  0, -1 });
    plane_rotation("plane_N_bankLeft180", quatD(0, rad2_2, 0, -rad2_2), {  0,  0, -1,   0, -1,  0,  -1,  0,  0 });
    plane_rotation("plane_W_bankLeft180", quatD(0, 0, 0, 1),            { -1,  0,  0,   0, -1,  0,   0,  0,  1 });
    plane_rotation("plane_S_bankLeft180", quatD(0, rad2_2, 0, rad2_2),  {  0,  0,  1,   0, -1,  0,   1,  0,  0 });
    std::cout << "ROLL 180 degrees assertions passed..............................................\n";

    std::cout << "\n\nBanking/Rolling 90 degrees right...\n";
    plane_rotation("plane_E_bankRight90", quatD(rad2_2, -rad2_2, 0, 0), {  1,  0,  0,   0,  0,  1,   0, -1,  0 });
    plane_rotation("plane_N_bankRight90", quatD(0.5, -0.5, 0.5, 0.5),   {  0, -1,  0,   0,  0,  1,  -1,  0,  0 });
    plane_rotation("plane_W_bankRight90", quatD(0, 0, rad2_2, rad2_2),  { -1,  0,  0,   0,  0,  1,   0,  1,  0 });
    plane_rotation("plane_S_bankRight90", quatD(0.5, -0.5, -0.5, -0.5), {  0,  1,  0,   0,  0,  1,   1,  0,  0 });
    std::cout << "ROLL 90 deg right assertions passed..............................................\n";
    std::cout << "\n -------- end BANK/ROLL ----------------\n";

    std::cout << "\nALL PLANE ROTATION ASSERTIONS PASSED ............................................\n\n";


    print("...test_quarternion assertions passed");
    print("====================  FINISHED testing quarternions  ========================");

}

int main(int argc, const char * argv[]) {
 //test_vectors();
 //test_matrices();
 //test_matrices_and_vectors();
 test_quarternion();
    
  return 0;
}


