// Dafty is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2010 The Dafty Authors
//
// Dafty is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with dafty; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// Please see the file AUTHORS for a list of authors.

#include <stdio.h>
#include <eigen2/Eigen/Core>

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class Foo : public VectorXd {
  int xxx;
public:
  explicit Foo(int num) : VectorXd(num) { xxx = 4; }
};

int main() {
  printf("Hello world!\n");
  VectorXd v(7);
  v << 1,2,3,4,5,6,7;
  std::cout << v << std::endl;
  Foo myfoo(7);
  // The following fails...
  // myfoo = v;
  std::cout << "myfoo is " << myfoo << std::endl;  
  return 0;
}
