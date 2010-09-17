// Deft is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2010 The Deft Authors
//
// Deft is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with deft; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// Please see the file AUTHORS for a list of authors.

#include <stdio.h>
#include <string.h>
#include "Expression.h"

int retval = 0;

void test_expression(const Expression &e, const char *expected) {
  printf("\n**************");
  for (unsigned i=0;i<strlen(expected);i++) printf("*");
  printf("\n* Testing \"%s\" *\n", expected);
  for (unsigned i=0;i<strlen(expected);i++) printf("*");
  printf("**************\n");
  int retval = 0;

  if (e.printme() != expected) {
    printf("FAIL: Got instead: \"%s\"\n", e.printme().c_str());
    retval++;
  }
}


int main(int, char **argv) {
  Expression foobar("foobar"), foo("foo"), bar("bar");
  test_expression(foobar, "foobar");
  test_expression(0, "0");
  test_expression(1 + foobar, "1 + foobar");
  test_expression(2 * foobar, "2*foobar");
  test_expression(2 / foobar, "2/foobar");
  test_expression(1 * foobar, "foobar");
  test_expression(2 * Expression(3), "6");
  test_expression(2 + Expression(3), "5");
  test_expression(0 * foobar, "0");

  // Test that functions work
  test_expression(funexpr("f", 3), "f(3)");
  test_expression(funexpr("f", foobar + 3), "f(foobar + 3)");

  // Test that methods work
  test_expression(foobar.method("foo"), "foobar.foo()");
  test_expression(foobar.method("foo", 3), "foobar.foo(3)");
  test_expression(foobar.method("foo", bar, bar + foo), "foobar.foo(bar, bar + foo)");
  test_expression((foo + bar).method("foo", bar, bar + foo), "(foo + bar).foo(bar, bar + foo)");

  // Test associativity...
  test_expression(1 + foo * bar, "1 + foo*bar");
  test_expression((1 + foo) * bar, "(1 + foo)*bar");
  test_expression((bar + 1) + foo, "bar + 1 + foo");
  test_expression(bar + (1 + foo), "bar + (1 + foo)");
  test_expression(1 + 1 + foo, "2 + foo");
  test_expression(1 * (1 * foo), "foo");
  test_expression(foobar * (bar * foo), "foobar*(bar*foo)");

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
