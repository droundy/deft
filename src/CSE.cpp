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

#include "Expression.h"
#include <sstream>
#include <iostream>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>

bool Expression::operator==(const Expression &e) const {
  if (type != e.type) return false; // pointer comparison for speed!
  if (kind != e.kind) return false; // pointer comparison for speed!
  //if (!typeIs(e.type)) return false;
  //if (!kindIs(e.kind)) return false;
  if (name != e.name) return false;
  if (arg3) return *arg1 == *e.arg1 && *arg2 == *e.arg2 && *arg3 == *e.arg3;
  if (arg2) return *arg1 == *e.arg1 && *arg2 == *e.arg2;
  if (arg1) return *arg1 == *e.arg1;
  return true;
}

void Expression::EliminateThisDouble(const Expression &c, const std::string name) {
  if (c.typeIs(type)) {
    // If the type is right, then it's possible that we *are* the subexpression.
    // First check if we are the same as the subexpression we're trying to eliminate.
    if (c == *this) {
      Expression n(name);
      n.type = "double";
      n.alias = c.alias;
      *this = n;
      return;
    }
    // Now check if perhaps we're the same thing up to a scalar prefactor.
    //Expression tmp(*this);
    //Expression sca = tmp.ScalarFactor();
    //if (c == tmp) {
    //  *this = sca * n;
    //  return true;
    //}
  }
  // Try to recursively eliminate this subexpression in our children.
  if (arg1) arg1->EliminateThisDouble(c, name);
  if (arg2) arg2->EliminateThisDouble(c, name);
  if (arg3) arg3->EliminateThisDouble(c, name);
}

void Expression::EliminateThisSubexpression(const Expression &c, const std::string name) {
  if (c.typeIs(type)) {
    // If the type is right, then it's possible that we *are* the subexpression.
    // First check if we are the same as the subexpression we're trying to eliminate.
    if (c.kindIs(kind) && c == *this) {
      Expression n(name);
      n.type = c.type;
      n.alias = c.alias;
      *this = n;
      return;
    }
    // Now check if perhaps we're the same thing up to a scalar prefactor.
    //Expression tmp(*this);
    //Expression sca = tmp.ScalarFactor();
    //if (c == tmp) {
    //  *this = sca * n;
    //  return true;
    //}
  }
  // Try to recursively eliminate this subexpression in our children.
  if (arg1) arg1->EliminateThisSubexpression(c, name);
  if (arg2) arg2->EliminateThisSubexpression(c, name);
  if (arg3) arg3->EliminateThisSubexpression(c, name);
}

int Expression::CountThisSubexpression(const Expression &c) const {
  // First check if we are the same as the subexpression we're trying to eliminate.
  if (typeIs(c.type)) {
    // If the type is right, then it's possible that we *are* the subexpression.
    Expression n(name);
    n.type = c.type;
    // First check if we are the same as the subexpression we're trying to eliminate.
    if (c == *this) return 1;
    // Now check if perhaps we're the same thing up to a scalar prefactor.
    //Expression tmp(*this);
    //tmp.ScalarFactor();
    //if (c == tmp) return 1;
  }
  // Recursively count this subexpression in all our children.
  int num = 0;
  if (arg1) num += arg1->CountThisSubexpression(c);
  if (arg2) num += arg2->CountThisSubexpression(c);
  if (arg3) num += arg3->CountThisSubexpression(c);
  return num;
}

bool Expression::FindVariable(const std::string n) const {
  if (name == n && kindIs("variable")) {
    return true;
  } else {
    if (arg1 && arg1->FindVariable(n)) return true;
    if (arg2 && arg2->FindVariable(n)) return true;
    if (arg3 && arg3->FindVariable(n)) return true;
  }
  return false;
}

Expression Expression::FindCommonSubexpression() const {
  Expression cs;
  if (arg2) {
    int d1 = arg1->depth;
    int d2 = arg2->depth;
    if (d1 > d2) {
      cs = arg1->FindCommonSubexpression();
      if (cs.unlazy) return cs;
      cs = arg2->FindCommonSubexpression();
      if (cs.unlazy) return cs;
    } else {
      cs = arg2->FindCommonSubexpression();
      if (cs.unlazy) return cs;
      cs = arg1->FindCommonSubexpression();
      if (cs.unlazy) return cs;
    }
  } else if (arg1) {
    cs = arg1->FindCommonSubexpression();
    if (cs.unlazy) return cs;
  }
  if (arg3) {
    cs = arg3->FindCommonSubexpression();
    if (cs.unlazy) return cs;
  }
  return *this;
}

int Expression::Size() const {
  int depth = 1;
  if (arg1) depth += arg1->Size();
  if (arg2) depth += arg2->Size();
  if (arg3) depth += arg3->Size();
  return depth;
}

bool Expression::IsUnlazy() const {
  if (unlazy) return true;
  if (arg1 && arg1->IsUnlazy()) return true;
  if (arg2 && arg2->IsUnlazy()) return true;
  if (arg3 && arg3->IsUnlazy()) return true;
  return false;
}

bool Expression::IsUnlazyApartFrom(const Expression &c, std::set<std::string> important) const {
  if (important.count(alias)) return true; // This isn't technically unlazy...
  if (*this == c) return false;
  if (unlazy) return true;
  if (arg1 && arg1->IsUnlazyApartFrom(c, important)) return true;
  if (arg2 && arg2->IsUnlazyApartFrom(c, important)) return true;
  if (arg3 && arg3->IsUnlazyApartFrom(c, important)) return true;
  return false;
}

Expression Expression::EasyParentOfThisSubexpression(const Expression &e, std::set<std::string> important) const {
  const int ncount = CountThisSubexpression(e);
  if (ncount == 0 || e == *this) return e;
  if (!unlazy && (typeIs(e.type) || !e.kindIs("variable"))) {
    bool unlazychildren = false;
    if (arg1 && arg1->IsUnlazyApartFrom(e, important)) {
      //printf("My child %s is unlazy.\n", arg1->printme().c_str());
      unlazychildren = true;
    }
    if (arg2 && arg2->IsUnlazyApartFrom(e, important)) unlazychildren = true;
    if (arg3 && arg3->IsUnlazyApartFrom(e, important)) unlazychildren = true;
    if (!unlazychildren) {
      //printf("Found one!!! %s\n", printme().c_str());
      return *this;
    }
  }
  //printf("Looking in my children, says %s\n", printme().c_str());

  // If we ourselves aren't the easy parent, let's see if our children
  // have such a possibility...
  if (arg1) {
    Expression easyparent = arg1->EasyParentOfThisSubexpression(e, important);
    // The following is a bit too all-or-nothing: we ought to look for
    // a less "agressive" easy parent (which is a child of this
    // parent), but that'd be more work to code...
    if (easyparent != e && CountThisSubexpression(easyparent) == ncount) {
      return easyparent;
    }
    if (arg2) {
      easyparent = arg2->EasyParentOfThisSubexpression(e, important);
      if (easyparent != e && CountThisSubexpression(easyparent) == ncount) {
        return easyparent;
      }
      if (arg3) {
        easyparent = arg3->EasyParentOfThisSubexpression(e, important);
        if (easyparent != e && CountThisSubexpression(easyparent) == ncount) {
          return easyparent;
        }
      }
    }
  }
  return e;
}
