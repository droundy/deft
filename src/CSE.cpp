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
  if (type != e.type) return false;
  if (name != e.name) return false;
  if (kind != e.kind) return false;
  if (arg3) return *arg1 == *e.arg1 && *arg2 == *e.arg2 && *arg3 == *e.arg3;
  if (arg2) return *arg1 == *e.arg1 && *arg2 == *e.arg2;
  if (arg1) return *arg1 == *e.arg1;
  return true;
}

bool Expression::EliminateThisSubexpression(const Expression &c, const std::string name) {
  if (c == *this) {
    *this = Expression(name);
    type = c.type;
    return true;
  } else {
    bool changed = false;
    if (arg1) changed = changed || arg1->EliminateThisSubexpression(c, name);
    if (arg2) changed = changed || arg2->EliminateThisSubexpression(c, name);
    if (arg3) changed = changed || arg3->EliminateThisSubexpression(c, name);
    return changed;
  }
}

Expression Expression::FindCommonSubexpression() const {
  Expression cs;
  if (arg1) {
    cs = arg1->FindCommonSubexpression();
    if (cs.unlazy) return cs;
  }
  if (arg2) {
    cs = arg2->FindCommonSubexpression();
    if (cs.unlazy) return cs;
  }
  if (arg3) {
    cs = arg3->FindCommonSubexpression();
    if (cs.unlazy) return cs;
  }
  return *this;
}
