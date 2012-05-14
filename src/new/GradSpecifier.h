// -*- mode: C++; -*-

#pragma once

#include <cassert>

enum Inclusion { skip, include };

// The GradSpecifier type allows us to specify which parts of a Vector (that is conceptually made of parts) are included in a derivative.

class GradSpecifier {
public:
  GradSpecifier() : inc(skip), next(0) {}
  GradSpecifier(const GradSpecifier &a) : inc(a.inc), next(0) {
    if (a.next) next = new GradSpecifier(*a.next);
  }
  void operator=(const GradSpecifier &a) {
    delete next;
    inc = a.inc;
    if (a.next) next = new GradSpecifier(*a.next);
    else next = 0;
  }

  // operator[] is a pair of "safe" operator for indexing a vector
  Inclusion operator[](int i) const {
    if (i <= 0) return inc;
    if (!next) return include; // We include everything beyond our bounds!
    return (*next)[i-1];
  };
  Inclusion &operator[](int i) {
    if (i <= 0) return inc;
    assert(next);
    return (*next)[i-1];
  };
private:
  Inclusion inc;
  GradSpecifier *next;
};
