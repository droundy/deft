// -*- mode: C++; -*-

#pragma once

enum Verbosity { silent, quiet, verbose, chatty, gossipy, argumentative };


inline Verbosity quieter(Verbosity v) {
  return Verbosity(v-1);
}

inline Verbosity quietest(Verbosity v) {
  return Verbosity(v-2);
}

inline Verbosity louder(Verbosity v) {
  return Verbosity(v+1);
}

inline Verbosity loudest(Verbosity v) {
  return Verbosity(v+2);
}
