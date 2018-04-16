// xorshift.h

// This is the xorshift1024* random number generated, with the
// essential code taken from Wikipedia, with some code to initialize,
// etc.

#ifndef XORSHIFT_H
#define XORSHIFT_H

#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

class Rand {
public:
  uint64_t s[16];
  int p;
public:
	Rand( const uint64_t oneSeed ) {  // initialize with a simple uint64_t
    seed(oneSeed);
  }
	Rand() {  // auto-initialize with /dev/urandom or time() and clock()
    seed();
  }

	Rand( const Rand& o ) {  // copy constructor
    for (int i=0;i<16;i++) s[i] = o.s[i];
    p = o.p;
  }
	Rand& operator=( const Rand& o ) {
    for (int i=0;i<16;i++) s[i] = o.s[i];
    p = o.p;
    return *this;
  }

  void dump_resume_info(FILE *f) {
    fprintf(f, "# ");
    for (int i=0;i<16;i++) {
      fprintf(f, "%016lx ", s[i]);
    }
    fprintf(f, "%x\n", p);
  }
  void resume_from_dump(FILE *f) {
    if (fscanf(f, "# %lx", &s[0]) != 1) {
      fprintf(stderr, "Unable to read random resume information!\n");
      exit(1);
    }
    for (int i=1; i<16; i++) {
      if (fscanf(f, " %lx", &s[i]) != 1) {
        fprintf(stderr, "Unable to read random resume information!\n");
        exit(1);
      }
    }
    if (fscanf(f, " %x\n", &p) != 1) {
      fprintf(stderr, "Unable to read random resume information!\n");
      exit(1);
    }
  }

  uint64_t rand64() {
    uint64_t s0 = s[ p ];
    uint64_t s1 = s[ p = (p+1) & 15 ];
    s1 ^= s1 << 31; // a
    s1 ^= s1 >> 11; // b
    s0 ^= s0 >> 30; // c
    return ( s[p] = s0 ^ s1 ) * 1181783497276652981ULL;
  }
  double rand() {
    return ldexp(rand64(), -64);
  }

	// Re-seeding functions with same behavior as initializers
	void seed( uint64_t oneSeed ) {
    for (int i=0; i<16; i++) {
      oneSeed = (oneSeed << 3) + (oneSeed >> 61); // a terribly crude PRNG
      oneSeed ^= 5; // to avoid getting all zeros when given a seed of 0
      s[i] = oneSeed;
    }
    p = oneSeed & 15;
  }
	void seed() {
    // Seed the generator with an array from /dev/urandom if available
    // Otherwise use a combination of time(), pid, a pointer, and
    // clock() values
    FILE* urandom = fopen( "/dev/urandom", "rb" );
    if (urandom) {
      // we don't check the return value of fread, but instead just
      // xor with our low-quality random data regardless.
      fread(s, sizeof(uint64_t), 16, urandom );
      fclose(urandom);
    }
    // Now use time(), clock(), etc, to get some additional
    // low-quality randomness, mostly in case we can't open or read
    // from /dev/urandom.
    for (int i=0; i<16; i += 4) {
      s[(i+0)&15] ^= uint64_t(time(0));
      s[(i+1)&15] ^= uint64_t(this);
      s[(i+2)&15] ^= uint64_t(clock());
      s[(i+3)&15] ^= uint64_t(getpid());
      s[(i+4)&15] ^= uint64_t(getpid());
      s[(i+5)&15] ^= uint64_t(urandom);
      s[(i+5)&15] ^= uint64_t(rand());
    }
    p = clock() & 15;
  }
};

#endif  // XORSHIFT_H
