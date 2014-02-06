#include "vector3d.h"
#pragma once

struct ball {
  vector3d pos;
  double R;
  int *neighbors;
  int num_neighbors;
  vector3d neighbor_center;

  ball();
  ball(const ball &p);

  ball operator=(const ball &p);
};

// Modulates v to within the periodic boundaries of the cell
vector3d sw_fix_periodic(vector3d v, const double len[3]);

// Return the vector pointing from a to b, accounting for periodic boundaries
vector3d periodic_diff(const vector3d &a, const vector3d  &b,
                          const double len[3], int num_walls);

// Create and initialize the neighbor tables for all balls (p).
// Returns the maximum number of neighbors that any ball has,
// or -1 if that number is larger than max_neighbors.
int initialize_neighbor_tables(ball *p, int N, double neighborR,
                               int max_neighbors, const double len[3],
                               int num_walls);

// Find's the neighbors of a by comparing a's position to the center of
// everyone else's neighborsphere, where id is the index of a in p.
void update_neighbors(ball &a, int id, const ball *p, int N,
                      double neighborR, const double len[3], int num_walls);

// Add ball new_n to the neighbor table of ball id
void add_neighbor(int new_n, ball *p, int id);

// Remove ball old_n from the neighbor table of ball id
void remove_neighbor(int old_n, ball *p, int id);

// Removes p from the neighbor table of anyone neighboring old_p.
// Adds p to the neighbor table of anyone neighboring new_p.
void inform_neighbors(const ball &new_p, const ball &old_p, ball *p);

// Check whether two balls overlap
bool overlap(const ball &a, const ball &b, const double len[3],
             int num_walls, double dr = 0);

// Check whether ball a overlaps with any of its neighbors in p.
// If count is true, it will return the total number of overlaps, otherwise constant
// it returns 1 if there is at least one overlap, 0 if there are none.
// If dr is nonzero, then each ball is treated as having a radius R + dr
int overlaps_with_any(const ball &a, const ball *p, const double len[3],
                      int num_walls, double dr=0);

// Return true if p doesn't intersect walls
bool in_cell(const ball &p, const double len[3], const int num_walls,
             double dr = 0);

// Move the ball by a random amount, in a gaussian distribution with
// respective standard deviations dist and angwidth
ball random_move(const ball &original, double size, const double len[3]);

// Attempt to move ball of id in p, while paying attention to collisions and walls
// Returns the sum of:
// 1 if the move was successful
// 2 if the neighbor table was updated
// 4 if, after updating the neighbor table, neighbors were informed
int move_one_ball(int id, ball *p, int N, double len[3],
                  int num_walls, double neighborR,
                  double dist, int max_neighbors, double dr);

// Count the number of interactions a given ball has
int count_interactions(int id, ball *p, double interaction_distance,
                       double len[3], int num_walls);
