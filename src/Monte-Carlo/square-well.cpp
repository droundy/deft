#include <stdlib.h>
#include "Monte-Carlo/square-well.h"
#include "handymath.h"

ball::ball() {
  pos = vector3d();
  R = 1;
  neighbors = new int[0];
  num_neighbors = 0;
  neighbor_center = vector3d();
}

ball::ball(const ball &p) {
  pos = p.pos;
  R = p.R;
  neighbors = p.neighbors;
  num_neighbors = p.num_neighbors;
  neighbor_center = p.neighbor_center;
}

ball ball::operator=(const ball &p) {
  pos = p.pos;
  R = p.R;
  neighbors = p.neighbors;
  num_neighbors = p.num_neighbors;
  neighbor_center = p.neighbor_center;
  return *this;
}

vector3d tmp_fix_periodic(vector3d v, const double len[3]) {
  for (int i = 0; i < 3; i++) {
    while (v[i] > len[i])
      v[i] -= len[i];
    while (v[i] < 0.0)
      v[i] += len[i];
  }
  return v;
}

vector3d tmp_periodic_diff(const vector3d &a, const vector3d  &b,
                       const double periodic[3]) {
  vector3d v = b - a;
  for (int i = 0; i < 3; i++) {
    if (periodic[i] > 0) {
      while (v[i] > periodic[i]/2.0)
        v[i] -= periodic[i];
      while (v[i] < -periodic[i]/2.0)
        v[i] += periodic[i];
    }
  }
  return v;
}

int initialize_neighbor_tables(ball *p, int N, double neighborR,
                               int max_neighbors,
                               const double periodic[3]) {
  int most_neighbors = 0;
  for (int i = 0; i < N; i++) {
    p[i].neighbor_center = p[i].pos;
  }
  for(int i = 0; i < N; i++) {
    p[i].neighbors = new int[max_neighbors];
    p[i].num_neighbors = 0;
    for (int j = 0; j < N; j++) {
      const bool is_neighbor = (i != j) &&
        (tmp_periodic_diff(p[i].pos, p[j].pos, periodic).normsquared() <
         sqr(p[i].R + p[j].R + neighborR));
      if (is_neighbor) {
        const int index = p[i].num_neighbors;
        p[i].num_neighbors++;
        if (p[i].num_neighbors > max_neighbors) return -1;
        p[i].neighbors[index] = j;
      }
    }
    most_neighbors = max(most_neighbors, p[i].num_neighbors);
  }
  return most_neighbors;
}

void update_neighbors(ball &a, int n, const ball *bs, int N,
                      double neighborR, const double periodic[3]) {
  a.num_neighbors = 0;
  for (int i = 0; i < N; i++) {
    if ((i != n) &&
        (tmp_periodic_diff(a.pos, bs[i].neighbor_center,
                       periodic).normsquared()
         < sqr(a.R + bs[i].R + neighborR))) {
      a.neighbors[a.num_neighbors] = i;
      a.num_neighbors++;
    }
  }
}

inline void add_neighbor(int new_n, ball *p, int id) {
  int i = p[id].num_neighbors;
  while (i > 0 && p[id].neighbors[i-1] > new_n) {
    p[id].neighbors[i] = p[id].neighbors[i-1];
    i --;
  }
  p[id].neighbors[i] = new_n;
  p[id].num_neighbors ++;
}

inline void remove_neighbor(int old_n, ball *p, int id) {
  int i = p[id].num_neighbors - 1;
  int temp = p[id].neighbors[i];
  while (temp != old_n) {
    i --;
    const int temp2 = temp;
    temp = p[id].neighbors[i];
    p[id].neighbors[i] = temp2;
  }
  p[id].num_neighbors --;
}

void inform_neighbors(const ball &new_p, const ball &old_p, ball *p, int n) {
  int new_index = 0, old_index = 0;
  while (true) {
    if (new_index == new_p.num_neighbors) {
      for(int i = old_index; i < old_p.num_neighbors; i++)
        remove_neighbor(n, p, old_p.neighbors[i]);
      return;
    }
    if (old_index == old_p.num_neighbors) {
      for(int i = new_index; i < new_p.num_neighbors; i++)
        add_neighbor(n, p, new_p.neighbors[i]);
      return;
    }
    if (new_p.neighbors[new_index] < old_p.neighbors[old_index]) {
      add_neighbor(n, p, new_p.neighbors[new_index]);
      new_index ++;
    } else if (old_p.neighbors[old_index] < new_p.neighbors[new_index]) {
      remove_neighbor(n, p, old_p.neighbors[old_index]);
      old_index ++;
    } else {
      new_index ++;
      old_index ++;
    }
  }
}

bool overlap(const ball &a, const ball &b, const double periodic[3], double dr) {
  const vector3d ab = tmp_periodic_diff(a.pos, b.pos, periodic);
  (ab.normsquared() > sqr(a.R + b.R + 2*dr)) ? false : true;
}


int overlaps_with_any(const ball &a, const ball *p,
                      const double periodic[3], double dr) {
  for (int i = 0; i < a.num_neighbors; i++) {
    if (overlap(a,p[a.neighbors[i]],periodic,dr))
      return true;
  }
  return false;
}

bool in_cell(const ball &p, const double walls[3], bool has_walls, double dr) {
  if(has_walls) {
    for (int i = 0; i < 3; i++) {
      if (walls[i] > 0) {
        if (p.pos[i]-p.R-dr < 0.0 || p.pos[i]+p.R+dr > walls[i])
          return false;
      }
    }
  }
  return true;
}

ball random_move(const ball &p, double size, const double len[3]) {
  ball temp = p;
  temp.pos = tmp_fix_periodic(temp.pos + vector3d::ran(size), len);
  return temp;
}

int move_one_ball(int id, ball *p, int N, const double periodic[3],
                  const double walls[3], bool real_walls, double neighborR,
                  double dist, int max_neighbors, double dr) {
  const double len[3] = {periodic[0]+walls[0], periodic[1]+walls[1], \
                         periodic[2]+walls[2]};
  ball temp = random_move(p[id], dist, len);
  int return_val = 0;
  if (in_cell(temp, walls, real_walls)) {
    bool overlaps = overlaps_with_any(temp, p, periodic);
    if (!overlaps) {
      const bool get_new_neighbors =
        (tmp_periodic_diff(temp.pos, temp.neighbor_center,
                           periodic).normsquared() > sqr(neighborR/2.0));
      if (get_new_neighbors) {
        // If we've moved too far, then the overlap test may have given a false
        // negative. So we'll find our new neighbors, and check against them.
        // If we still don't overlap, then we'll have to update the tables
        // of our neighbors that have changed.
        temp.neighbors = new int[max_neighbors];
        update_neighbors(temp, id, p, N, neighborR + 2*dr, periodic);
        return_val += 2;
        // However, for this check (and this check only), we don't need to
        // look at all of our neighbors, only our new ones.
        // fixme: do this!
        //int *new_neighbors = new int[max_neighbors];

        overlaps = overlaps_with_any(temp, p, periodic);
        if (!overlaps) {
          // Okay, we've checked twice, just like Santa Clause, so we're definitely
          // keeping this move and need to tell our neighbors where we are now.
          temp.neighbor_center = temp.pos;
          inform_neighbors(temp, p[id], p, id);
          return_val += 4;
          delete[] p[id].neighbors;
        }
        else delete[] temp.neighbors;
      }
      if (!overlaps) {
        p[id] = temp;
        return return_val + 1; //move successful
      }
    }
  }
  return return_val; //move unsucessful
}

int fcc_faces(int nx, int ny, int nz) {
  return 3*nx*ny*nz;
}

int fcc_inner_corner_spots(int nx, int ny, int nz) {
  return (nx-1)*(ny-1)*(nz-1);
}

int fcc_boundary_corners(int nx, int ny, int nz) {
  return nx*(ny-1)+ny*(nz-1)+nz*(nx-1)+1;
}

int fcc_total(int nx, int ny, int nz) {
  return fcc_faces(nx,ny,nz) \
    + fcc_inner_corner_spots(nx,ny,nz) \
    + fcc_boundary_corners(nx,ny,nz);
}
