#include <stdlib.h>
#include "Monte-Carlo/polyhedra.h"
#include "handymath.h"

const poly_shape empty_shape;


void update_neighbors(polyhedron &a, int n, const polyhedron *bs, int N,
                      double neighborR, const double periodic[3]) {
  a.num_neighbors = 0;
  for(int i=0; i<N; i++) {
    if ((i!=n) &&
        (periodic_diff(a.pos, bs[i].neighbor_center, periodic).normsquared()
         < sqr(a.R + bs[i].R + neighborR))) {
      a.neighbors[a.num_neighbors] = i;
      a.num_neighbors ++;
    }
  }
}

inline void add_new_neighbor(int new_n, polyhedron *p, int id) {
  int pindex = p[id].num_neighbors - 1;
  p[id].num_neighbors ++;
  while (pindex >= 0 && p[id].neighbors[pindex] > new_n) {
    p[id].neighbors[pindex + 1] = p[id].neighbors[pindex];
    pindex --;
  }
  p[id].neighbors[pindex+1] = new_n;
}

inline void remove_old_neighbor(int old_n, polyhedron *p, int id) {
  int pindex = p[id].num_neighbors - 1;
  int temp = p[id].neighbors[pindex];
  while (temp != old_n) {
    pindex --;
    const int temp2 = temp;
    temp = p[id].neighbors[pindex];
    p[id].neighbors[pindex] = temp2;
  }
  p[id].num_neighbors --;
}

void inform_neighbors(const polyhedron &newp, const polyhedron &oldp, int n,
                      polyhedron *p) {
  int new_index = 0, old_index = 0;
  while (true) {
    if (new_index >= newp.num_neighbors) {
      for(int i=old_index; i<oldp.num_neighbors; i++)
        remove_old_neighbor(n, p, oldp.neighbors[i]);
      return;
    }
    if (old_index >= oldp.num_neighbors) {
      for(int i=new_index; i<newp.num_neighbors; i++)
        add_new_neighbor(n, p, newp.neighbors[i]);
      return;
    }
    if (newp.neighbors[new_index] < oldp.neighbors[old_index]) {
      add_new_neighbor(n, p, newp.neighbors[new_index]);
      new_index ++;
    } else if (oldp.neighbors[old_index] < newp.neighbors[new_index]) {
      remove_old_neighbor(n, p, oldp.neighbors[old_index]);
      old_index ++;
    } else { // (newp.neighbors[new_index] == oldp.neighbors[old_index])
      new_index ++;
      old_index ++;
    }
  }
}

vector3d fix_periodic(vector3d v, const double len[3]) {
  for (int i=0; i<3; i++) {
    while (v[i] > len[i])
      v[i] -= len[i];
    while (v[i] < 0.0)
      v[i] += len[i];
  }
  return v;
}

vector3d periodic_diff(const vector3d &a, const vector3d  &b, const double periodic[3]) {
  vector3d v = b - a;
  for (int i=0; i<3; i++) {
    if (periodic[i] > 0) {
      while (v[i] > periodic[i]/2.0)
        v[i] -= periodic[i];
      while (v[i] < -periodic[i]/2.0)
        v[i] += periodic[i];
    }
  }
  return v;
}

bool overlap(const polyhedron &a, const polyhedron &b, const double periodic[3]) {
  const vector3d ab = periodic_diff(a.pos, b.pos, periodic);
  if (ab.normsquared() > (a.R + b.R)*(a.R + b.R))
    return false;
  // construct axes from a
  // project a and b to each axis
  for (int i=0; i<a.mypoly->nfaces; i++) {
    const vector3d axis = a.rot.rotate_vector(a.mypoly->faces[i]);
    double projection = axis.dot(a.rot.rotate_vector((a.mypoly->vertices[0])*a.R));
    double mina = projection, maxa = projection;
    for (int j=1; j<a.mypoly->nvertices; j++) {
      projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[j]*a.R));
      if (projection < mina) mina = projection;
      else if (projection > maxa) maxa = projection;
    }
    projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[0]*b.R) + ab);
    double minb = projection, maxb = projection;
    for (int j=1; j<b.mypoly->nvertices; j++) {
      projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[j]*b.R) + ab);
      if (projection < minb) minb = projection;
      else if (projection > maxb) maxb = projection;
    }
    if (mina > maxb || minb > maxa) {
      return false;
    }
  }
  // construct axes from b
  // project a and b to each axis
  for (int i=0; i<b.mypoly->nfaces; i++) {
    const vector3d axis = b.rot.rotate_vector(b.mypoly->faces[i]);
    double projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[0]*a.R));
    double mina = projection, maxa = projection;
    for (int j=1; j<a.mypoly->nvertices; j++) {
      projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[j]*a.R));
      if (projection < mina) mina = projection;
      else if (projection > maxa) maxa = projection;
    }
    projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[0]*b.R) + ab);
    double minb = projection, maxb = projection;
    for (int j=1; j<b.mypoly->nvertices; j++) {
      projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[j]*b.R) + ab);
      if (projection < minb) minb = projection;
      else if (projection > maxb) maxb = projection;
    }
    if (mina > maxb || minb > maxa) {
      return false;
    }
  }
  return true;
}


int overlaps_with_any(const polyhedron &a, const polyhedron *bs,
                      const double periodic[3], bool count, double dr) {
  // construct axes from a and a's projection onto them
  vector3d aaxes[a.mypoly->nfaces];
  double amins[a.mypoly->nfaces], amaxes[a.mypoly->nfaces];
  for (int i=0; i<a.mypoly->nfaces; i++) {
    aaxes[i] = a.rot.rotate_vector(a.mypoly->faces[i]);
    double projection = aaxes[i].dot(a.rot.rotate_vector(a.mypoly->vertices[0]*a.R));
    amins[i] = projection, amaxes[i] = projection;
    for (int j=1; j< a.mypoly->nvertices; j++) {
      projection = aaxes[i].dot(a.rot.rotate_vector(a.mypoly->vertices[j]*a.R));
      amins[i] = min(projection, amins[i]);
      amaxes[i] = max(projection, amaxes[i]);
    }
  }
  int num_overlaps = 0;
  for (int l=0; l<a.num_neighbors; l++) {
    const int k = a.neighbors[l];
    const vector3d ab = periodic_diff(a.pos, bs[k].pos, periodic);
    if (ab.normsquared() < (a.R + bs[k].R)*(a.R + bs[k].R)) {
      bool overlap = true; // assume overlap until we prove otherwise or fail to.
      // check projection of b against a's axes
      for (int i=0; i<a.mypoly->nfaces; i++) {
        double projection = aaxes[i].dot
          (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[0]*bs[k].R) + ab);
        double bmin = projection, bmax = projection;
        for (int j=1; j<bs[k].mypoly->nvertices; j++) {
          projection = aaxes[i].dot
            (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[j]*bs[k].R) + ab);
          bmin = min(projection, bmin);
          bmax = max(projection, bmax);
        }
        if (amins[i] > bmax || bmin > amaxes[i]) {
          overlap = false;
          i = a.mypoly->nfaces; // no overlap, move on to next
        }
      }
      if (overlap) { // still need to check against b's axes
        for (int i=0; i<bs[k].mypoly->nfaces; i++) {
          const vector3d axis = bs[k].rot.rotate_vector(bs[k].mypoly->faces[i]);
          double projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[0]*a.R));
          double amin = projection, amax = projection;
          for (int j=1; j<a.mypoly->nvertices; j++) {
            projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[j]*a.R));
            amin = min(projection, amin);
            amax = max(projection, amax);
          }
          projection = axis.dot
            (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[0]*bs[k].R) + ab);
          double bmin = projection, bmax = projection;
          for (int j=1; j<bs[k].mypoly->nvertices; j++) {
            projection = axis.dot
              (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[j]*bs[k].R) + ab);
            bmin = min(projection, bmin);
            bmax = max(projection, bmax);

          }
          if (amin > bmax || bmin > amax) {
            overlap = false;
            i = bs[k].mypoly->nfaces; //no overlap, move on to next
          }
        }
      }
      if (overlap) {
        if(!count)
          return 1;
        num_overlaps ++;
      }
    }
  }
  return num_overlaps;
}

bool in_cell(const polyhedron &p, const double walls[3], bool real_walls) {
  if(real_walls) {
    for (int i=0; i<3; i++) {
      if (walls[i] > 0) {
        if (p.pos[i] - p.R > 0.0 && p.pos[i] + p.R < walls[i]) {
          continue;
        }
        double coord = (p.rot.rotate_vector(p.mypoly->vertices[0]*p.R) + p.pos)[i];
        double pmin = coord, pmax = coord;
        for (int j=1; j<p.mypoly->nvertices; j++) {
          coord = (p.rot.rotate_vector(p.mypoly->vertices[j]*p.R) + p.pos)[i];
          pmin = min(coord, pmin);
          pmax = max(coord, pmax);
        }
        if (pmin < 0.0 || pmax > walls[i])
          return false;
      }
    }
  }
  return true;
}


polyhedron random_move(const polyhedron &original, double size,
                                    double angwidth, const double len[3]) {
  polyhedron temp = original;
  temp.pos = fix_periodic(temp.pos + vector3d::ran(size), len);
  temp.rot = rotation::ran(angwidth)*original.rot;
  return temp;
}

bool move_one_polyhedron(int id, polyhedron *p, int N, const double periodic[3],
                         const double walls[3], bool real_walls, double neighborR,
                         double dist, double angwidth, int max_neighbors) {
  const double len[3] = {periodic[0]+walls[0], periodic[1]+walls[1], periodic[2]+walls[2]};
  polyhedron temp = random_move(p[id], dist, angwidth, len);
  if (in_cell(temp, walls, real_walls)) {
    bool overlaps = overlaps_with_any(temp, p, periodic);
    if (!overlaps) {
      const bool get_new_neighbors =
        (periodic_diff(temp.pos, temp.neighbor_center, periodic).normsquared() >
         sqr(neighborR/2.0));
      if (get_new_neighbors) {
        // If we've moved too far, then the overlap test may have given a false
        // negative. So we'll find our new neighbors, and check against them.
        // If we still don't overlap, then we'll have to update the tables
        // of our neighbors that have changed.
        temp.neighbors = new int[max_neighbors];
        update_neighbors(temp, id, p, N, neighborR, periodic);
        // However, for this check (and this check only), we don't need to
        // look at all of our neighbors, only our new ones.
        // fixme: do this!
        //int *new_neighbors = new int[max_neighbors];

        overlaps = overlaps_with_any(temp, p, periodic);
        if (!overlaps) {
          // Okay, we've checked twice, just like Santa Clause, so we're definitely
          // keeping this move and need to tell our neighbors where we are now.
          temp.neighbor_center = temp.pos;
          inform_neighbors(temp, p[id], id, p);
          delete[] p[id].neighbors;
        }
        else delete[] temp.neighbors;
      }
      if (!overlaps) {
        p[id] = temp;
        return true;
      }
    }
  }
  return false;
}

int initialize_neighbor_tables(polyhedron *p, int N, double neighborR,
                               int max_neighbors, const double periodic[3]) {
  int most_neighbors = 0;
  for(int i=0; i<N; i++) {
    p[i].neighbors = new int[max_neighbors];
    p[i].num_neighbors = 0;
    for(int j=0; j<N; j++) {
      const bool is_neighbor = (i != j) &&
        (periodic_diff(p[i].pos, p[j].pos, periodic).normsquared() <
         uipow(p[i].R + p[j].R + neighborR, 2));
      if (is_neighbor) {
        const int index = p[i].num_neighbors;
        p[i].num_neighbors ++;
        if (p[i].num_neighbors > max_neighbors) return -1;
        p[i].neighbors[index] = j;
      }
    }
    most_neighbors = max(most_neighbors, p[i].num_neighbors);
  }
  return most_neighbors;
}

poly_shape::poly_shape() {
  nvertices = 0;
  nfaces = 0;
  volume = 0;
  vertices = NULL;
  faces = NULL;
  name = new char[6];
  sprintf(name, "empty");
}

poly_shape::poly_shape(const char *set_name) {
  if (strcmp(set_name, "cube") == 0) {
    nvertices = 8;
    nfaces = 3;
    volume = 1.0;
    vertices = new vector3d[nvertices];
    faces = new vector3d[nfaces];
    name = new char[5];
    sprintf(name, "cube");

    const double v_cube = 1.0/sqrt(3.0);
    vertices[0] = vector3d( v_cube,  v_cube,  v_cube);
    vertices[1] = vector3d(-v_cube,  v_cube,  v_cube);
    vertices[2] = vector3d(-v_cube, -v_cube,  v_cube);
    vertices[3] = vector3d( v_cube, -v_cube,  v_cube);
    vertices[4] = vector3d( v_cube,  v_cube, -v_cube);
    vertices[5] = vector3d(-v_cube,  v_cube, -v_cube);
    vertices[6] = vector3d(-v_cube, -v_cube, -v_cube);
    vertices[7] = vector3d( v_cube, -v_cube, -v_cube);

    faces[0] = vector3d(1, 0, 0);
    faces[1] = vector3d(0, 1, 0);
    faces[2] = vector3d(0, 0, 1);
  }
  else if (strcmp(set_name, "tetrahedron") == 0) {
    nvertices = 4;
    nfaces = 4;
    volume = 1.0/3.0;
    vertices = new vector3d[nvertices];
    faces = new vector3d[nfaces];
    name = new char[13];
    sprintf(name, "tetrahedron");

    vertices[0] = vector3d( sqrt(2.0/3.0),    -sqrt(2.0)/3.0, -1.0/3.0);
    vertices[1] = vector3d(-sqrt(2.0/3.0),    -sqrt(2.0)/3.0, -1.0/3.0);
    vertices[2] = vector3d(             0, 2.0*sqrt(2.0)/3.0, -1.0/3.0);
    vertices[3] = vector3d(             0,                 0,      1.0);

    faces[0] = vector3d(         0, -sqrt(2.0),           0.5);
    faces[1] = vector3d( sqrt(3.0),        1.0, 1.0/sqrt(2.0));
    faces[2] = vector3d(-sqrt(3.0),        1.0, 1.0/sqrt(2.0));
    faces[3] = vector3d(         0,          0,             1);
  }
  else if (strcmp(set_name, "truncated_tetrahedron") == 0) {
    nvertices = 18;
    nfaces = 4;
    volume = 0; //fixme
    vertices = new vector3d[nvertices];
    faces = new vector3d[nfaces];
    name = new char[22];
    sprintf(name, "truncated tetrahedron");

    faces[0] = vector3d(         0, -sqrt(2.0),           0.5);
    faces[1] = vector3d( sqrt(3.0),        1.0, 1.0/sqrt(2.0));
    faces[2] = vector3d(-sqrt(3.0),        1.0, 1.0/sqrt(2.0));
    faces[3] = vector3d(         0,          0,             1);
  }
  else {
    nvertices = 0;
    nfaces = 0;
    volume = 0;
    vertices = NULL;
    faces = NULL;
    name = new char[14];
    sprintf(name, "invalid shape");
  }
}

poly_shape::~poly_shape() {
  delete[] vertices;
  delete[] faces;
  delete[] name;
}

polyhedron::polyhedron() {
  pos = vector3d();
  rot = rotation();
  R = 0;
  mypoly = &empty_shape;
  neighbors = new int[0];
  num_neighbors = 0;
  neighbor_center = vector3d();
}

polyhedron::polyhedron(const polyhedron &p) {
  pos = p.pos;
  rot = p.rot;
  R = p.R;
  mypoly = p.mypoly;
  neighbors = p.neighbors;
  num_neighbors = p.num_neighbors;
  neighbor_center = p.neighbor_center;
}

polyhedron polyhedron::operator=(const polyhedron &p) {
  pos = p.pos;
  rot = p.rot;
  R = p.R;
  mypoly = p.mypoly;
  neighbors = p.neighbors;
  num_neighbors = p.num_neighbors;
  neighbor_center = p.neighbor_center;
  return *this;
}

// polyhedron::~polyhedron() {
//   delete[] neighbors;
// }
