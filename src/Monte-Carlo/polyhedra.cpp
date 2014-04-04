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

bool overlap(const polyhedron &a, const polyhedron &b, const double periodic[3], double dr) {
  const vector3d ab = periodic_diff(a.pos, b.pos, periodic);
  if (ab.normsquared() > sqr(a.R + b.R + 2*dr))
    return false;
  // construct axes from a
  // project a and b to each axis
  for (int i=0; i<a.mypoly->nfaces; i++) {
    const vector3d axis = a.rot.rotate_vector(a.mypoly->faces[i]);
    double projection = axis.dot(a.rot.rotate_vector((a.mypoly->vertices[0])*(a.R+dr)));
    double mina = projection, maxa = projection;
    for (int j=1; j<a.mypoly->nvertices; j++) {
      projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[j]*(a.R+dr)));
      if (projection < mina) mina = projection;
      else if (projection > maxa) maxa = projection;
    }
    projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[0]*(b.R+dr)) + ab);
    double minb = projection, maxb = projection;
    for (int j=1; j<b.mypoly->nvertices; j++) {
      projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[j]*(b.R+dr)) + ab);
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
    double projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[0]*(a.R+dr)));
    double mina = projection, maxa = projection;
    for (int j=1; j<a.mypoly->nvertices; j++) {
      projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[j]*(a.R+dr)));
      if (projection < mina) mina = projection;
      else if (projection > maxa) maxa = projection;
    }
    projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[0]*(b.R+dr)) + ab);
    double minb = projection, maxb = projection;
    for (int j=1; j<b.mypoly->nvertices; j++) {
      projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[j]*(b.R+dr)) + ab);
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
    double projection = aaxes[i].dot(a.rot.rotate_vector(a.mypoly->vertices[0]*(a.R+dr)));
    amins[i] = projection, amaxes[i] = projection;
    for (int j=1; j< a.mypoly->nvertices; j++) {
      projection = aaxes[i].dot(a.rot.rotate_vector(a.mypoly->vertices[j]*(a.R+dr)));
      amins[i] = min(projection, amins[i]);
      amaxes[i] = max(projection, amaxes[i]);
    }
  }
  int num_overlaps = 0;
  for (int l=0; l<a.num_neighbors; l++) {
    const int k = a.neighbors[l];
    const vector3d ab = periodic_diff(a.pos, bs[k].pos, periodic);
    if (ab.normsquared() < sqr(a.R + bs[k].R + 2*dr)) {
      bool overlap = true; // assume overlap until we prove otherwise or fail to.
      // check projection of b against a's axes
      for (int i=0; i<a.mypoly->nfaces; i++) {
        double projection = aaxes[i].dot
          (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[0]*(bs[k].R+dr)) + ab);
        double bmin = projection, bmax = projection;
        for (int j=1; j<bs[k].mypoly->nvertices; j++) {
          projection = aaxes[i].dot
            (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[j]*(bs[k].R+dr)) + ab);
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
          double projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[0]*(a.R+dr)));
          double amin = projection, amax = projection;
          for (int j=1; j<a.mypoly->nvertices; j++) {
            projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[j]*(a.R+dr)));
            amin = min(projection, amin);
            amax = max(projection, amax);
          }
          projection = axis.dot
            (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[0]*(bs[k].R+dr)) + ab);
          double bmin = projection, bmax = projection;
          for (int j=1; j<bs[k].mypoly->nvertices; j++) {
            projection = axis.dot
              (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[j]*(bs[k].R+dr)) + ab);
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

bool in_cell(const polyhedron &p, const double walls[3], bool real_walls, double dr) {
  if(real_walls) {
    for (int i=0; i<3; i++) {
      if (walls[i] > 0) {
        if (p.pos[i] - p.R - dr > 0.0 && p.pos[i] + p.R + dr < walls[i]) {
          continue;
        }
        double coord = (p.rot.rotate_vector(p.mypoly->vertices[0]*(p.R+dr)) + p.pos)[i];
        double pmin = coord, pmax = coord;
        for (int j=1; j<p.mypoly->nvertices; j++) {
          coord = (p.rot.rotate_vector(p.mypoly->vertices[j]*(p.R+dr)) + p.pos)[i];
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

counter move_one_polyhedron(int id, polyhedron *p, int N, const double periodic[3],
                         const double walls[3], bool real_walls, double neighborR,
                        double dist, double angwidth, int max_neighbors, double dr) {
  const double len[3] = {periodic[0]+walls[0], periodic[1]+walls[1], periodic[2]+walls[2]};
  polyhedron temp = random_move(p[id], dist, angwidth, len);
  counter move;
  move.totalmoves ++;
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
        update_neighbors(temp, id, p, N, neighborR + 2*dr, periodic);
        move.updates ++;
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
          move.informs ++;
          delete[] p[id].neighbors;
        }
        else delete[] temp.neighbors;
      }
      if (!overlaps) {
        p[id] = temp;
        move.workingmoves ++; // move sucessful
        return move;
      }
    }
  }
  return move; // move unsucessful
}

int initialize_neighbor_tables(polyhedron *p, int N, double neighborR,
                               int max_neighbors, const double periodic[3]) {
  int most_neighbors = 0;
  for(int i=0; i<N; i++) {
    p[i].neighbor_center = p[i].pos;
  }
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
  type = NONE;
  nvertices = 0;
  nfaces = 0;
  volume = 0;
  vertices = NULL;
  faces = NULL;
  name = new char[6];
  sprintf(name, "empty");
}

poly_shape::poly_shape(const char *set_name, double ratio) {
  // vertices and edges should be in an order such that the second
  // half are a reflection of the first half for polyhedra with
  // inversion symmetry (i.e. not tetrahedra)
  //
  // Also, the first two vertices should be along the same edge, so
  // that edge length can be easily calculated
  //
  // All face vectors should be normalized
  //
  // For tetrahedra, the faces should all be oriented so as to point
  // outward And for truncated tetrahedra, outward from the large face
  //
  // This is for constructing distribution functions nicely
  if (strcmp(set_name, "cube") == 0 || strcmp(set_name, "cuboid") == 0) {
    if (strcmp(set_name, "cube") == 0) {
      type = CUBE;
      ratio = 1;
      name = new char[5];
      sprintf(name, "cube");
    }
    else {
      type = CUBOID;
      name = new char[13];
      sprintf(name, "cuboid_%05.2f", ratio);
    }

    const double x = 1/sqrt(2 + sqr(ratio));
    const double y = x;
    const double z = ratio*x;

    nvertices = 8;
    nfaces = 3;
    nedges = 12;
    volume = 8.0*x*y*z;
    vertices = new vector3d[nvertices];
    faces = new vector3d[nfaces];
    edges = new vector3d[nedges];

    vertices[0] = vector3d( x,  y,  z);
    vertices[1] = vector3d(-x,  y,  z);
    vertices[2] = vector3d( x, -y,  z);
    vertices[3] = vector3d( x,  y, -z);

    vertices[4] = vector3d(-x, -y,  z);
    vertices[5] = vector3d(-x,  y, -z);
    vertices[6] = vector3d(-x, -y, -z);
    vertices[7] = vector3d( x, -y, -z);

    faces[0] = vector3d(1, 0, 0).normalized();
    faces[1] = vector3d(0, 1, 0).normalized();
    faces[2] = vector3d(0, 0, 1).normalized();


    edges[0]  = vector3d( 1,  1,  0).normalized();
    edges[1]  = vector3d( 1, -1,  0).normalized();
    edges[2]  = vector3d( 1,  0,  1).normalized();
    edges[3]  = vector3d( 1,  0, -1).normalized();
    edges[4]  = vector3d( 0,  1,  1).normalized();
    edges[5]  = vector3d( 0,  1, -1).normalized();

    edges[6]  = vector3d(-1,  1,  0).normalized();
    edges[7]  = vector3d(-1, -1,  0).normalized();
    edges[8]  = vector3d(-1,  0,  1).normalized();
    edges[9]  = vector3d(-1,  0, -1).normalized();
    edges[10] = vector3d( 0, -1,  1).normalized();
    edges[11] = vector3d( 0, -1, -1).normalized();
  }
  else if (strcmp(set_name, "tetrahedron") == 0) {
    type = TETRAHEDRON;
    nvertices = 4;
    nfaces = 4;
    nedges = 6;
    volume = 8.0/9.0/sqrt(3.0);
    vertices = new vector3d[nvertices];
    faces = new vector3d[nfaces];
    edges = new vector3d[nedges];
    name = new char[13];
    sprintf(name, "tetrahedron");

    vertices[0] = vector3d( sqrt(2.0/3.0),    -sqrt(2.0)/3.0, -1.0/3.0);
    vertices[1] = vector3d(-sqrt(2.0/3.0),    -sqrt(2.0)/3.0, -1.0/3.0);
    vertices[2] = vector3d(             0, 2.0*sqrt(2.0)/3.0, -1.0/3.0);
    vertices[3] = vector3d(             0,                 0,      1.0);

    faces[0] = vector3d(         0, -sqrt(2.0),           0.5).normalized();
    faces[1] = vector3d( sqrt(3.0),        1.0, 1.0/sqrt(2.0)).normalized();
    faces[2] = vector3d(-sqrt(3.0),        1.0, 1.0/sqrt(2.0)).normalized();
    faces[3] = vector3d(         0,          0,            -1).normalized();
  }
  else if (strcmp(set_name, "truncated_tetrahedron") == 0) {
    type = TRUNCATED_TETRAHEDRON;
    nvertices = 12;
    nfaces = 4;
    volume = 184.0/33.0/sqrt(11.0);
    vertices = new vector3d[nvertices];
    faces = new vector3d[nfaces];
    name = new char[22];
    sprintf(name, "truncated_tetrahedron");

    vertices[0]  = vector3d( 3,  1,  1);
    vertices[1]  = vector3d( 1,  3,  1);
    vertices[2]  = vector3d( 1,  1,  3);

    vertices[3]  = vector3d(-3, -1,  1);
    vertices[4]  = vector3d(-1, -3,  1);
    vertices[5]  = vector3d(-1, -1,  3);

    vertices[6]  = vector3d(-3,  1, -1);
    vertices[7]  = vector3d(-1,  3, -1);
    vertices[8]  = vector3d(-1,  1, -3);

    vertices[9]  = vector3d( 3, -1, -1);
    vertices[10] = vector3d( 1, -3, -1);
    vertices[11] = vector3d( 1, -1, -3);

    for(int i=0; i<nvertices; i++) vertices[i]/=sqrt(11.0);

    faces[0] = vector3d(-1, -1, -1).normalized();
    faces[1] = vector3d( 1,  1, -1).normalized();
    faces[2] = vector3d( 1, -1,  1).normalized();
    faces[3] = vector3d(-1,  1,  1).normalized();
  }
  else {
    type=NONE;
    nvertices = 0;
    nfaces = 0;
    nedges = 0;
    volume = 0;
    vertices = NULL;
    faces = NULL;
    edges = NULL;
    name = new char[14];
    sprintf(name, "invalid shape");
  }
  // Define edges:
  // fixme: nfaces isn't actual number of faces so can't use this
  nedges = nfaces + nvertices - 2;
  edges = new vector3d[nedges];

  
}

poly_shape::~poly_shape() {
  delete[] vertices;
  delete[] faces;
  delete[] edges;
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
