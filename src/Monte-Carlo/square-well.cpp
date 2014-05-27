#include <stdlib.h>
#include "Monte-Carlo/square-well.h"
#include "handymath.h"

ball::ball(){
  pos = vector3d();
  R = 1;
  neighbors = new int[0];
  num_neighbors = 0;
  neighbor_center = vector3d();
}

ball::ball(const ball &p){
  pos = p.pos;
  R = p.R;
  neighbors = p.neighbors;
  num_neighbors = p.num_neighbors;
  neighbor_center = p.neighbor_center;
}

ball ball::operator=(const ball &p){
  pos = p.pos;
  R = p.R;
  neighbors = p.neighbors;
  num_neighbors = p.num_neighbors;
  neighbor_center = p.neighbor_center;
  return *this;
}

move_info::move_info(){
  total_old = 0;
  total = 0;
  working = 0;
  working_old = 0;
  updates = 0;
  informs = 0;
}

vector3d sw_fix_periodic(vector3d v, const double len[3]){
  for (int i = 0; i < 3; i++){
    while (v[i] > len[i])
      v[i] -= len[i];
    while (v[i] < 0.0)
      v[i] += len[i];
  }
  return v;
}

vector3d periodic_diff(const vector3d &a, const vector3d  &b,
                          const double len[3], const int walls){
  vector3d v = b - a;
  for (int i = 0; i < 3; i++){
    if (i >= walls && len[i] > 0){
      while (v[i] > len[i]/2.0)
        v[i] -= len[i];
      while (v[i] < -len[i]/2.0)
        v[i] += len[i];
    }
  }
  return v;
}

int initialize_neighbor_tables(ball *p, int N, double neighbor_R,
                               int max_neighbors, const double len[3],
                               int walls){
  int most_neighbors = 0;
  for (int i = 0; i < N; i++){
    p[i].neighbor_center = p[i].pos;
  }
  for(int i = 0; i < N; i++){
    p[i].neighbors = new int[max_neighbors];
    p[i].num_neighbors = 0;
    for (int j = 0; j < N; j++){
      const bool is_neighbor = (i != j) &&
        (periodic_diff(p[i].pos, p[j].pos, len, walls).normsquared() <
         sqr(p[i].R + p[j].R + neighbor_R));
      if (is_neighbor){
        const int index = p[i].num_neighbors;
        p[i].num_neighbors++;
        if (p[i].num_neighbors > max_neighbors) {
          printf("Found too many neighbors: %d > %d\n", p[i].num_neighbors, max_neighbors);
          return -1;
        }
        p[i].neighbors[index] = j;
      }
    }
    most_neighbors = max(most_neighbors, p[i].num_neighbors);
  }
  return most_neighbors;
}

void update_neighbors(ball &a, int n, const ball *bs, int N,
                      double neighbor_R, const double len[3], int walls){
  a.num_neighbors = 0;
  for (int i = 0; i < N; i++){
    if ((i != n) &&
        (periodic_diff(a.pos, bs[i].neighbor_center, len,
                          walls).normsquared()
         < sqr(a.R + bs[i].R + neighbor_R))){
      a.neighbors[a.num_neighbors] = i;
      a.num_neighbors++;
    }
  }
}

inline void add_neighbor(int new_n, ball *p, int id){
  int i = p[id].num_neighbors;
  while (i > 0 && p[id].neighbors[i-1] > new_n){
    p[id].neighbors[i] = p[id].neighbors[i-1];
    i --;
  }
  p[id].neighbors[i] = new_n;
  p[id].num_neighbors ++;
}

inline void remove_neighbor(int old_n, ball *p, int id){
  int i = p[id].num_neighbors - 1;
  int temp = p[id].neighbors[i];
  while (temp != old_n){
    i --;
    const int temp2 = temp;
    temp = p[id].neighbors[i];
    p[id].neighbors[i] = temp2;
  }
  p[id].num_neighbors --;
}

void inform_neighbors(const ball &new_p, const ball &old_p, ball *p, int n){
  int new_index = 0, old_index = 0;
  while (true){
    if (new_index == new_p.num_neighbors){
      for(int i = old_index; i < old_p.num_neighbors; i++)
        remove_neighbor(n, p, old_p.neighbors[i]);
      return;
    }
    if (old_index == old_p.num_neighbors){
      for(int i = new_index; i < new_p.num_neighbors; i++)
        add_neighbor(n, p, new_p.neighbors[i]);
      return;
    }
    if (new_p.neighbors[new_index] < old_p.neighbors[old_index]){
      add_neighbor(n, p, new_p.neighbors[new_index]);
      new_index ++;
    } else if (old_p.neighbors[old_index] < new_p.neighbors[new_index]){
      remove_neighbor(n, p, old_p.neighbors[old_index]);
      old_index ++;
    } else {
      new_index ++;
      old_index ++;
    }
  }
}

bool overlap(const ball &a, const ball &b, const double len[3], int walls,
             double dr){
  const vector3d ab = periodic_diff(a.pos, b.pos, len, walls);
  return (ab.normsquared() < sqr(a.R + b.R + 2*dr));
}

int overlaps_with_any(const ball &a, const ball *p,
                      const double len[3], int walls, double dr){
  for (int i = 0; i < a.num_neighbors; i++){
    if (overlap(a,p[a.neighbors[i]],len,walls,dr))
      return true;
  }
  return false;
}

bool in_cell(const ball &p, const double len[3], const int walls,
             double dr){
  for (int i = 0; i < walls; i++){
    if (p.pos[i]-p.R-dr < 0.0 || p.pos[i]+p.R+dr > len[i])
      return false;
  }
  return true;
}

ball random_move(const ball &p, double move_scale, const double len[3]){
  ball temp = p;
  temp.pos = sw_fix_periodic(temp.pos + vector3d::ran(move_scale), len);
  return temp;
}

void move_one_ball(int id, ball *p, int N, double len[3], int walls,
                   double neighbor_R, double translation_distance,
                   double interaction_distance, int max_neighbors, double dr,
                   move_info *moves, int interactions, double *ln_energy_weights){
  moves->total++;
  moves->old_count =
    count_interactions(id, p, interaction_distance, len, walls);
  moves->new_count = moves->old_count;
  ball temp = random_move(p[id], translation_distance, len);
  if (!in_cell(temp, len, walls, dr)) return; // bad move!
  if (overlaps_with_any(temp, p, len, walls, dr)) return; // bad move!
  const bool get_new_neighbors =
    (periodic_diff(temp.pos, temp.neighbor_center, len, walls).normsquared()
     > sqr(neighbor_R/2.0));
  if (get_new_neighbors){
    // If we've moved too far, then the overlap test may have given a false
    // negative. So we'll find our new neighbors, and check against them.
    // If we still don't overlap, then we'll have to update the tables
    // of our neighbors that have changed.
    temp.neighbors = new int[max_neighbors];
    update_neighbors(temp, id, p, N, neighbor_R + 2*dr, len, walls);
    moves->updates++;
    // However, for this check (and this check only), we don't need to
    // look at all of our neighbors, only our new ones.
    // fixme: do this!
    //int *new_neighbors = new int[max_neighbors];

    if (overlaps_with_any(temp, p, len, walls, dr)) {
      // turns out we overlap after all.  :(
      delete[] temp.neighbors;
      return;
    }
  }
  // Now that we know that we are keeping the new move, and after we
  // have updated the neighbor tables if needed, we can compute the
  // new interaction count.
  ball pid = p[id]; // save a copy
  p[id] = temp; // temporarily update the position
  moves->new_count = count_interactions(id, p, interaction_distance, len, walls);
  p[id] = pid;
  // Now we can check if we actually want to do this move based on the
  // new energy.
  const double lnPmove =
    ln_energy_weights[moves->new_count - moves->old_count + interactions]
    - ln_energy_weights[interactions];
  if (lnPmove < 0) {
    const double Pmove = exp(lnPmove);
    if (random::ran() > Pmove) {
      // We want to reject this move because it is too improbable
      // based on our weights.
      moves->new_count = moves->old_count; // undo the energy change
      return;
    }
  }
  if (get_new_neighbors) {
    // Okay, we've checked twice, just like Santa Clause, so we're definitely
    // keeping this move and need to tell our neighbors where we are now.
    temp.neighbor_center = temp.pos;
    inform_neighbors(temp, p[id], p, id);
    moves->informs++;
    delete[] p[id].neighbors;
  }
  p[id] = temp; // Yay, we have a successful move!
  moves->working++;
}

int count_interactions(int id, ball *p, double interaction_distance,
                       double len[3], int walls){
  int interactions = 0;
  for(int i = 0; i < p[id].num_neighbors; i++){
    if(periodic_diff(p[id].pos, p[p[id].neighbors[i]].pos,
                     len, walls).normsquared()
       <= uipow(interaction_distance,2))
      interactions++;
  }
  return interactions;
}

int count_all_interactions(ball *balls, int N, double interaction_distance,
                           double len[3], int walls) {
  // Count initial number of interactions
  // Sum over i < k for all |ball[i].pos - ball[k].pos| < interaction_distance
  int interactions = 0;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < balls[i].num_neighbors; j++) {
      if(i < balls[i].neighbors[j]
         && periodic_diff(balls[i].pos,
                          balls[balls[i].neighbors[j]].pos,
                          len, walls).norm()
         <= interaction_distance)
        interactions++;
    }
  }
  return interactions;
}

int max_entropy_index(long *energy_histogram, double *ln_energy_weights,
                      int energy_levels){
  int max_entropy = 0;
  for (int i = 0; i < energy_levels; i++) {
    if (log(energy_histogram[i]) - ln_energy_weights[i]
        > (log(energy_histogram[max_entropy])
           - ln_energy_weights[max_entropy])) {
      max_entropy = i;
    }
  }
  return max_entropy;
}

void flush_arrays(long *energy_histogram, double *ln_energy_weights,
                  int energy_levels, bool flush_weights = true){
  int max_entropy =
    max_entropy_index(energy_histogram, ln_energy_weights, energy_levels);
  for (int i = 0; i < max_entropy; i++)
    ln_energy_weights[i] = ln_energy_weights[max_entropy];
  if(flush_weights){
    for (int i = 0; i < energy_levels; i++)
      energy_histogram[i] = 0;
  }
  return;
}

void flat_hist(long *energy_histogram, double *ln_energy_weights,
               int energy_levels){
  int max_entropy =
    max_entropy_index(energy_histogram, ln_energy_weights, energy_levels);
  for (int i = max_entropy; i < energy_levels; i++) {
    ln_energy_weights[i] -= log(energy_histogram[i] > 0 ?
                                energy_histogram[i] : 0.01);
  }
  flush_arrays(energy_histogram, ln_energy_weights, energy_levels);

  return;
}

void gaussian_hist(long *energy_histogram, double *ln_energy_weights,
                   int energy_levels){
  // a gaussian dos takes the form dos(E) = a*exp(-(E-m)^2/(2 s^2))
  // neglecting constant terms, ln_weights = (i-max_entropy)^2/(2 s^2)
  // here the variance s^2 is given by s = 2*sqrt(2 ln2) / FWHM(dos(E))
  // for FWHM we also need to find i satisfying dos[i] = dos[max_entropy] / 2
  // this function assumess no weighing has been used, so dos = energy_histogram
  int max_entropy =
    max_entropy_index(energy_histogram, ln_energy_weights, energy_levels);
  int target_value = energy_histogram[max_entropy] / 2;
  int smallest_diff = target_value;
  int current_diff = 0;
  int best_halfway[2] = {0,0};
  for(int i = 0; i < energy_levels; i++){
    if(i == max_entropy) smallest_diff = target_value;
    current_diff = abs(energy_histogram[i] - target_value);
    if(current_diff < smallest_diff){
      smallest_diff = current_diff;
      best_halfway[i > max_entropy] = i;
    }
  }
  double s = (best_halfway[1]-best_halfway[0])/(2*sqrt(2*log(2)));
  // set actual energy weights
  for(int i = max_entropy; i < energy_levels; i++)
    ln_energy_weights[i] = uipow(i-max_entropy,2)/(2*s*s);

  flush_arrays(energy_histogram, ln_energy_weights, energy_levels);
  return;
}

void walker_hist(long *energy_histogram, double *ln_energy_weights,
                 int energy_levels, long *walkers_plus,
                 long *walkers_total, move_info *moves){
  int max_entropy =
    max_entropy_index(energy_histogram, ln_energy_weights, energy_levels);
  for(int i = max_entropy; i < energy_levels; i++){
    int top = i < energy_levels-1 ? i+1 : i;
    int bottom = i > 0 ? i-1 : i;
    int dE = bottom-top; // energy = -interactions
    double df = double(walkers_plus[top]) / walkers_total[top]
      - (double(walkers_plus[bottom]) / walkers_total[bottom]);
    double df_dE = (df != 0 && !isnan(df)) ? df/dE : 1;
    double walker_density =
      double(walkers_total[i] != 0 ? walkers_total[i] : 0.01)/moves->total;
    ln_energy_weights[i] += 0.5*(log(df_dE) - log(walker_density));
  }
  for (int i = 0; i < energy_levels; i++)
    walkers_total[i] = 0;

  flush_arrays(energy_histogram, ln_energy_weights, energy_levels);
  return;
}

double count_variation(long *energy_histogram, double *ln_energy_weights,
                     int energy_levels){
  int max_entropy =
    max_entropy_index(energy_histogram, ln_energy_weights, energy_levels);
  double mean_counts = 0;
  for(int i = max_entropy; i < energy_levels; i++)
    mean_counts += energy_histogram[i];
  mean_counts /= energy_levels;
  double variation = 0;
  for(int i = max_entropy; i < energy_levels; i++)
    variation += sqr(energy_histogram[i] - mean_counts);
  variation = sqrt(variation/energy_levels)/mean_counts;
  return variation;
}

double fcc_dist(int n, int m, int l, double x, double y, double z){
  return sqrt(sqr(n+x/2) + sqr(m+y/2) + sqr(l+z/2));
}

vector3d pos_at(int n, int m, int l, double x, double y, double z, double a){
  return a*vector3d(n+x/2,m+y/2,l+z/2);
}

int max_balls_within(double distance){ // distances are all normalized to ball radius
  distance += 1e-10; // add a tiny, but necessary margin of error
  double a = 2*sqrt(2); // fcc lattice constant
  int c = int(ceil(distance/a)); // number of cubic fcc cells to go out from center
  int num = -1; // number of balls within a given radius; don't count the center ball
  for(int n  = -c; n <= c; n++){
    for(int m = -c; m <= c; m++){
      for(int l = -c; l <= c; l++){
        num += (a*fcc_dist(n,m,l,0,0,0) <= distance)
          + (a*fcc_dist(n,m,l,1,1,0) <= distance)
          + (a*fcc_dist(n,m,l,1,0,1) <= distance)
          + (a*fcc_dist(n,m,l,0,1,1) <= distance);
      }
    }
  }
  return num;
}

int maximum_interactions(int N, double interaction_distance, double neighbor_R,
                         int max_neighbors, double len[3]){
  double a = 2*sqrt(2); // fcc lattice constant in terms of ball radius
  double droplet_radius = pow(double(N),1./3); // lower bound for spherical droplet radius
  while(max_balls_within(droplet_radius) > N)
    droplet_radius -= 1; // decrease droplet radius until it is too small
  while(max_balls_within(droplet_radius) < N)
    droplet_radius += 0.01; // increase droplet size slowly to fit all balls

  int num_balls = max_balls_within(droplet_radius)+1; // number of balls in droplet
  int c = int(ceil(droplet_radius/a)); // number of cubic fcc cells to go out from center
  ball *balls = new ball[num_balls];

  int i = 0;
  for(int n  = -c; n <= c; n++){
    for(int m = -c; m <= c; m++){
      for(int l = -c; l <= c; l++){
        if(a*fcc_dist(n,m,l,0,0,0) <= droplet_radius){
          balls[i].pos = pos_at(n,m,l,0,0,0,a);
          i++;
        }
        if(a*fcc_dist(n,m,l,1,1,0) <= droplet_radius){
          balls[i].pos = pos_at(n,m,l,1,1,0,a);
          i++;
        }
        if(a*fcc_dist(n,m,l,1,0,1) <= droplet_radius){
          balls[i].pos = pos_at(n,m,l,1,0,1,a);
          i++;
        }
        if(a*fcc_dist(n,m,l,0,1,1) <= droplet_radius){
          balls[i].pos = pos_at(n,m,l,0,1,1,a);
          i++;
        }
      }
    }
  }

  while(num_balls > N){
    double max_distance = 0;
    int max_index = -1; // index of ball that is the furthest out from the droplet
    for(int i = 0; i < num_balls; i++){
      if(balls[i].pos.norm() > max_distance){
        max_distance = balls[i].pos.norm();
        max_index = i;
      }
    }
    balls[max_index].pos = balls[num_balls-1].pos;
    num_balls--;
  }

  ball *droplet_balls = new ball[N];
  for(int i = 0; i < N; i++)
    droplet_balls[i].pos = balls[i].pos;

  initialize_neighbor_tables(droplet_balls,N,neighbor_R,max_neighbors,len,0);
  int interactions = count_all_interactions(droplet_balls,N,interaction_distance,len,0);
  delete[] balls;
  delete[] droplet_balls;
  return interactions;
}
