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
                  int energy_levels, bool flush_histogram = true){
  // make weights flat at energies above the maximum entropy state
  // (which correspond to boring "negative" temperatures.
  int max_entropy =
    max_entropy_index(energy_histogram, ln_energy_weights, energy_levels);
  for (int i = 0; i < max_entropy; i++)
    ln_energy_weights[i] = ln_energy_weights[max_entropy];
  // Now let's just add the minimum, so our weights don't get
  // out of hand (which could lead to increased roundoff errors).
  double min_weight = ln_energy_weights[0];
  for (int i = 0; i < max_entropy; i++)
    min_weight = min(min_weight, ln_energy_weights[i]);
  for (int i = 0; i < energy_levels; i++)
    ln_energy_weights[i] -= min_weight;
  // Finally, if requested, zero out the histogram.
  if (flush_histogram) {
    for (int i = 0; i < energy_levels; i++) {
      if (energy_histogram[i] > 0) energy_histogram[i] = 1;
    }
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
  double lowest = 1e200, highest = 0;
  int maxE = 0; // tracks the lowest (largest magnitude) energy yet seen
  int highest_i = 0; // the most commonly visited energy
  int lowest_i = 0; // the least commonly visited energy
  int num_visited = 0; // the total number of distinct energies visited.
  double total_counts = 0;
  for(int i = max_entropy; i < energy_levels; i++) {
    if (energy_histogram[i] > 0) {
      num_visited++;
      total_counts += energy_histogram[i];
      if (energy_histogram[i] > highest) {
        highest = energy_histogram[i];
        highest_i = i;
      }
      if (energy_histogram[i] < lowest) {
        lowest = energy_histogram[i];
        lowest_i = i;
      }
      maxE = i;
    }
  }
  double mean_counts = total_counts / num_visited;
  double variation = 0;
  for(int i = max_entropy; i < energy_levels; i++)
    variation += sqr(energy_histogram[i] - mean_counts);
  variation = sqrt(variation/energy_levels)/mean_counts;
  printf("Highest/lowest = %.2g (%d) / %.2g (%d) mean %.2g\n",
         highest, highest_i, lowest, lowest_i, mean_counts);
  printf("    (ratio = %g, maxE = %d, visited %d, max_entropy %d)\n",
         mean_counts/lowest, maxE, num_visited, max_entropy);
  return mean_counts/lowest - 1; // or alternatively variation;
}

vector3d fcc_pos(int n, int m, int l, double x, double y, double z, double a){
  return a*vector3d(n+x/2,m+y/2,l+z/2);
}

int max_balls_within(double distance){ // distances are all normalized to ball radius
  distance += 1e-10; // add a tiny, but necessary margin of error
  double a = 2*sqrt(2); // fcc lattice constant
  int c = int(ceil(distance/a)); // number of cubic fcc cells to go out from center
  int num = -1; // number of balls within a given radius; don't count the center ball

  int xs[] = {0,1,1,0};
  int ys[] = {0,1,0,1};
  int zs[] = {0,0,1,1};
  for(int n  = -c; n <= c; n++){
    for(int m = -c; m <= c; m++){
      for(int l = -c; l <= c; l++){
        for(int k = 0; k < 4; k++)
          num += (fcc_pos(n,m,l,xs[k],ys[k],zs[k],a).norm() <= distance);
      }
    }
  }
  return num;
}

int maximum_interactions(int N, double interaction_distance, double neighbor_R,
                         int max_neighbors, double len[3]){
  // Find size of smallest droplet with N balls
  double droplet_radius_lower = pow(double(N),1./3); // lower bound for droplet radius
  double droplet_radius_upper = droplet_radius_lower*sqrt(3); // upper bound
  double droplet_radius = droplet_radius_lower;
  bool too_small = true;

  while(droplet_radius_upper - droplet_radius_lower > 1e-3){
    if(too_small)
      droplet_radius_lower = droplet_radius;
    else
      droplet_radius_upper = droplet_radius;
    droplet_radius = (droplet_radius_lower + droplet_radius_upper)/2;
    too_small = max_balls_within(droplet_radius) < N;
  }
  if(too_small)
    droplet_radius = droplet_radius_upper;

  double a = 2*sqrt(2); // fcc lattice constant in terms of ball radius
  int num_balls = max_balls_within(droplet_radius)+1; // number of balls in droplet
  int c = int(ceil(droplet_radius/a)); // number of cubic fcc cells to go out from center
  ball *balls = new ball[num_balls];

  int i = 0;
  int xs[] = {0,1,1,0};
  int ys[] = {0,1,0,1};
  int zs[] = {0,0,1,1};
  for(int n  = -c; n <= c; n++){
    for(int m = -c; m <= c; m++){
      for(int l = -c; l <= c; l++){
        for(int k = 0; k < 4; k++){
          if(fcc_pos(n,m,l,xs[k],ys[k],zs[k],a).norm() <= droplet_radius){
            balls[i].pos = fcc_pos(n,m,l,xs[k],ys[k],zs[k],a);
            i++;
          }
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
    // place the ball farthest from the center in place of that with the highest index,
    //   and forget about the ball with the highest index
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


// sw_simulation methods

void sw_simulation::move_a_ball() {
  move_one_ball(iteration % N, balls, N, len, walls, neighbor_R, translation_distance,
                interaction_distance, max_neighbors, dr, &moves, interactions,
                ln_energy_weights);
  interactions += moves.new_count - moves.old_count;
  energy_histogram[interactions]++;
  iteration++; // increment the number of "iterations"
}

void sw_simulation::initialize_max_entropy_and_translation_distance(double acceptance_goal) {
  int num_up = 0;
  const int num_moves = 500;
  const double deviation_allowed = 0.1/sqrt(num_moves);
  double dscale = 0.1;
  const int starting_iterations = iteration;
  int attempts_to_go = 4; // always try this many times, to get translation_distance right.
  double variance;
  while (fabs(num_up/double(num_moves) - 0.5) > deviation_allowed || attempts_to_go-- >= 0) {
    num_up = 0;
    int num_moved = 0;
    for (int i=0; i<energy_levels; i++) energy_histogram[i] = 0;
    while (num_moved < num_moves) {
      int old_interactions = interactions;
      move_a_ball();
      if (interactions != old_interactions) {
        num_moved++;
        if (interactions > old_interactions) num_up++;
      }
    }
    double mean = 0;
    double counted_in_mean = 0;
    for (int i=0; i<energy_levels; i++) {
      counted_in_mean += energy_histogram[i];
      mean += i*energy_histogram[i];
    }
    mean /= counted_in_mean;
    variance = 0;
    for (int i=0; i<energy_levels; i++) variance += (i-mean)*(i-mean)*energy_histogram[i];
    variance /= counted_in_mean;
    state_of_max_entropy = int(mean + 0.5);

    // ---------------------------------------------------------------
    // Fine-tune translation scale to reach acceptance goal
    // ---------------------------------------------------------------
    const double acceptance_rate =
      double(moves.working-moves.working_old)/(moves.total-moves.total_old);
    moves.working_old = moves.working;
    moves.total_old = moves.total;
    if (acceptance_rate < acceptance_goal)
      translation_distance /= 1+dscale;
    else
      translation_distance *= 1+dscale;
    // hokey heuristic for tuning dscale
    const double closeness = fabs(acceptance_rate - acceptance_goal)/acceptance_rate;
    if(closeness > 0.5) dscale *= 2;
    else if(closeness < dscale*2 && dscale > 0.01) dscale/=2;

    // printf("Figure of merit:  %g (state %d)\n",
    //        (num_up/double(num_moves) - 0.5)/deviation_allowed,
    //        state_of_max_entropy);
  }
  printf("Took %ld iterations to find max entropy state:  %d with width %g\n",
         iteration - starting_iterations, state_of_max_entropy, sqrt(variance));
}
