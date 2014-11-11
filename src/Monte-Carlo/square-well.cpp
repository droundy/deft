#include <stdlib.h>
#include "Monte-Carlo/square-well.h"
#include "handymath.h"

ball::ball(){
  pos = vector3d();
  R = 1;
  neighbors = 0;
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

vector3d periodic_diff(const vector3d &a, const vector3d  &b, const double len[3],
                       const int walls){
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

int initialize_neighbor_tables(ball *p, int N, double neighbor_R, int max_neighbors,
                               const double len[3], int walls){
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
                  int energy_levels, bool flush_histogram = true, int max_entropy = 0){
  // make weights flat at energies above the maximum entropy state
  // (which correspond to boring "negative" temperatures.
  if (max_entropy == 0)
    max_entropy = max_entropy_index(energy_histogram, ln_energy_weights, energy_levels);
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

void walker_hist(long *energy_histogram, double *ln_energy_weights,
                 int energy_levels, long *walkers_up,
                 long *walkers_total, move_info *moves){
  int max_entropy =
    max_entropy_index(energy_histogram, ln_energy_weights, energy_levels);
  for(int i = max_entropy; i < energy_levels; i++){
    int top = i < energy_levels-1 ? i+1 : i;
    int bottom = i > 0 ? i-1 : i;
    int dE = bottom-top; // energy = -interactions
    double df = double(walkers_up[top]) / walkers_total[top]
      - (double(walkers_up[bottom]) / walkers_total[bottom]);
    double df_dE = (df != 0 && df == df) ? df/dE : 1;
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
                       int energy_levels, int max_entropy){
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
  printf("Highest/lowest = %.2g (%d) / %.2g (%d) mean %.2g    total %g\n",
         highest, highest_i, lowest, lowest_i, mean_counts, total_counts);
  printf("    (maxE = %d, visited %d, max_entropy %d)\n",
         maxE, num_visited, max_entropy);
  return mean_counts/lowest - 1;
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

// sw_simulation methods

void sw_simulation::move_a_ball() {
  int id = moves.total % N;
  moves.total++;
  moves.old_count = count_interactions(id, balls, interaction_distance, len, walls);
  moves.new_count = moves.old_count;
  ball temp = random_move(balls[id], translation_distance, len);
  if (!in_cell(temp, len, walls, dr)) return; // bad move!
  if (overlaps_with_any(temp, balls, len, walls, dr)) return; // bad move!
  const bool get_new_neighbors =
    (periodic_diff(temp.pos, temp.neighbor_center, len, walls).normsquared()
     > sqr(neighbor_R/2.0));
  if (get_new_neighbors){
    // If we've moved too far, then the overlap test may have given a false
    // negative. So we'll find our new neighbors, and check against them.
    // If we still don't overlap, then we'll have to update the tables
    // of our neighbors that have changed.
    temp.neighbors = new int[max_neighbors];
    update_neighbors(temp, id, balls, N, neighbor_R + 2*dr, len, walls);
    moves.updates++;
    // However, for this check (and this check only), we don't need to
    // look at all of our neighbors, only our new ones.
    // fixme: do this!
    //int *new_neighbors = new int[max_neighbors];

    if (overlaps_with_any(temp, balls, len, walls, dr)) {
      // turns out we overlap after all.  :(
      delete[] temp.neighbors;
      return;
    }
  }
  // Now that we know that we are keeping the new move, and after we
  // have updated the neighbor tables if needed, we can compute the
  // new interaction count.
  ball pid = balls[id]; // save a copy
  balls[id] = temp; // temporarily update the position
  moves.new_count = count_interactions(id, balls, interaction_distance, len, walls);
  balls[id] = pid;
  // Now we can check if we actually want to do this move based on the
  // new energy.
  const double lnPmove =
    ln_energy_weights[moves.new_count - moves.old_count + interactions]
    - ln_energy_weights[interactions];
  if (lnPmove < 0) {
    const double Pmove = exp(lnPmove);
    if (random::ran() > Pmove) {
      // We want to reject this move because it is too improbable
      // based on our weights.
      moves.new_count = moves.old_count; // undo the energy change
      if (get_new_neighbors) delete[] temp.neighbors;
      return;
    }
  }
  if (get_new_neighbors) {
    // Okay, we've checked twice, just like Santa Clause, so we're definitely
    // keeping this move and need to tell our neighbors where we are now.
    temp.neighbor_center = temp.pos;
    inform_neighbors(temp, balls[id], balls, id);
    moves.informs++;
    delete[] balls[id].neighbors;
  }
  balls[id] = temp; // Yay, we have a successful move!
  moves.working++;

  // Update interaction count and energy histogram
  interactions += moves.new_count - moves.old_count;
  energy_histogram[interactions]++;

  // increment iteraction count
  if(moves.total % N == 0) iteration++;

  update_state_search_info();
}

void sw_simulation::update_state_search_info() {
  // update round trip observations
  if(interactions <= state_of_max_entropy){
    for(int i = state_of_max_entropy; i < energy_levels; i++)
      seeking_energy[i] = true;
  }
  else if(seeking_energy[interactions]) {
    seeking_energy[interactions] = false;
    round_trips[interactions]++;
  }
  // update walker counts
  if(interactions > max_observed_interactions){
    if(interactions > max_observed_interactions) {
      max_observed_interactions = interactions;
      for(int i = state_of_max_entropy; i < energy_levels; i++)
        walkers_up[i] = 0;
    }
  }
  if(!seeking_energy[max_observed_interactions])
    walkers_up[interactions]++;
  walkers_total[interactions]++;
}

// iterate enough times for the energy to change n times.  Return the
// number of "up" moves.
int sw_simulation::simulate_energy_changes(int num_moves) {
  int num_moved = 0, num_up = 0;
  while (num_moved < num_moves) {
    int old_interactions = interactions;
    move_a_ball();
    if (interactions != old_interactions) {
      num_moved++;
      if (interactions > old_interactions) num_up++;
    }
  }
  return num_up;
}

int sw_simulation::initialize_max_entropy_and_translation_distance(double acceptance_goal) {
  int num_moves = 500;
  const double mean_allowance = 1.0;
  double dscale = 0.1;
  const int starting_iterations = iteration;
  int attempts_to_go = 4; // always try this many times, to get translation_distance right.
  double variance = 0;
  double last_energy = 0, mean = 0;
  int counted_in_mean;
  while (fabs(last_energy - interactions) > max(2,mean_allowance*sqrt(variance)) ||
         attempts_to_go >= 0) {
    attempts_to_go -= 1;
    last_energy = interactions;
    simulate_energy_changes(num_moves);
    mean = 0;
    counted_in_mean = 0;
    for (int i = 0; i < energy_levels; i++) {
      counted_in_mean += energy_histogram[i];
      mean += i*energy_histogram[i];
    }
    mean /= counted_in_mean;
    variance = 0;
    for (int i = 0; i < energy_levels; i++) {
      variance += (i-mean)*(i-mean)*energy_histogram[i];
      energy_histogram[i] = 0; // clear it out for the next attempt
    }
    variance /= counted_in_mean;

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

    // printf("Figure of merit:  %4.2g (state %.1f from %g(%d), stdev %.1f) [acceptance rate %.2f, scale %.2g]\n",
    //        (last_energy - interactions)/max(2,mean_allowance*sqrt(variance)), mean, last_energy,
    //        interactions, sqrt(variance), acceptance_rate, translation_distance);
    num_moves *= 2;
  }
  printf("Took %ld iterations to find max entropy state:  %d with width %.2g\n",
         iteration - starting_iterations, state_of_max_entropy, sqrt(variance));
  return int(mean + 0.5);
}

// initialize the weight array using the Gaussian approximation.
double sw_simulation::initialize_gaussian(double scale) {
  double variance = 0, old_variance;
  double mean = 0, old_mean;
  int counted_in_mean;
  int num_energy_moves = 200;
  const int starting_iterations = iteration;
  const double fractional_precision_required = 0.2;
  do {
    simulate_energy_changes(num_energy_moves);
    num_energy_moves *= 2; // if we haven't run long enough, run longer next time!

    old_variance = variance;
    old_mean = mean;
    mean = 0;
    counted_in_mean = 0;
    for (int i = 0; i < energy_levels; i++) {
      counted_in_mean += energy_histogram[i];
      mean += i*energy_histogram[i];
    }
    mean /= counted_in_mean;
    variance = 0;
    for (int i = 0; i < energy_levels; i++)
      variance += (i-mean)*(i-mean)*energy_histogram[i];
    variance /= counted_in_mean;

    // Keep simulating until the mean has not moved by more than a
    // tiny fraction of the standard deviation.  Also we keep going if
    // the standard deviation has changed much.  This is an attempt to
    // avoid scenarios where we haven't run long enough to adequately
    // sample the variance.
    //printf("mean is %4.1f with stdev %5.4g\tdifference %.2g vs %.2g\n", mean, sqrt(variance),
    //       variance - old_variance, fractional_precision_required*sqrt(variance));
  } while (fabs(old_mean - mean) > fractional_precision_required*sqrt(variance) ||
           fabs(sqrt(old_variance) - sqrt(variance))
             > fractional_precision_required*sqrt(variance));

  for (int i = state_of_max_entropy; i < energy_levels; i++)
    ln_energy_weights[i] -= scale*exp(-uipow(i-mean,2)/(scale*2*variance));
  for (int i = 0; i < state_of_max_entropy; i++)
    ln_energy_weights[i] -= scale*exp(-uipow(state_of_max_entropy-mean,2)/(scale*2*variance));

  printf("Took %ld iterations to find mean and variance for gaussian:  %g with width %g\n",
         iteration - starting_iterations, mean, sqrt(variance));
  return sqrt(variance);
}

// initialize the weight array using the specified temperature.
void sw_simulation::initialize_canonical(double kT) {
  for(int i=0; i < energy_levels; i++){
    ln_energy_weights[i] = i/kT;
  }
}

// improve the weight array using the Wang-Landau method.
void sw_simulation::initialize_wang_landau(double wl_factor, double wl_threshold,
                                           double wl_cutoff) {
  int weight_updates = 0;
  const int starting_iterations = iteration;
  while (1) {
    for (int i=0; i < 100*N*energy_levels; i++) {
      move_a_ball();
      ln_energy_weights[interactions] -= wl_factor;
    }

    // check whether our histogram is flat enough to update wl_factor
    const double variation = count_variation(energy_histogram, ln_energy_weights,
                                             energy_levels, state_of_max_entropy);

    // print status text for testing purposes
    printf("\nweight update: %i\n",weight_updates++);
    printf("  WL factor: %g\n",wl_factor);
    printf("  count variation: %g\n", variation);
    fflush(stdout);

    if (variation > 0 && variation < max(wl_threshold, exp(wl_factor))) {
      wl_factor /= 2;
      // for wang-landau, only flush energy histogram; keep weights
      flush_arrays(energy_histogram, ln_energy_weights, energy_levels, true,
                   state_of_max_entropy);

      // repeat until terminal condition is met
      if (wl_factor < wl_cutoff) {
        printf("Took %ld iterations to initialize with Wang-Landau method\n",
               iteration - starting_iterations);
        return;
      }
    }
  }
}
