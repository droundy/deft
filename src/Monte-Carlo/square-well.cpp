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
  // for (int i = 0; i < 3; i++){
  //   while (v[i] > len[i])
  //     v[i] -= len[i];
  //   while (v[i] < 0.0)
  //     v[i] += len[i];
  // }
  if (v.x > len[0]) v.x -= len[0];
  else if (v.x < 0) v.x += len[0];
  if (v.y > len[1]) v.y -= len[1];
  else if (v.y < 0) v.y += len[1];
  if (v.z > len[2]) v.z -= len[2];
  else if (v.z < 0) v.z += len[2];
  return v;
}

vector3d periodic_diff(const vector3d &a, const vector3d  &b, const double len[3],
                       const int walls){
  vector3d v = b - a;
  // for (int i = walls; i < 3; i++){
  //   if (len[i] > 0){
  //     while (v[i] > len[i]/2.0)
  //       v[i] -= len[i];
  //     while (v[i] < -len[i]/2.0)
  //       v[i] += len[i];
  //   }
  // }
  if (2 >= walls) {
    if (v.z > 0.5*len[2]) v.z -= len[2];
    else if (v.z < -0.5*len[2]) v.z += len[2];
    if (1 >= walls) {
      if (v.y > 0.5*len[1]) v.y -= len[1];
      else if (v.y < -0.5*len[1]) v.y += len[1];
      if (0 >= walls) {
        if (v.x > 0.5*len[0]) v.x -= len[0];
        else if (v.x < -0.5*len[0]) v.x += len[0];
      }
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

bool overlap(const ball &a, const ball &b, const double len[3], int walls){
  const vector3d ab = periodic_diff(a.pos, b.pos, len, walls);
  return (ab.normsquared() < sqr(a.R + b.R));
}

int overlaps_with_any(const ball &a, const ball *p, const double len[3], int walls){
  for (int i = 0; i < a.num_neighbors; i++){
    if (overlap(a,p[a.neighbors[i]],len,walls))
      return true;
  }
  return false;
}

bool in_cell(const ball &p, const double len[3], const int walls){
  for (int i = 0; i < walls; i++){
    if (p.pos[i]-p.R < 0.0 || p.pos[i]+p.R > len[i])
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

int new_max_entropy_state(long *energy_histogram, double *ln_energy_weights,
                      int energy_levels){
  int max_entropy_state = 0;
  for (int i = 0; i < energy_levels; i++) {
    if (log(energy_histogram[i]) - ln_energy_weights[i]
        > (log(energy_histogram[max_entropy_state])
           - ln_energy_weights[max_entropy_state])) {
      max_entropy_state = i;
    }
  }
  return max_entropy_state;
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
  const int old_interaction_count =
    count_interactions(id, balls, interaction_distance, len, walls);
  ball temp = random_move(balls[id], translation_scale, len);
  // If we're out of the cell or we overlap, this is a bad move!
  if (!in_cell(temp, len, walls) || overlaps_with_any(temp, balls, len, walls)){
    end_move_updates();
    return;
  }
  const bool get_new_neighbors =
    (periodic_diff(temp.pos, temp.neighbor_center, len, walls).normsquared()
     > sqr(neighbor_R/2.0));
  if (get_new_neighbors){
    // If we've moved too far, then the overlap test may have given a false
    // negative. So we'll find our new neighbors, and check against them.
    // If we still don't overlap, then we'll have to update the tables
    // of our neighbors that have changed.
    temp.neighbors = new int[max_neighbors];
    update_neighbors(temp, id, balls, N, neighbor_R, len, walls);
    moves.updates++;
    // However, for this check (and this check only), we don't need to
    // look at all of our neighbors, only our new ones.
    // fixme: do this!
    //int *new_neighbors = new int[max_neighbors];

    if (overlaps_with_any(temp, balls, len, walls)) {
      // turns out we overlap after all.  :(
      delete[] temp.neighbors;
      end_move_updates();
      return;
    }
  }
  // Now that we know that we are keeping the new move (unless the
  // weights say otherwise), and after we have updated the neighbor
  // tables if needed, we can compute the new interaction count.
  ball pid = balls[id]; // save a copy
  balls[id] = temp; // temporarily update the position
  const int new_interaction_count =
    count_interactions(id, balls, interaction_distance, len, walls);
  balls[id] = pid;
  // Now we can check whether we actually want to do this move based on the
  // new energy.
  const int energy_change = new_interaction_count - old_interaction_count;
  transitions(energy, energy_change) += 1; // update the transition histogram
  const double lnPmove =
    ln_energy_weights[energy + energy_change] - ln_energy_weights[energy];
  if (lnPmove < 0) {
    const double Pmove = exp(lnPmove);
    if (random::ran() > Pmove) {
      // We want to reject this move because it is too improbable
      // based on our weights.
      if (get_new_neighbors) delete[] temp.neighbors;
      end_move_updates();
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
  energy += energy_change;

  if(energy_change != 0) energy_change_updates();

  end_move_updates();
}

void sw_simulation::end_move_updates(){
   // update iteration counter, energy histogram, and walker counters
  if(moves.total % N == 0) iteration++;
  energy_histogram[energy]++;
  if(energy_observed[min_energy_state])
    walkers_up[energy]++;
  walkers_total[energy]++;
}

void sw_simulation::energy_change_updates(){
  // If we're at or above the state of max entropy, we have not yet observed any energies
  if(energy <= max_entropy_state){
    for(int i = max_entropy_state; i < energy_levels; i++)
      energy_observed[i] = false;
  }
  // If we have not yet observed this energy since the last time we were at max entropy,
  //   now we have!
  else if(!energy_observed[energy]){
    energy_observed[energy] = true;
    samples[energy]++;
  }
  // If we have observed a new lowest energy, remember it; furthermore, this means we can't
  //   have had any walkers going up from the lowest energy, so we reset walkers_up
  if(energy > min_energy_state){
    min_energy_state = energy;
    for(int i = max_entropy_state; i < energy_levels; i++)
      walkers_up[i] = 0;
  }
}

int sw_simulation::simulate_energy_changes(int num_moves) {
  int num_moved = 0, num_up = 0;
  while (num_moved < num_moves) {
    int old_energy = energy;
    move_a_ball();
    if (energy != old_energy) {
      num_moved++;
      if (energy > old_energy) num_up++;
    }
  }
  return num_up;
}

void sw_simulation::flush_weight_array(){
  for (int i = 0; i < max_entropy_state; i++)
    ln_energy_weights[i] = ln_energy_weights[max_entropy_state];
  double min_weight = ln_energy_weights[0];
  for (int i = 1; i < max_entropy_state; i++)
    min_weight = min(min_weight, ln_energy_weights[i]);
  for (int i = 0; i < energy_levels; i++)
    ln_energy_weights[i] -= min_weight;
  return;
}

int sw_simulation::initialize_max_entropy_and_translation_distance(double acceptance_goal) {
  printf("Moving to most probable state and tuning translation distance.\n");
  int num_moves = 500;
  const double mean_allowance = 1.0;
  double dscale = 0.1;
  const int starting_iterations = iteration;
  int attempts_to_go = 4; // always try this many times, to get translation_scale right.
  double mean = 0, variance = 0;
  int last_energy = 0;
  int counted_in_mean;
  while (abs(last_energy-energy) > max(2,mean_allowance*sqrt(variance)) ||
         attempts_to_go >= 0) {
    attempts_to_go -= 1;
    last_energy = energy;
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
      translation_scale /= 1+dscale;
    else
      translation_scale *= 1+dscale;
    // hokey heuristic for tuning dscale
    const double closeness = fabs(acceptance_rate - acceptance_goal)/acceptance_rate;
    if(closeness > 0.5) dscale *= 2;
    else if(closeness < dscale*2 && dscale > 0.01) dscale/=2;

    num_moves *= 2;
  }
  printf("Took %ld iterations to find most probable state: %d with width %.2g\n",
         iteration - starting_iterations, max_entropy_state, sqrt(variance));
  return int(mean + 0.5);
}

void sw_simulation::initialize_translation_distance(double acceptance_goal) {
  printf("Tuning translation distance.\n");
  const int max_tries = 5;
  const int num_moves = 100*N*N;
  const int starting_iterations = iteration;
  double dscale = 0.1;
  double acceptance_rate = 0;
  for (int num_tries=0;
       fabs(acceptance_rate - acceptance_goal) > 0.05 && num_tries < max_tries;
       num_tries++) {
    // ---------------------------------------------------------------
    // Fine-tune translation scale to reach acceptance goal
    // ---------------------------------------------------------------
    for (int i=0;i<num_moves;i++) move_a_ball();
    acceptance_rate =
      double(moves.working-moves.working_old)/(moves.total-moves.total_old);
    moves.working_old = moves.working;
    moves.total_old = moves.total;
    if (acceptance_rate < acceptance_goal)
      translation_scale /= 1+dscale;
    else
      translation_scale *= 1+dscale;
    // hokey heuristic for tuning dscale
    const double closeness = fabs(acceptance_rate - acceptance_goal)/acceptance_rate;
    if(closeness > 0.5) dscale *= 2;
    else if(closeness < dscale*2 && dscale > 0.01) dscale/=2;
  }
  printf("Took %ld iterations to find acceptance rate of %.2g with translation scale %.2g\n",
         iteration - starting_iterations,
         acceptance_rate, translation_scale);
}

// initialize the weight array using the Gaussian method
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
  } while (fabs(old_mean - mean) > fractional_precision_required*sqrt(variance) ||
           fabs(sqrt(old_variance) - sqrt(variance))
             > fractional_precision_required*sqrt(variance));

  for (int i = max_entropy_state; i < energy_levels; i++)
    ln_energy_weights[i] -= scale*exp(-uipow(i-mean,2)/(scale*2*variance));
  for (int i = 0; i < max_entropy_state; i++)
    ln_energy_weights[i] -= scale*exp(-uipow(max_entropy_state-mean,2)/(scale*2*variance));

  printf("Took %ld iterations to find gaussian: mean %g; width %g\n",
         iteration - starting_iterations, mean, sqrt(variance));
  return sqrt(variance);
}

// initialize the weight array using the specified temperature.
void sw_simulation::initialize_canonical(double kT) {
  for(int i=0; i < energy_levels; i++){
    ln_energy_weights[i] = i/kT;
  }
}

// initialize the weight array using the Wang-Landau method.
void sw_simulation::initialize_wang_landau(double wl_factor, double wl_fmod,
                                           double wl_threshold, double wl_cutoff) {
  int weight_updates = 0;
  bool finished_initializing = false;
  int last_min_energy_state = min_energy_state;
  const int starting_iterations = iteration;
  while (!finished_initializing) {
    for (int i=0; i < N*energy_levels; i++) {
      move_a_ball();
      ln_energy_weights[energy] -= wl_factor;
    }

    // compute variation in energy histogram
    int highest_hist_i = 0; // the most commonly visited energy
    int lowest_hist_i = 0; // the least commonly visited energy
    double highest_hist = 0; // highest histogram value
    double lowest_hist = 1e200; // lowest hitsogram value
    double total_counts = 0; // total counts in energy histogram
    for(int i = max_entropy_state+1; i < energy_levels; i++){
      if(energy_histogram[i] > 0){
        total_counts += energy_histogram[i];
        if(energy_histogram[i] > highest_hist){
          highest_hist = energy_histogram[i];
          highest_hist_i = i;
        }
        if(energy_histogram[i] < lowest_hist){
          lowest_hist = energy_histogram[i];
          lowest_hist_i = i;
        }
      }
    }
    double hist_mean = total_counts / (min_energy_state-max_entropy_state);
    const double variation = hist_mean/lowest_hist - 1;

    // print status text for testing purposes
    if(printing_allowed()){
      printf("\nWL weight update: %i\n",weight_updates++);
      printf("  WL factor: %g\n",wl_factor);
      printf("  count variation: %g\n", variation);
      printf("  highest/lowest histogram energies (values): %d (%.2g) / %d (%.2g)\n",
             highest_hist_i, highest_hist, lowest_hist_i, lowest_hist);
      printf("  minimum_energy_state: %d,  max_entropy_state: %d,  energies visited: %d\n",
             min_energy_state, max_entropy_state, min_energy_state-max_entropy_state);
      printf("  current energy: %d\n", energy);
      printf("  hist_mean: %g,  total_counts: %g\n", hist_mean, total_counts);
    }

    // check whether our histogram is flat enough to update wl_factor
    if (variation > 0 && variation < max(wl_threshold, exp(wl_factor))) {
      wl_factor /= wl_fmod;
      flush_weight_array();
      for (int i = 0; i < energy_levels; i++) {
        if (energy_histogram[i] > 0) energy_histogram[i] = 1;
      }

      // repeat until terminal condition is met,
      // and make sure we're not stuck at a newly introduced minimum energy state
      if (wl_factor < wl_cutoff && last_min_energy_state == min_energy_state) {
        printf("Took %ld iterations and %i updates to initialize with Wang-Landau method.\n",
               iteration - starting_iterations, weight_updates);
        finished_initializing = true;
      }
      last_min_energy_state = min_energy_state;
    }
  }
}

// initialize the weight array using the optimized ensemble method.
void sw_simulation::initialize_walker_optimization(int first_update_iterations,
                                                   int init_min_energy_samples){
  int weight_updates = 0;
  while(iteration <= first_update_iterations ||
        samples[min_energy_state] < init_min_energy_samples) {
    // run an iteration
    for(int i = 0; i < N; i++) move_a_ball();
    assert(energy ==
           count_all_interactions(balls, N, interaction_distance, len, walls));
    // update weights when we've run long enough
    if((iteration - first_update_iterations)
       % int(first_update_iterations*uipow(2,weight_updates)) == 0){
      for(int i = max_entropy_state; i < energy_levels; i++){
        int top = i < energy_levels-1 ? i+1 : i;
        int bottom = i > 0 ? i-1 : i;
        int dE = bottom-top;
        double df = double(walkers_up[top]) / walkers_total[top]
          - (double(walkers_up[bottom]) / walkers_total[bottom]);
        double df_dE = (df < 0 && df == df) ? df/dE : 1;
        double walker_density =
          double(walkers_total[i] != 0 ? walkers_total[i] : 0.01)/moves.total;
        ln_energy_weights[i] += 0.5*(log(df_dE) - log(walker_density));
      }
      flush_weight_array();
      for(int i = 0; i < energy_levels; i++){
        walkers_up[i] = 0;
        walkers_total[i] = 0;
      }
      printf("Weight update: %i. min_energy_state: %i. samples: %li\n",
             weight_updates, min_energy_state, samples[min_energy_state]);
      weight_updates++;
    }
  }
}

void sw_simulation::initialize_robustly_optimistic(double robust_update_scale,
                                                   double robust_cutoff){
  bool done;
  double mean_hist = 0;
  do {
    done = true; // Let's be optimistic!
    // First, let's reset our weights based on what we already know!
    for (int e = max_entropy_state; e < energy_levels; e++) {
      if (energy_histogram[e])
        ln_energy_weights[e] -= robust_update_scale*log(energy_histogram[e]);
    }
    flush_weight_array();
    /* Check whether we have *ever* seen the minimum energy in the
       last time through.  If not, we absolutely shouldn't be so
       optimistic as to believe that this time we'll have everything
       right. */
    if(energy_histogram[min_energy_state] == 0) done = false;
    // Now reset the calculation!
    for (int e = 0; e < energy_levels; e++) {
      energy_histogram[e] = 0;
      samples[e] = 0;
      energy_observed[e] = false;
    }
    // Now let's run a while to see if we can find another answer.
    // fixme: SUPER SCARY WARNING:
    // N^3 scaling is probably VERY BAD!!! but it works really well for now...
    const int test_iterations = energy_levels*uipow(N,3);
    for (int i = 0; i < test_iterations; i++){
      for (int i = 0; i < N; i++)
        move_a_ball();
    }
    // Dump energy histogram to the console so that we can see how we're doing
    if(printing_allowed()){
      for(int i = max_entropy_state; i <= min_energy_state; i++)
        printf("%i  %li\n",i,energy_histogram[i]);
    }
    // Check whether our histogram is sufficiently flat; if not, we have to restart!
    mean_hist = test_iterations*N/double(min_energy_state - max_entropy_state);
    for (int e = max_entropy_state; e <= min_energy_state; e++) {
      if (energy_histogram[e] < N
          || energy_histogram[e] < robust_cutoff*mean_hist) {
        printf("After %d iterations, we fail at at energy %d widh hist = %ld"
               " vs. mean hist %g.\n",
               test_iterations*N, e, energy_histogram[e], mean_hist);
        done = false;
        break;
      }
    }
  } while (!done);
}


void sw_simulation::initialize_bubble_suppression(double bubble_scale, double bubble_cutoff){
  double width;
  double range;
  do {
    initialize_max_entropy_and_translation_distance(); // get to point of max entropy
    width = initialize_gaussian(bubble_scale); // subtract off gaussian
    range = min_energy_state - max_entropy_state;
    printf("*** Gaussian has width %.1f and range %.0f (ratio %.2f)\n",
           width, range, width/range);
  } while (width < bubble_cutoff*range);
}

/* Update the weight array using the transition matrix to make the
   downward rates equal, which should optimize the round-trip
   time. */
void sw_simulation::update_weights_using_transition_flux(double fractional_precision) {
  update_weights_using_transitions(fractional_precision);
  const int energies_observed = min_energy_state+1;

  double *wanted_n = new double[energies_observed];

  for (int i = max_entropy_state; i < energies_observed; i++) {
    wanted_n[i] = 1.0;
    for (int j=i+1; j < energies_observed; j++) {
      for (int k=0; k < i; k++) {
        /* counting transitions from the high energy j to a lower
           energy j.  Ignoring any transitions from energies above the
           maximum entropy point.  It would be nice (in principle) to
           include those, but I don't think it matters. */
        wanted_n[i] -= transition_matrix(j, k)*wanted_n[k];
      }
    }
    double norm = 0;
    for (int j=i+1; j < energies_observed; j++) {
      norm += transition_matrix(j, i);
    }
    wanted_n[i] /= norm;
  }

  for (int i=max_entropy_state; i<energies_observed;i++) {
    /* Assuming we have already updated the weights for a flat
       histogram, this should adjust the weights to give a flat
       downward transition rate. */
    ln_energy_weights[i] += log(wanted_n[i]);
  }
}

// update the weight array using transitions
void sw_simulation::update_weights_using_transitions(double fractional_precision) {
  // first, we find the range of energies for which we have data
  const int energies_observed = min_energy_state+1;

  // now we create two density of states vectors
  /* FIXME Can we work with ln_dos rather than dos to handle larger
     systems without underflow? */
  double *old_dos = new double[energies_observed];
  double *dos = new double[energies_observed];
  double dos_magnitude_squared;
  double dos_magnitude;

  // initialize a uniform density of states with unit norm
  const double dos_init = 1/sqrt(energies_observed);
  for (int i = 0; i < energies_observed; i++) dos[i] = dos_init;

  // now find the eigenvector of our transition matrix via the power iteration method
  /* fixme: I was thinking of implementing the inverse iteration method, but this
     method seems to converge quite nicely and very quickly, so it probably doesn't matter */
  bool done = false;
  int iters = 0;
  while(!done){
    iters++;
    for (int i = 0; i < energies_observed; i++) old_dos[i] = dos[i];

    // compute D_n = T*D_{n-1}
    dos_magnitude_squared = 0;
    for (int i = 0; i < energies_observed; i++){
      dos[i] = 0;
      for(int e = -biggest_energy_transition; e <= biggest_energy_transition; e++)
        /* FIXME the following looks buggy should be old_dos[i +/- e]; */
        dos[i] += transitions(i,e)*old_dos[i];
      dos_magnitude_squared += dos[i]*dos[i];
    }

    // normalize D_n
    dos_magnitude = sqrt(dos_magnitude_squared);
    for (int i = 0; i < energies_observed; i++) dos[i] /= dos_magnitude;

    // check whether D_n is close enough to D_{n-1} for us to quit
    done = true;
    for (int i = 0; i < energies_observed; i++){
      if (fabs(1-dos[i]/old_dos[i]) > fractional_precision){
        done = false;
        break;
      }
    }
  }

  // compute the weights w(E) = 1/D(E)
  for(int i = 0; i < energies_observed; i++)
    ln_energy_weights[i] = 1/dos[i];

  printf("Computed weights from transition matrix with dominant eigenvalue %g\n",
         dos_magnitude);
}

void sw_simulation::initialize_transitions(int max_iter) {
  /* We will simulate one single energy at a time. */
}

bool sw_simulation::printing_allowed(){
  const double time_skip = 1; // seconds
  const clock_t time_now = clock();
  if(time_now-last_print_time > time_skip*double(CLOCKS_PER_SEC)){
    last_print_time = time_now;
    return true;
  }
  else return false;
}
