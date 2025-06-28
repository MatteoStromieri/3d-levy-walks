#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h> // Per strcmp e per gestire le stringhe
#include "func.h"

double a_coeffs[MAX_MU_INDEX][MAX_LMAX_VALUE + 1];
double expected_lengths[MAX_MU_INDEX][MAX_LMAX_VALUE + 1];

int get_mu_index(double mu_val) {
    // This needs to map 1.0 -> 0, 1.2 -> 1, ..., 3.0 -> 10
    // (mu_val - 1.0) / 0.2
    return (int)round((mu_val - 1.0) / 0.2);
}

// --- Function to calculate toroidal distance squared ---
double toroidal_distance_squared(double p1[3], double p2[3], double side_length) {
    double delta[3];
    delta[0] = fabs(p1[0] - p2[0]);
    delta[1] = fabs(p1[1] - p2[1]);
    delta[2] = fabs(p1[2] - p2[2]);

    delta[0] = fmin(delta[0], side_length - delta[0]);
    delta[1] = fmin(delta[1], side_length - delta[1]);
    delta[2] = fmin(delta[2], side_length - delta[2]);

    return delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
}

// --- Function to calculate toroidal distance squared for 2D ---
// Aggiunta per risolvere il warning stringop-overflow
double toroidal_distance_squared_2D(double p1[2], double p2[2], double side_length) {
    double delta[2];
    delta[0] = fabs(p1[0] - p2[0]);
    delta[1] = fabs(p1[1] - p2[1]);

    delta[0] = fmin(delta[0], side_length - delta[0]);
    delta[1] = fmin(delta[1], side_length - delta[1]);

    return delta[0]*delta[0] + delta[1]*delta[1];
}


// --- Function to calculate toroidal distance ---
double toroidal_distance(double p1[3], double p2[3], double side_length) {
    return sqrt(toroidal_distance_squared(p1, p2, side_length));
}

// --- Function to generate a Levy step ---
int Levy(double mu, int lmax) {
    double x = (double)rand() / RAND_MAX; // random.uniform(0, 1)
    double s = 0;
    int l = 0;
    int mu_idx = get_mu_index(mu); // Get the index for mu

    while (s < x) {
        l++;
        if (l > lmax) {
            return lmax; // Prevent infinite loop
        }
        // Access pre-calculated 'a' coefficient
        // Assicurati che lmax sia un indice valido per a_coeffs e non vada fuori bounds
        if (lmax >= 0 && lmax <= MAX_LMAX_VALUE) { // Aggiunto un controllo di sicurezza per lmax
             s += a_coeffs[mu_idx][lmax] / pow((double)l, mu);
        } else {
            // Gestisci l'errore o usa un valore predefinito
            // Questo scenario dovrebbe essere evitato con MAX_LMAX_VALUE corretto
            fprintf(stderr, "Errore: lmax_val (%d) supera MAX_LMAX_VALUE (%d) in Levy.\n", lmax, MAX_LMAX_VALUE);
            return lmax; // O un altro valore sensato per evitare crash
        }
    }
    return l;
}



double LevySearch3D_MultiWalker(int n_walkers, const char* initialization, double n_volume, double mu, int lmax,
                                int D, const char* TargetShape, int num_targets_to_generate,
                                double target_distance_from_origin, double p) {
    double cube_side = cbrt(n_volume); // Equivalent to n**(1/3)

    // Allocate memory for walkers
    double (*walkers)[3] = (double (*)[3])malloc(n_walkers * sizeof(double[3]));
    if (walkers == NULL) {
        fprintf(stderr, "Memory allocation failed for walkers\n");
        exit(EXIT_FAILURE);
    }

    // Initialize walkers
    if (strcmp(initialization, "independent") == 0) {
        for (int i = 0; i < n_walkers; ++i) {
            walkers[i][0] = (double)rand() / RAND_MAX * cube_side;
            walkers[i][1] = (double)rand() / RAND_MAX * cube_side;
            walkers[i][2] = (double)rand() / RAND_MAX * cube_side;
        }
    } else if (strcmp(initialization, "nest") == 0) {
        double random_point[3];
        random_point[0] = (double)rand() / RAND_MAX * cube_side;
        random_point[1] = (double)rand() / RAND_MAX * cube_side;
        random_point[2] = (double)rand() / RAND_MAX * cube_side;
        for (int i = 0; i < n_walkers; ++i) {
            walkers[i][0] = random_point[0];
            walkers[i][1] = random_point[1];
            walkers[i][2] = random_point[2];
        }
    }

    // Allocate memory for target positions
    double (*target_positions)[3] = (double (*)[3])malloc(num_targets_to_generate * sizeof(double[3]));
    if (target_positions == NULL) {
        fprintf(stderr, "Memory allocation failed for target_positions\n");
        exit(EXIT_FAILURE);
    }

    // Initialize target positions
    if (target_distance_from_origin >= 0) { // Fixed distance
        for (int i = 0; i < num_targets_to_generate; ++i) {
            double theta = (double)rand() / RAND_MAX * 2 * M_PI;
            double phi = (double)rand() / RAND_MAX * M_PI;

            double x = target_distance_from_origin * sin(phi) * cos(theta);
            double y = target_distance_from_origin * sin(phi) * sin(theta);
            double z = target_distance_from_origin * cos(phi);

            target_positions[i][0] = fmod(x + cube_side / 2.0, cube_side);
            target_positions[i][1] = fmod(y + cube_side / 2.0, cube_side);
            target_positions[i][2] = fmod(z + cube_side / 2.0, cube_side);
            if (target_positions[i][0] < 0) target_positions[i][0] += cube_side;
            if (target_positions[i][1] < 0) target_positions[i][1] += cube_side;
            if (target_positions[i][2] < 0) target_positions[i][2] += cube_side;
        }
    } else { // Random placement
        for (int i = 0; i < num_targets_to_generate; ++i) {
            target_positions[i][0] = (double)rand() / RAND_MAX * cube_side;
            target_positions[i][1] = (double)rand() / RAND_MAX * cube_side;
            target_positions[i][2] = (double)rand() / RAND_MAX * cube_side;
        }
    }

    double* walker_times = (double*)calloc(n_walkers, sizeof(double)); // Initialize to 0.0
    double* discovery_times = (double*)malloc(n_walkers * sizeof(double));
    if (walker_times == NULL || discovery_times == NULL) {
        fprintf(stderr, "Memory allocation failed for times arrays\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n_walkers; ++i) {
        discovery_times[i] = HUGE_VAL; // Initialize to infinity
    }

    int any_walker_found_target = 0;

    while (1) {
        double min_overall_discovery_time = HUGE_VAL;
        for (int i = 0; i < n_walkers; ++i) {
            if (discovery_times[i] < min_overall_discovery_time) {
                min_overall_discovery_time = discovery_times[i];
            }
        }

        if (any_walker_found_target) {
            int all_remaining_walkers_past_min_time = 1;
            for (int i = 0; i < n_walkers; ++i) {
                if (discovery_times[i] == HUGE_VAL && walker_times[i] < min_overall_discovery_time) {
                    all_remaining_walkers_past_min_time = 0;
                    break;
                }
            }
            if (all_remaining_walkers_past_min_time) {
                free(walkers);
                free(target_positions);
                free(walker_times);
                free(discovery_times);
                return min_overall_discovery_time;
            }
        }

        for (int i = 0; i < n_walkers; ++i) {
            if (discovery_times[i] != HUGE_VAL && walker_times[i] >= min_overall_discovery_time) {
                continue;
            }

            int l = Levy(mu, lmax);
            walker_times[i] += l;

            double theta = (double)rand() / RAND_MAX * 2 * M_PI;
            double phi = (double)rand() / RAND_MAX * M_PI;
            double direction[3];
            direction[0] = sin(phi) * cos(theta);
            direction[1] = sin(phi) * sin(theta);
            direction[2] = cos(phi);

            walkers[i][0] += direction[0] * l;
            walkers[i][1] += direction[1] * l;
            walkers[i][2] += direction[2] * l;

            // Apply toroidal boundary conditions
            walkers[i][0] = fmod(walkers[i][0], cube_side);
            walkers[i][1] = fmod(walkers[i][1], cube_side);
            walkers[i][2] = fmod(walkers[i][2], cube_side);
            if (walkers[i][0] < 0) walkers[i][0] += cube_side;
            if (walkers[i][1] < 0) walkers[i][1] += cube_side;
            if (walkers[i][2] < 0) walkers[i][2] += cube_side;

            int found_this_walker_in_this_step = 0;

            for (int j = 0; j < num_targets_to_generate; ++j) {
                double *target_center = target_positions[j];

                int inside_target = 0;

                if (strcmp(TargetShape, "Ball") == 0) {
                    if (toroidal_distance(walkers[i], target_center, cube_side) <= 0.5 * D + 1) {
                        inside_target = 1;
                    }
                } else if (strcmp(TargetShape, "Line") == 0) {
                    double dx_torus = fmin(fabs(walkers[i][0] - target_center[0]), cube_side - fabs(walkers[i][0] - target_center[0]));
                    double dz_torus = fmin(fabs(walkers[i][2] - target_center[2]), cube_side - fabs(walkers[i][2] - target_center[2]));

                    double p1_line_y = target_center[1] - D / 2.0;
                    double p2_line_y = target_center[1] + D / 2.0;

                    double closest_y = fmax(p1_line_y, fmin(walkers[i][1], p2_line_y));
                    double dy = walkers[i][1] - closest_y;

                    double dist_to_line_torus_squared = dx_torus*dx_torus + dy*dy + dz_torus*dz_torus;

                    if (sqrt(dist_to_line_torus_squared) <= 1) {
                        inside_target = 1;
                    }

                } else if (strcmp(TargetShape, "Square") == 0) { // Cube
                    double half_side_target = D / 2.0 + 1;

                    double delta_x = fmin(fabs(target_center[0] - walkers[i][0]), cube_side - fabs(target_center[0] - walkers[i][0]));
                    double delta_y = fmin(fabs(target_center[1] - walkers[i][1]), cube_side - fabs(target_center[1] - walkers[i][1]));
                    double delta_z = fmin(fabs(target_center[2] - walkers[i][2]), cube_side - fabs(target_center[2] - walkers[i][2]));

                    if (delta_x < half_side_target &&
                        delta_y < half_side_target &&
                        delta_z < half_side_target) {
                        inside_target = 1;
                    }
                } else if (strcmp(TargetShape, "Disk") == 0) {
                    double walker_xy[2] = {walkers[i][0], walkers[i][1]};
                    double target_xy[2] = {target_center[0], target_center[1]};
                    double dist_xy_squared = toroidal_distance_squared_2D(walker_xy, target_xy, cube_side);
                    double dist_xy = sqrt(dist_xy_squared);

                    double dist_z = fmin(fabs(walkers[i][2] - target_center[2]), cube_side - fabs(walkers[i][2] - target_center[2]));

                    if (dist_xy <= 0.5 * D + 1 && dist_z <= 1) {
                        inside_target = 1;
                    }
                }

                if (inside_target) {
                    double r = (double)rand() / RAND_MAX;
                    if (r <= p) {
                        found_this_walker_in_this_step = 1;
                        break;
                    }
                }
            }

            if (found_this_walker_in_this_step) {
                if (discovery_times[i] == HUGE_VAL) {
                    discovery_times[i] = walker_times[i];
                    any_walker_found_target = 1;
                }
            }
        }
    }
}
