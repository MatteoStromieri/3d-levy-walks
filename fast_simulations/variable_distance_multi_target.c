#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h> // Per strcmp e per gestire le stringhe
#include "func.h"

int main(int argc, char *argv[]) {
    const char *output_filename = "simulazioni_multiwalker_fixed_distance.csv";
    int num_max_trials = 100;
    double fixed_target_dist = -1.0;  // Default: target in posizione casuale

    if (argc < 2 || argc > 4) {
        printf("Utilizzo: %s <nome_file_csv> [numero_iterazioni] [distanza_fissa_target]\n", argv[0]);
        printf("  <nome_file_csv>: Nome del file CSV per i risultati.\n");
        printf("  [numero_iterazioni]: Numero di prove (default: %d).\n", num_max_trials);
        printf("  [distanza_fissa_target]: Distanza fissa dal centro per i target (opzionale).\n");
        return EXIT_FAILURE;
    }

    output_filename = argv[1];

    if (argc >= 3) {
        char *endptr;
        long trials = strtol(argv[2], &endptr, 10);
        if (*endptr != '\0' || trials <= 0 || trials > 2000000000) {
            fprintf(stderr, "Errore: 'numero_iterazioni' non valido. Deve essere un intero positivo.\n");
            return EXIT_FAILURE;
        }
        num_max_trials = (int)trials;
    }

    if (argc == 4) {
        char *endptr;
        double fixed_dist = strtod(argv[3], &endptr);
        if (*endptr != '\0' || fixed_dist < 0.0) {
            fprintf(stderr, "Errore: 'distanza_fissa_target' non valida. Deve essere >= 0.\n");
            return EXIT_FAILURE;
        }
        fixed_target_dist = fixed_dist;
    }    // Seed the random number generator
    srand((unsigned int)time(NULL));

    // Define ranges
    double rangemu_LevyDistrib[] = {2}; //{1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
    int len_rangemu_LevyDistrib = sizeof(rangemu_LevyDistrib) / sizeof(rangemu_LevyDistrib[0]);

    int range_diam[] = {0, 1, 2, 4, 8, 16, 32, 64};
    int len_range_diam = sizeof(range_diam) / sizeof(range_diam[0]);

    int range_side[] = {32, 64, 128, 256};
    int len_range_side = sizeof(range_side) / sizeof(range_side[0]);

    const char* list_shapes[] = {"Ball"};
    int len_list_shapes = sizeof(list_shapes) / sizeof(list_shapes[0]);

    int range_nwalkers[] = {1, 2, 4, 8, 16, 32, 64, 128};
    int len_range_nwalkers = sizeof(range_nwalkers) / sizeof(range_nwalkers[0]);

    int range_ntargets[] = {1, 2, 4, 8, 16};
    int len_range_ntargets = sizeof(range_ntargets) / sizeof(range_ntargets[0]);

    int range_distance[] = {2, 4, 8, 16, 32, 64, 128};
    int len_range_distance = sizeof(range_distance) / sizeof(range_distance[0]);

    // --- Pre-calculation of normalization constant and expected length ---
    printf("Pre-calcolando costanti di Lévy...\n");

    // Determine unique lmax values from range_side
    int unique_lmax_values[len_range_side]; // Max possible unique lmax values
    int current_unique_lmax_count = 0;
    for (int i = 0; i < len_range_side; ++i) {
        int lmax_val = range_side[i] / 2;
        int found = 0;
        for (int j = 0; j < current_unique_lmax_count; ++j) {
            if (unique_lmax_values[j] == lmax_val) {
                found = 1;
                break;
            }
        }
        if (!found) {
            unique_lmax_values[current_unique_lmax_count++] = lmax_val;
        }
    }

    for (int i = 0; i < len_rangemu_LevyDistrib; ++i) {
        double mu_val = rangemu_LevyDistrib[i];
        int mu_idx = get_mu_index(mu_val);
        if (mu_idx < 0 || mu_idx >= MAX_MU_INDEX) {
            fprintf(stderr, "Errore: mu_idx (%d) fuori dai limiti per mu_val %.1f. Aumentare MAX_MU_INDEX.\n", mu_idx, mu_val);
            return EXIT_FAILURE;
        }

        for (int k = 0; k < current_unique_lmax_count; ++k) {
            int lmax_val = unique_lmax_values[k];
            if (lmax_val < 0 || lmax_val > MAX_LMAX_VALUE) {
                fprintf(stderr, "Errore: lmax_val (%d) fuori dai limiti per MAX_LMAX_VALUE. Aumentare MAX_LMAX_VALUE.\n", lmax_val);
                return EXIT_FAILURE;
            }

            double sum_prob = 0;
            double sum_expected_length = 0;
            for (int l_step = 1; l_step <= lmax_val; ++l_step) {
                double prob_density = 1.0 / pow((double)l_step, mu_val);
                sum_prob += prob_density;
                sum_expected_length += l_step * prob_density;
            }
            if (sum_prob > 0) {
                a_coeffs[mu_idx][lmax_val] = 1.0 / sum_prob;
                expected_lengths[mu_idx][lmax_val] = a_coeffs[mu_idx][lmax_val] * sum_expected_length;
            } else {
                a_coeffs[mu_idx][lmax_val] = 0;
                expected_lengths[mu_idx][lmax_val] = 0;
            }
        }
    }
    printf("Pre-calcolo completato.\n");

    printf("## Starting Levy Search Simulations on 3D Torus\n");

    FILE *output_file = fopen(output_filename, "w");
    if (output_file == NULL) {
        perror("Error opening output file");
        return EXIT_FAILURE;
    }
    fprintf(output_file, "n_walkers,n_volume,mu,lmax,D,TargetShape,n_targets,fixed_target_dist,detection_time,probability\n");

    long start_total_time = time(NULL);

    long total_inner_iterations = (long)len_range_side * len_list_shapes * len_range_diam *
                                   len_range_ntargets * len_range_nwalkers * len_rangemu_LevyDistrib * len_range_distance;

    // Ora si usa num_max_trials preso da riga di comando o il default
    for (int current_trial = 0; current_trial < num_max_trials; ++current_trial) {
        printf("\n--- Overall Trial: %d/%d ---\n", current_trial + 1, num_max_trials);
        long elapsed_seconds = time(NULL) - start_total_time;
        printf("    %d trials completed, total elapsed time: %.2f minutes\n", current_trial + 1, (double)elapsed_seconds / 60.0);
        if (current_trial > 0) {
            double estimated_remaining_time = ((double)elapsed_seconds / (current_trial + 1)) * (num_max_trials - (current_trial + 1));
            printf("    Estimated time remaining: %.2f minutes\n", estimated_remaining_time / 60.0);
        }

        long pbar_counter = 0;

        for (int i_side = 0; i_side < len_range_side; ++i_side) {
            int side = range_side[i_side];
            double n_volume = pow(side, 3);
            int lmax_current = side / 2;

            for (int i_shape = 0; i_shape < len_list_shapes; ++i_shape) {
                const char* TargetShape = list_shapes[i_shape];
                for (int i_D = 0; i_D < len_range_diam; ++i_D) {
                    int D = range_diam[i_D];
                    for (int i_ntargets = 0; i_ntargets < len_range_ntargets; ++i_ntargets) {
                        int n_targets = range_ntargets[i_ntargets];
                        for (int i_nwalkers = 0; i_nwalkers < len_range_nwalkers; ++i_nwalkers) {
                            int n_walkers = range_nwalkers[i_nwalkers];
                            for (int i_mu = 0; i_mu < len_rangemu_LevyDistrib; ++i_mu) {
                                double mu = rangemu_LevyDistrib[i_mu];
                                for (int i_distance = 0; i_distance < len_range_distance; ++i_distance) {
                                    fixed_target_dist = range_distance[i_distance];
                                    double detection_time = LevySearch3D_MultiWalker(n_walkers, "nest", n_volume, mu, lmax_current,
                                                                                    D, TargetShape, n_targets, fixed_target_dist, 1.0);

                                    fprintf(output_file, "%d,%.0f,%.1f,%d,%d,%s,%d,%.1f,%.2f,1\n",
                                            n_walkers, n_volume, mu, lmax_current, D, TargetShape, n_targets, fixed_target_dist, detection_time);

                                    pbar_counter++;
                                    if (total_inner_iterations > 0 && pbar_counter % (total_inner_iterations / 100 + 1) == 0) {
                                        printf("\rTrial %d/%d Progress: %.2f%%", current_trial + 1, num_max_trials,
                                            (double)pbar_counter * 100.0 / total_inner_iterations);
                                        fflush(stdout);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        printf("\rTrial %d/%d Progress: 100.00%%\n", current_trial + 1, num_max_trials);
    }

    long end_total_time = time(NULL);
    printf("\nTotal simulation time: %.2f minutes\n", (double)(end_total_time - start_total_time) / 60.0);

    fclose(output_file);
    printf("Saving results to %s\n", output_filename);
    printf("Save complete.\n");

    return 0;
}