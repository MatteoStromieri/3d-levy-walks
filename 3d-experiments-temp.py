import random
import numpy as np
import time
import pickle

## LÉVY DISTRIBUTION with exponent mu and cut-off lmax

rangemu_LevyDistrib = [round(1 + 0.2 * i, 2) for i in range(21)]
range_lmax = [16, 112, 150,
              187,
              262,
              375,
              525,
              750,
              900,
              1050,
              1200]
range_diam = [0, 1, 2, 4, 8, 16, 32, 64]
range_side = [32, 64, 128, 256, 512]
list_shapes = ['Ball', 'Line', 'Disk']
range_nwalkers = [1, 2, 4, 8, 16]
range_ntargets = [1, 2, 4, 8, 16]

# --- Pre-calculation of normalization constant and expected length ---
a = dict()
ExpectedLength = dict()
for lmax_val in range_lmax:
    for mu_val in rangemu_LevyDistrib:
        sum_prob = 0
        sum_expected_length = 0
        # Include lmax_val in the loop if it represents the maximum valid length.
        # Used range(1, lmax_val + 1) to include lmax_val.
        for l_step in range(1, lmax_val + 1):
            prob_density = 1.0 / (l_step**mu_val)
            sum_prob += prob_density
            sum_expected_length += l_step * prob_density

        a[(mu_val, lmax_val)] = sum_prob**(-1)
        ExpectedLength[(mu_val, lmax_val)] = a[(mu_val, lmax_val)] * sum_expected_length

# --- Function to generate a Levy step ---
def Levy(mu, lmax):
    """
    Generates a step length 'l' following a discretized Levy distribution
    with exponent 'mu' and a cut-off 'lmax'.
    """
    x = random.uniform(0, 1)
    s = 0
    l = 0
    while s < x:
        l += 1
        # If l exceeds lmax, it should be handled based on the definition of the truncated distribution.
        # Here, it's assumed that the distribution is only defined up to lmax.
        if l > lmax: # Added a check to prevent infinite loops if x is close to 1 and lmax is small
            return lmax
        s += a[(mu, lmax)] / (l**mu)
    return l

# --- Functions for calculating toroidal distance ---
def toroidal_distance_squared(p1, p2, side_length):
    """
    Calculates the squared minimum distance between two points on a 3D torus.
    p1, p2: numpy arrays [x, y, z]
    side_length: length of the side of the toroidal cube
    """
    delta = np.abs(p1 - p2)
    # Apply periodicity condition: min(direct_dist, side_length - direct_dist)
    delta = np.minimum(delta, side_length - delta)
    return np.sum(delta**2)

def toroidal_distance(p1, p2, side_length):
    """
    Calculates the minimum distance between two points on a 3D torus.
    """
    return np.sqrt(toroidal_distance_squared(p1, p2, side_length))

# The distance_point_to_segment_euclidean function is for Euclidean space.
# For a torus, the concept is more complex (periodic images).
# For axis-aligned or point targets, we will use the toroidal_distance functions.
# If TargetShape == 'Line' is a generic line, this function alone is not sufficient.
# I've adapted the Line detection assuming it's centered and axis-aligned.
def distance_point_to_segment_euclidean(point, p1, p2):
    """
    Calculates the shortest distance from a point to a line segment in 3D (Euclidean space).
    point: numpy array [x, y, z] of the walker
    p1: numpy array [x, y, z] of the segment's start point
    p2: numpy array [x, y, z] of the segment's end point
    """
    if np.array_equal(p1, p2): # Degenerate case: segment is a point
        return np.linalg.norm(point - p1)

    v = p2 - p1 # Vector of segment P1P2
    w = point - p1 # Vector from P1 to the point

    t = np.dot(w, v) / np.dot(v, v)

    if t < 0.0: # Projection falls before P1
        distance = np.linalg.norm(point - p1)
    elif t > 1.0: # Projection falls after P2
        distance = np.linalg.norm(point - p2)
    else: # Projection falls within the segment
        projection = p1 + t * v
        distance = np.linalg.norm(point - projection)
    return distance

# --- Multi-walker and multi-target search function ---
def LevySearch3D_MultiWalker(n_walkers, initialization, n, mu, lmax, D, TargetShape, num_targets_to_generate):
    cube_side = n**(1/3)
    
    if initialization == 'independent':
        walkers = np.array([[random.uniform(0, cube_side),
                             random.uniform(0, cube_side),
                             random.uniform(0, cube_side)] for _ in range(n_walkers)])
    elif initialization == 'nest': # All walkers start from the same random point
        random_point = np.array([random.uniform(0, cube_side),
                                 random.uniform(0, cube_side),
                                 random.uniform(0, cube_side)])
        walkers = np.array([random_point.copy() for _ in range(n_walkers)]) # Use .copy() to avoid references to the same array
    
    target_positions = np.array([[random.uniform(0, cube_side),
                                  random.uniform(0, cube_side),
                                  random.uniform(0, cube_side)] for _ in range(num_targets_to_generate)])

    walker_times = np.zeros(n_walkers) # Accumulated time for each walker
    discovery_times = np.full(n_walkers, float('inf')) # Discovery time for each walker, initialized to infinity

    any_walker_found_target = False

    while True:
        min_overall_discovery_time = np.min(discovery_times) # The lowest discovery time among all walkers so far

        # Termination criterion: If at least one walker has found the target
        # AND all other walkers that haven't found the target yet have
        # a cumulative time already greater than the minimum discovery time.
        if any_walker_found_target:
            all_remaining_walkers_past_min_time = True
            for i in range(n_walkers):
                if discovery_times[i] == float('inf') and walker_times[i] < min_overall_discovery_time:
                    all_remaining_walkers_past_min_time = False
                    break
            
            if all_remaining_walkers_past_min_time:
                return min_overall_discovery_time # Returns the time of the first detection

        for i in range(n_walkers):
            # If this walker has already found the target and its discovery time is recorded,
            # or if its cumulative time has already exceeded the minimum discovery time,
            # there's no point in making it walk further in this loop to find the *first* target.
            if discovery_times[i] != float('inf') and walker_times[i] >= min_overall_discovery_time:
                continue # This walker has already contributed to the minimum time or is too slow

            l = Levy(mu, lmax)
            walker_times[i] += l
            
            theta = random.uniform(0, 2 * np.pi)
            phi = random.uniform(0, np.pi)
            direction = np.array([
                np.sin(phi) * np.cos(theta),
                np.sin(phi) * np.sin(theta),
                np.cos(phi)
            ])
            
            walkers[i] += direction * l
            walkers[i] %= cube_side # Apply toroidal boundary conditions

            found_this_walker_in_this_step = False
            
            # Target detection check for each target
            for target_center in target_positions:
                if TargetShape == 'Ball':
                    # Use toroidal distance
                    if toroidal_distance(walkers[i], target_center, cube_side) <= 0.5 * D + 1:
                        found_this_walker_in_this_step = True
                        break
                
                elif TargetShape == 'Line':
                    # Assume a line centered along the Y-axis for simpler toroidal distance.
                    # D is the length of the line.
                    # Calculate the toroidal distance in X and Z dimensions relative to the line's center.
                    dx_torus = np.minimum(np.abs(walkers[i][0] - target_center[0]), cube_side - np.abs(walkers[i][0] - target_center[0]))
                    dz_torus = np.minimum(np.abs(walkers[i][2] - target_center[2]), cube_side - np.abs(walkers[i][2] - target_center[2]))
                    
                    # For the Y-coordinate, consider the Euclidean distance from the line segment.
                    # If the line is very long and could wrap around the torus, this part would be more complex.
                    # For a "local" (non-wrapping) line, Euclidean distance along the line's axis is fine.
                    # Here, we use Euclidean distance for the y-part relative to the segment.
                    p1_line_y = target_center[1] - D / 2.0
                    p2_line_y = target_center[1] + D / 2.0
                    
                    # Find the closest point along the Y-projection of the segment
                    closest_y = max(p1_line_y, min(walkers[i][1], p2_line_y))
                    dy = walkers[i][1] - closest_y

                    # The overall distance is the square root of the sum of squared minimum distances
                    dist_to_line_torus_squared = dx_torus**2 + dy**2 + dz_torus**2
                    
                    if np.sqrt(dist_to_line_torus_squared) <= 1: # Tolerance of 1
                        found_this_walker_in_this_step = True
                        break

                elif TargetShape == 'Square': # Intended as an axis-aligned Cube
                    # D is interpreted as the side length of the cube.
                    half_side_target = D / 2.0 + 1 # Your effective detection radius (D/2 + tolerance)
                    
                    # Calculate the toroidal difference for each axis
                    delta_x = np.minimum(np.abs(target_center[0] - walkers[i][0]), cube_side - np.abs(target_center[0] - walkers[i][0]))
                    delta_y = np.minimum(np.abs(target_center[1] - walkers[i][1]), cube_side - np.abs(target_center[1] - walkers[i][1]))
                    delta_z = np.minimum(np.abs(target_center[2] - walkers[i][2]), cube_side - np.abs(target_center[2] - walkers[i][2]))

                    # Check if the walker is inside the toroidal "detection cube"
                    if (delta_x < half_side_target and
                        delta_y < half_side_target and
                        delta_z < half_side_target):
                        found_this_walker_in_this_step = True
                        break
                
                elif TargetShape == 'Disk':
                    # Toroidal distance for XY coordinates, toroidal distance for Z.
                    # Assumes a disk parallel to the XY plane, centered at target_center.
                    
                    # Distance in the XY plane (toroidal)
                    dist_xy_squared = toroidal_distance_squared(walkers[i][:2], target_center[:2], cube_side)
                    dist_xy = np.sqrt(dist_xy_squared)

                    # Distance along the Z axis (toroidal)
                    dist_z = np.minimum(np.abs(walkers[i][2] - target_center[2]), cube_side - np.abs(walkers[i][2] - target_center[2]))

                    if dist_xy <= 0.5 * D + 1 and dist_z <= 1: # D is the disk diameter, 1 is Z tolerance
                        found_this_walker_in_this_step = True
                        break
            
            if found_this_walker_in_this_step:
                if discovery_times[i] == float('inf'): # Only record the first discovery for this walker
                    discovery_times[i] = walker_times[i]
                    any_walker_found_target = True


if __name__ == "__main__":
    ###### COMPUTATIONS

    print("## Starting Levy Search Simulations on 3D Torus")
    TimesLevyProbaDetect = dict()

    start_total_time = time.time()
    
    num_trials_per_config = 100 # Number of simulations to run for each parameter combination

    # Reorganizing loops to run multiple trials for each configuration
    for side in range_side:
        n = side**3 # Volume of the toroidal cube
        lmax_current = side / 2 # Lmax proportional to the cube side

        for TargetShape in list_shapes:
            for D in range_diam:
                for n_targets in range_ntargets:
                    for n_walkers in range_nwalkers:
                        for mu in rangemu_LevyDistrib:
                            # 'p' was removed from parameters as it was not used
                            
                            config_key = (n_walkers, n, mu, lmax_current, D, TargetShape, n_targets)
                            
                            if config_key not in TimesLevyProbaDetect:
                                TimesLevyProbaDetect[config_key] = []

                            print(f"\nSimulating Config: Side={side}, Shape={TargetShape}, D={D}, Targets={n_targets}, Walkers={n_walkers}, Mu={mu}")
                            
                            for trial in range(num_trials_per_config):
                                # Using 'nest' as initialization, as specified in your original code
                                detection_time = LevySearch3D_MultiWalker(n_walkers, 'nest', n, mu, lmax_current, D, TargetShape, n_targets)
                                TimesLevyProbaDetect[config_key].append(detection_time)
                                
                                # Print occasional progress
                                if (trial + 1) % 10 == 0 or trial == num_trials_per_config - 1:
                                    print(f"  Trial {trial + 1}/{num_trials_per_config} completed. Last time: {detection_time:.2f}")

    end_total_time = time.time()
    print(f"\nTotal simulation time: {(end_total_time - start_total_time) / 60.0:.2f} minutes")

# save
print("Saving results...")
filehandler = open('simulazioni_multiwalker_multitarget_toroidal.obj','wb')
pickle.dump(TimesLevyProbaDetect,filehandler)
filehandler.close()
print("Save complete.")