import random
import numpy as np
import time
import pickle

##LÉVY DISTRIBUTION with exponent mu and cut-off lmax

rangemu_LevyDistrib = [round(1 + 0.1 * i, 2) for i in range(21)]
range_lmax = [112, 150,
              187,
              262,
              375,
              525,
              750,
              900,
              1050,
              1200]

a = dict()
ExpectedLength = dict()
for lmax in range_lmax:
    for mu in rangemu_LevyDistrib:
        s = 0
        for l in range(1, lmax): # Calcola la costante di normalizzazione per la distribuzione di Levy troncata
            s += 1.0 / (l**mu)
        a[(mu, lmax)] = s**(-1)  #########BEWARE NOW IT IS a[lmax,mu] and NOT a[n,mu]
    for mu in rangemu_LevyDistrib:
        s = 0
        # Correzione: il ciclo qui dovrebbe iterare su tutte le possibili lunghezze di passo 'l'
        # fino a 'lmax' per calcolare il valore atteso della lunghezza di passo.
        for l_step in range(1, lmax): # 'l_step' rappresenta le singole lunghezze di passo
            s += a[(mu, lmax)] * l_step / (l_step**mu) # E[L] = sum(l * P(l))
        ExpectedLength[(mu, lmax)] = s

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
        # In Python 3, / performs float division, so explicit float() is often not needed,
        # but 1.0 ensures consistency if mu was an int and l**mu resulted in an int.
        # Since mu is a float, l**mu will be a float, so division will be float.
        s += a[(mu, lmax)] / (l**mu)
    return l

def Levy_continuous(mu, lmax):
    """
    Generates a step length 'l' following a continuous Levy distribution
    with exponent 'mu' and a cut-off 'lmax'.
    """
    x = random.uniform(0, 1)
    if mu == 1:
        l = lmax ** x
    else:
        base = lmax**(1 - mu) - 1
        l = (x * base + 1) ** (1 / (1 - mu))
    return l

# Intermittent Levy search (detection only in-between steps)

def LevySearch3D(n, mu, lmax, D, TargetShape):
    """
    Simulates a 3D Levy search for a target shape within a cubic environment.
    Detection occurs only at the end of each Levy step.

    Args:
        n (int): Defines the size of the 3D cube (n^(1/3) is the side length).
        mu (float): Exponent of the Levy distribution.
        lmax (int): Cut-off length for Levy steps.
        D (float): Diameter/size parameter for the target shape.
        TargetShape (str): Type of target ('Ball', 'Line', 'Square').

    Returns:
        float: The total time (sum of Levy step lengths) until the target is found.
    """
    cube_side = n**(1/3)
    walker = np.array([random.uniform(0, cube_side),
                       random.uniform(0, cube_side),
                       random.uniform(0, cube_side)])
    center = np.array([cube_side * 0.5, cube_side * 0.5, cube_side * 0.5])

    time = 0
    found = False

    while not found:
        l = Levy(mu, lmax)

        # Coordinate sferiche per direzione 3D
        theta = random.uniform(0, 2 * np.pi)  # angolo azimutale [0, 2pi]
        phi = random.uniform(0, np.pi)        # angolo polare [0, pi]

        direction = np.array([
            np.sin(phi) * np.cos(theta),
            np.sin(phi) * np.sin(theta),
            np.cos(phi)
        ])

        walker += direction * l
        walker %= cube_side      # toro 3D

        time += l

        if TargetShape == 'Ball':
            # distanza euclidea 3D dal centro
            if np.linalg.norm(walker - center) <= 0.5 * D + 1:
                found = True

        elif TargetShape == 'Line':
            # Linea lungo asse x (ad esempio)
            if abs(walker[0] - center[0]) < 1:
                cond1 = abs(walker[1] - center[1]) < D * 0.5
                cond2 = np.linalg.norm(walker[[0, 1]] - np.array([center[0], center[1] + 0.5 * D])) < 1
                cond3 = np.linalg.norm(walker[[0, 1]] - np.array([center[0], cube_side * 0.5 - 0.5 * D])) < 1
                if cond1 or cond2 or cond3:
                    found = True

        elif TargetShape == 'Square':
            # ora un cubo (ipercubo 3D) centrato
            # Note: For a square, D/2*sqrt(3) might be for a specific relation to an inscribed sphere.
            # Assuming 'Square' actually means a cube here, given it's 3D and compares all 3 axes.
            half_side_threshold = D / (2 * np.sqrt(3)) + 1
            if (abs(center[0] - walker[0]) < half_side_threshold and
                abs(center[1] - walker[1]) < half_side_threshold and
                abs(center[2] - walker[2]) < half_side_threshold):
                found = True

    return time

def LevySearchProbaDetect3D(n, mu, lmax, D, TargetShape, p=1):
    """
    Simulates a 3D Levy search with probabilistic detection at each sub-step.

    Args:
        n (int): Defines the size of the 3D cube (n^(1/3) is the side length).
        mu (float): Exponent of the Levy distribution.
        lmax (int): Cut-off length for Levy steps.
        D (float): Diameter/size parameter for the target shape.
        TargetShape (str): Type of target ('Ball', 'Line', 'Square').
        p (float): Probability of detection if the walker is within the target.

    Returns:
        float: The total time (sum of individual steps) until the target is found.
    """
    time = 0
    cube_side = n ** (1/3)

    walker = np.array([random.uniform(0, cube_side),
                       random.uniform(0, cube_side),
                       random.uniform(0, cube_side)])

    center = np.array([cube_side * 0.5, cube_side * 0.5, cube_side * 0.5])

    found = False

    while not found:
        l = Levy(mu, lmax)
        # Direzione 3D in coordinate sferiche
        theta = random.uniform(0, 2 * np.pi)
        phi = random.uniform(0, np.pi)
        direction = np.array([
            np.sin(phi) * np.cos(theta),
            np.sin(phi) * np.sin(theta),
            np.cos(phi)
        ])

        # Move step by step and check for detection
        for _ in range(int(l)):
            walker += direction
            walker %= cube_side
            time += 1

            if TargetShape == 'Ball':
                if np.linalg.norm(walker - center) <= 0.5 * D + 1:
                    if random.uniform(0, 1) < p:
                        found = True
                        break # Break out of the inner loop if found

            elif TargetShape == 'Line':
                # Qui assumiamo linea lungo asse x come prima
                if abs(walker[0] - center[0]) < 1:
                    cond1 = abs(walker[1] - center[1]) < D * 0.5
                    cond2 = np.linalg.norm(walker[[0, 1]] - np.array([center[0], center[1] + 0.5 * D])) < 1
                    cond3 = np.linalg.norm(walker[[0, 1]] - np.array([center[0], cube_side * 0.5 - 0.5 * D])) < 1
                    if cond1 or cond2 or cond3:
                        if random.uniform(0, 1) < p:
                            found = True
                            break # Break out of the inner loop if found

            elif TargetShape == 'Square':
                half_side = D / (2 * np.sqrt(3)) + 1
                if (abs(center[0] - walker[0]) < half_side and
                    abs(center[1] - walker[1]) < half_side and
                    abs(center[2] - walker[2]) < half_side):
                    if random.uniform(0, 1) < p:
                        found = True
                        break # Break out of the inner loop if found

    return time

if __name__ == "__main__":
    ######COMPUTATIONS

    print("##FIGURE 3: n=300**3")
    TimesLevyProbaDetect = dict()

    start = time.time()
    for side in [32]: #[32, 300]
        n = side**3

        tests = 0
        while tests < 100:
            if tests == 10 or tests == 100: #evaluate how much time the task will require to complete
                # Changed print statement syntax for Python 3
                print(str(tests) + ' tests, total ' + str((time.time() - start) / 60.0) + ' minutes')
            tests += 1

            for lmax in range_lmax:
                for TargetShape in ['Ball', 'Line']: # Original code commented out ['Ball','Line','Square']
                    for D in [0, 1, 2, 4, 8, 16]:
                        for mu in rangemu_LevyDistrib:
                            for p in [1]: # [0, 0.1, 1]
                                if (n, mu, lmax, D, TargetShape, p) not in TimesLevyProbaDetect:
                                    TimesLevyProbaDetect[(n, mu, lmax, D, TargetShape, p)] = []
                                # Corrected function call from LevySearch to LevySearchProbaDetect3D
                                # and added the 'p' argument, as it's iterated over.
                                TimesLevyProbaDetect[(n, mu, lmax, D, TargetShape, p)].append(
                                    LevySearch3D(n, mu, lmax, D, TargetShape))
    end = time.time()
    # Changed print statement syntax for Python 3
    print((end - start) / 60.0)
"""
    ##FIGURE 4: influence of the cut-off lmax
    print("##FIGURE 4: influence of the cut-off lmax")

    start = time.time()
    for side in [300]:
        n = side**3
        l_medium = int(0.5 * side)
        tests = 0

        while tests < 1000:
            if tests == 10 or tests == 100: #evaluate how much time the task will require to complete
                print(str(tests) + ' tests, total ' + str((time.time() - start) / 60.0) + ' minutes')
            tests += 1

            for lmax_val in range_lmax:
                for TargetShape in ['Ball']: # Original code commented out ['Ball','Line','Square']
                    for D in [0, 2, 4, 8, 16]:
                        for mu in [2.0]:
                            # p is hardcoded to 0 in this new block, so using LevySearch3D
                            # which does not take 'p' as an argument.
                            if (n, mu, lmax_val, D, TargetShape, 0) not in TimesLevyProbaDetect:
                                TimesLevyProbaDetect[(n, mu, lmax_val, D, TargetShape, 0)] = []
                            print(TimesLevyProbaDetect)
                            TimesLevyProbaDetect[(n, mu, lmax_val, D, TargetShape, 0)].append(
                                LevySearch3D(n, mu, lmax_val, D, TargetShape))
    end = time.time()
    print((end - start) / 60.0)
"""
# save
filehandler=open('HittingTimeLevyDifferentSizeTargetShapeCutOffProbaDetect.obj','wb')
pickle.dump(TimesLevyProbaDetect,filehandler)
filehandler.close()