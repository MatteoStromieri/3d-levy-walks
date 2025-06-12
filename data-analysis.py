import matplotlib.pyplot as plt
import numpy as np 
import pickle
import os

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
# load data
filehandler=open('data/HittingTimeLevyDifferentSizeTargetShapeCutOffProbaDetect-temp.obj','rb')
TimesLevyProbaDetect=pickle.load(filehandler)
filehandler.close()
# setup
output_dir = 'Levy_Plots'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

"""
p=1
for TargetShape in ['Ball', 'Line']:
    plt.figure()
    CauchyRatio=dict()
    Opt = dict()
    rangeD=[0,1,2,4,8, 16]
    for side in [32]:
        n = 32**3
        lmax=int(side/2)
        for D in rangeD:
            m=0
            for mu in rangemu_LevyDistrib:
                #Opt[(n,TargetShape,mu,lmax,D)]=np.average(TimesLevyProbaDetect[(n,mu,lmax,D,TargetShape,p)])
                cauchy = np.average(TimesLevyProbaDetect[(n,2,lmax,D,TargetShape,p)])
                CauchyRatio[(n,TargetShape,mu,lmax,D)]=np.average(TimesLevyProbaDetect[(n,mu,lmax,D,TargetShape,p)]) / cauchy
            plt.plot(rangemu_LevyDistrib,[CauchyRatio[(n, TargetShape, mu, lmax, D)] for mu in rangemu_LevyDistrib],'-', label='D=%s'%D)
            
        plt.plot([1.0,3.0],[1,1], '-.k')
        plt.legend()

        plt.xlabel('mu')
        plt.ylabel('t_detect(X^mu)')
        # Salva la figura
        # Costruisci un nome file significativo
        filename = os.path.join(output_dir, f'Levy_Performance_n{n}_cutoff{lmax}_{TargetShape}_D.png')
        plt.savefig(filename)
        plt.close() # Chiudi la figura per liberare memoria e non sovrapporre i plot
    """

p=1
mu = 2.0
side = 32
plt.figure()
for TargetShape in ['Ball', 'Line']:
    if TargetShape == 'Line':
        filehandler=open('data/simulazioni_linea.obj','rb')
        TimesLevyProbaDetect=pickle.load(filehandler)
        filehandler.close()
        
    time = dict()
    rangeD=[1,2,4,8, 16]
    n = 32**3
    lmax=int(side/2)
        
    plt.plot(rangeD,[np.average(TimesLevyProbaDetect[(n, mu, lmax, D, TargetShape, p)]) for D in rangeD],'-', label='%s'%TargetShape)

            
plt.plot([1.0,3.0],[1,1], '-.k')
plt.legend()

plt.xlabel('D')
plt.ylabel('t_detect(X^2)')
# Salva la figura
# Costruisci un nome file significativo
filename = os.path.join(output_dir, f'Ball_vs_Line.png')
plt.savefig(filename)
plt.close() # Chiudi la figura per liberare memoria e non sovrapporre i plot