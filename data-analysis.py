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
filehandler=open('data/simulazioni_multiwalker_multitarget.obj','rb')
TimesLevyProbaDetect=pickle.load(filehandler)
filehandler.close()
# setup
output_dir = 'Levy_Plots'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

filehandler=open('data/simulazioni_linea.obj','rb')
TimesLevyProbaDetect=pickle.load(filehandler)
filehandler.close()


# plot diametro d contro k figure di diametro D/k
p=1
plt.figure()
rangeD=[0,1,2,4,8,16]
n_targets = [1, 2, 4, 8]
side = 32
n = 32**3
lmax=int(side/2)
n_walkers = 1
D = 8
mu=2
TargetShape = 'Ball'
val = []
for i,k in enumerate(n_targets):
    val[i] = np.average(TimesLevyProbaDetect[(n_walkers, n, mu, lmax, D/k, TargetShape, k, p)])
    
plt.plot(n_targets, val,'-', label='D=%s'%D)
            
plt.plot([1.0,3.0],[1,1], '-.k')
plt.legend()

plt.xlabel('#targets')
plt.ylabel('t_detect')
# Salva la figura
# Costruisci un nome file significativo
filename = os.path.join(output_dir, f'multi_target.png')
plt.savefig(filename)
plt.close() # Chiudi la figura per liberare memoria e non sovrapporre i plot
    
# multi-agent vs one agent
p=1
plt.figure()
n_targets = 1
side = 32
n = 32**3
lmax=int(side/2)
n_walkers = [1, 2, 4, 8, 16]
D = 4
mu=2
TargetShape = 'Ball'
val = []
for i,k in enumerate(n_walkers):
    val[i] = np.average(TimesLevyProbaDetect[(k, n, mu, lmax, D, TargetShape, n_targets, p)])
    
plt.plot(n_targets, val,'-', label='D=%s'%D)
            
plt.plot([1.0,3.0],[1,1], '-.k')
plt.legend()

plt.xlabel('#agents')
plt.ylabel('t_detect')
# Salva la figura
# Costruisci un nome file significativo
filename = os.path.join(output_dir, f'multi_agent.png')
plt.savefig(filename)
plt.close() # Chiudi la figura per liberare memoria e non sovrapporre i plot


# line vs disk vs ball
p=1
plt.figure()
rangeD=[0,1,2,4,8,16]
n_targets = 1 #[1, 2, 4, 8]
side = 32
n = 32**3
lmax=int(side/2)
n_walkers = 1
mu=2
val = []
for TargetShape in ['Line', 'Disk', 'Ball']:
    for i,D in enumerate(rangeD):
        val[i] = np.average(TimesLevyProbaDetect[(n_walkers, n, mu, lmax, D, TargetShape, n_targets, p)])
    plt.plot(D, val,'-', label=TargetShape)
            
plt.plot([1.0,3.0],[1,1], '-.k')
plt.legend()

plt.xlabel('D')
plt.ylabel('t_detect')
# Salva la figura
# Costruisci un nome file significativo
filename = os.path.join(output_dir, f'line_disk_ball.png')
plt.savefig(filename)
plt.close() # Chiudi la figura per liberare memoria e non sovrapporre i plot