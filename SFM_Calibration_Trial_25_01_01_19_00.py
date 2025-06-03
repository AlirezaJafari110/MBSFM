import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML
import matplotlib as mpl
from Constants import *
from Intitialization import *
np.set_printoptions(threshold=np.inf)
from SFM_Trial_Monte import SFM_Trial_X as STX
import time
#Initialization
t_000=time.time()
W=1.0
N_a=11
N_b=21
N_s=11
N=1
Alpha_Set=np.linspace(50, 300, N_a)
Beta_Set=np.linspace(0.5, 2.5, N_b)
Sigma_Set=np.linspace(50, 150, N_s)
shape = (N_a, N_b, N_s, N)  
Q = np.zeros(shape)
I_t = np.zeros(shape)
I_v = np.zeros(shape)
Result=np.array([0.0, 0.0, 0.0, 0.0, 0.0], dtype=np.float32)
print(Result)
Max=-1.0
counter=0
rows=[]
for i_a in range(N_a):
    Alpha=Alpha_Set[i_a]
    for i_b in range(N_b):
        Beta=Beta_Set[i_b]
        for i_s in range(N_s):
            Sigma=Sigma_Set[i_s]    
            #P_ped_y0_Set=np.linspace(0.45, 0.55, N+2)[1:-1]
            for i in range(N):
                    #P_ped_y0=P_ped_y0_Set[i]
                    print(Alpha, Beta, Sigma)
                    I_v[i_a,i_b, i_s, i], I_t[i_a,i_b, i_s, i]=STX(Alpha, Beta, Sigma, P_ped_y0)
                    Q[i_a,i_b, i_s, i]=I_v[i_a,i_b, i_s, i]+W*I_t[i_a,i_b, i_s, i]
                    if Q[i_a,i_b, i_s, i]!=0 and Q[i_a,i_b, i_s, i]>Max and I_v[i_a,i_b, i_s, i]>0.5 and I_t[i_a,i_b, i_s, i]>0.5:
                        Max=Q[i_a,i_b, i_s, i]
                        counter+=1
                        Result=[Alpha, Beta, Sigma]
                        New_row=[Max, I_v[i_a,i_b, i_s, i], I_t[i_a,i_b, i_s, i], Alpha, Beta, Sigma]
                        rows.append(New_row)
                        print('Max=', Max, 'I_v=', I_v[i_a,i_b, i_s, i], 'I_t=',I_t[i_a,i_b, i_s, i], 'T_p=', -math.log(1-I_t[i_a,i_b, i_s, i]))
                        print(Result[0], Result[1], Result[2])
print('Result: ','Alpha=',Result[0], 'Beta=',Result[1],'Sigma=',Result[2], 'Max=', Max)
print(counter)
Result_array = np.array(rows)
duration=time.time()-t_000
print('Result Array:\n ', Result_array )
#print(duration)