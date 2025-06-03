import numpy as np
import math
from Constants import *
#Intitialization
#Simulation Constants
P_ped_y0=0.5
PTTC_max=5.0
dt=0.03
t_end=15
n=int(t_end/dt)+1
t=np.linspace(0,10,n)
v_0=0.0
V_ped_0=np.array([-1.33, 0.0])
u_V_des=np.array([1, 0])
V_max_vec=V_max*u_V_des
V_min_vec=V_min*u_V_des
Mn=[[m, 0], [0, m*pow(b,2)+I]]
Bn=[[1, 0, 0], [0, b, 1]]