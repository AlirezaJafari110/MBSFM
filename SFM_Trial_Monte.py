#Trial
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML
import matplotlib as mpl
np.set_printoptions(threshold=np.inf)
import math
from Constants import *
from Intitialization import *
import Myfunctions_Explicit
def SFM_Trial_X(Alpha, Beta, Sigma, P_ped_y0):
    P_ped_x0=10.0
    P_ped=np.hstack((P_ped_x0,P_ped_y0))
    P_rob=np.array([0.0,0.0])
    V_rob=np.array([0.0, 0.0])
    P_ped_History=P_ped
    P_rob_History=P_rob
    Speed_History=[np.linalg.norm(V_rob), 0]
    PTTC_History=np.array([PTTC_max, 0.0])
    Theta_d=0.0
    Theta=-0.0#CHeck The Theta measurements
    z_dot=np.array([float(0),  Theta_d])
    z=np.array([v_0,  float(0)])
    u_rob=np.array([math.cos(Theta), math.sin(Theta)])
    V_ped=V_ped_0
    for i in range(n):
        Cn=[[0.0, -m*b*Theta_d],[m*b*Theta_d, 0.0]]
        t_i=i*dt
        f_boundary=[0.0, Myfunctions_Explicit.f_bnd(Width, P_rob, Alpha_bnd, Beta_bnd)]
        f_social_ped=Myfunctions_Explicit.f_soc_ped(P_ped, P_rob, Alpha_ped, Beta_ped, Lambda_ped, V_rob, V_ped)
        #Pedestrian SFM-Start
        F_ped=Myfunctions_Explicit.f_des_ped(np.array([-1.33, 0]), V_ped, tau_ped, m_ped)+[0.0, Myfunctions_Explicit.f_bnd(Width, P_ped, Alpha_bnd_ped, Beta_bnd)]+f_social_ped
        V_ped+=F_ped/m_ped*dt
        P_ped+=V_ped*dt
        #Override
        P_ped[0]=P_ped_x0
        P_ped[1]=P_ped_y0
        #Pedestrian SFM-End
        f_social=Myfunctions_Explicit.f_soc(P_ped, P_rob, Alpha, Beta, V_ped, V_rob, Lambda)
        S=np.linalg.norm(f_boundary)+np.linalg.norm(f_social)
        F=Myfunctions_Explicit.f_des_rob(V_max_vec, V_min_vec, V_rob, tau, m, S, Sigma)+f_boundary+f_social
        z_dot=Myfunctions_Explicit.Robot_Dynamics(Theta, F, Mn, Cn, Bn, z)
        z+=z_dot*dt
        v=z[0]
        w=z[1]
        Theta+=w*dt
        u_rob=np.array([math.cos(Theta), math.sin(Theta)])
        V_rob=v*u_rob
        Theta_d=w   
        P_rob+=V_rob*dt
        Speed=[np.linalg.norm(V_rob), t_i]
        Speed_History=np.vstack((Speed_History,Speed))
        P_rob_History=np.vstack((P_rob_History,P_rob))
        P_ped_History=np.vstack((P_ped_History,P_ped))
        PTTC_index=[Myfunctions_Explicit.PTTC(P_ped, P_rob, V_ped, V_rob, PTTC_max), t_i]
        PTTC_History=np.vstack((PTTC_History,PTTC_index))
    T_p=min(PTTC_History[:,0])
    I_v=Myfunctions_Explicit.Vel_Index(P_ped_History, P_rob_History, Speed_History, V_max, dt)
    I_t=1-math.exp(-T_p)
    if I_v==0:
        I_t=0
    #Q=I_v+I_t
    return I_v, I_t