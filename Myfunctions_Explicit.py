
import numpy as np
import math

def f_des_ped(V_des, V, tau, m):
    f_des=m*(V_des-V)/tau
    return f_des

def f_des_rob(V_max_vec, V_min_vec, V, tau, m, S, sigma):
    V_des=math.exp(-S/sigma)*(V_max_vec-V_min_vec)+V_min_vec
    f_des=m*(V_des-V)/tau
    return f_des

def f_bnd(Width, P_rob, Alpha_bnd, Beta_bnd):
    d_ij_up=abs(Width/2-P_rob[1])
    d_ij_low=abs(-Width/2-P_rob[1])
    f_bnd_upward=Alpha_bnd*np.exp(-d_ij_low/Beta_bnd)
    f_bnd_downward=-Alpha_bnd*np.exp(-d_ij_up/Beta_bnd)
    f_bnd=f_bnd_upward+f_bnd_downward
    return f_bnd

def f_soc(P_ped, P_rob, Alpha, Beta, V_ped, V_rob, Lambda):
    dij=P_rob-P_ped
    norm_dij=np.linalg.norm(dij)
    u_dij=dij/norm_dij
    f_soc=Alpha*np.exp(-norm_dij/Beta)*u_dij*Gamma(P_ped, P_rob, V_ped, V_rob, Lambda)
    return f_soc

def f_soc_ped(P_ped, P_rob, Alpha_ped, Beta_ped, Lambda_ped, V_rob, V_ped):
    dij=P_ped-P_rob
    norm_dij=np.linalg.norm(dij)
    u_dij=dij/norm_dij
    f_soc_ped=Alpha_ped*np.exp(-norm_dij/Beta_ped)*u_dij*Gamma(P_rob, P_ped, V_rob, V_ped, Lambda_ped)
    return f_soc_ped

def Robot_Dynamics(Theta, F, Mn, Cn, Bn, z):
    u_rob=np.array([math.cos(Theta), math.sin(Theta)])
    u_n_rob=np.array([-math.sin(Theta), math.cos(Theta)])
    Fm=np.dot(F, u_rob)
    Fn=np.dot(F, u_n_rob)
    u=np.array([Fm, Fn, 0])
    z_dot=np.matmul(np.linalg.inv(Mn),(np.matmul(Bn,u)-np.matmul(Cn,z)))
    return z_dot

def Gamma(P_ped, P_rob, V_ped, V_rob, Lambda):
    dij=P_rob-P_ped
    norm_dij=np.linalg.norm(dij)
    vij=V_rob-V_ped
    norm_vij=np.linalg.norm(vij)
    cos=np.dot(vij,-dij)/(norm_dij*norm_vij)
    Gamma=Lambda+(1-Lambda)*(1+cos)/2
    return Gamma

def PTTC(P_ped, P_rob, V_ped, V_rob, PTTC_max):
    pij=P_rob-P_ped
    norm_pij=np.linalg.norm(pij)
    vij=V_rob-V_ped
    PTTC=-norm_pij*norm_pij/np.dot(pij,vij)
    if (PTTC>=PTTC_max) or (PTTC<0):
        PTTC=PTTC_max
    return PTTC
def Vel_Index(P_ped_History, P_rob_History, Speed_History, V_max, dt):
    P_rel_History=P_rob_History-P_ped_History
    n_steps=len(P_rel_History[:,0])
    i=0
    while (P_rel_History[i,0]<-5) and i<n_steps-1:
            i+=1
    t_s=i
    j=0
    while (P_rel_History[j,0]<2) and j<n_steps-1:
            j+=1
    t_e=j
    #print(t_e, t_s)
    V_vicinity=(P_rob_History[j,0]-P_rob_History[i,0])/(t_e-t_s)/dt
    Vel_Index=V_vicinity/V_max
    if t_e>n_steps-2:
        Vel_Index=0.0
    for k in range(n_steps):
        if Speed_History[k,0]<-0.01:
            Vel_Index=0.0
            break
    return Vel_Index
    

#def Trial(P_ped, Width, tau, Alpha, Beta, sigma, Lambda):
    PTTC_max=5.0
    dt=0.01
    t_end=10
    n=int(t_end/dt)+1
    t=np.linspace(0,10,n)
    m=20
    m_ped=75
    b=0.2
    I=1.216
    V_ped_0=np.array([-1.33, 0.0])
    V_max=2.5
    u_V_des=np.array([1, 0])
    V_max_vec=V_max*u_V_des
    V_des=V_max*u_V_des
    V_rob=np.array([0.0, 0.0])
    V_ped=V_ped_0
    P_rob=np.array([0.0,0.25])
    P_rob_History=P_rob
    P_ped_History=P_ped
    Speed_History=np.linalg.norm(V_rob)
    PTTC_History=np.array([PTTC_max, 0.0])
    Theta_d=0.0
    Theta=0.0
    Mn=[[m, 0], [0, m*pow(b,2)+I]]
    Bn=[[1, 0, 0], [0, b, 1]]
    v_0=0.0
    Lambda_ped=0.08
    Alpha_bnd=2.23*m_ped
    Beta_bnd=0.073
    Alpha_ped=2.65*m_ped
    Beta_ped=0.82
    tau_ped=1.62
    z_dot=np.array([float(0),  Theta_d])
    z=np.array([v_0,  float(0)])
    u_rob=np.array([math.cos(Theta), math.sin(Theta)])
    for i in range(n):
        Cn=[[0.0, -m*b*Theta_d],[m*b*Theta_d, 0.0]]
        t_i=i*dt
        f_boundary=[0.0, Myfunctions_Explicit.f_bnd(Width, P_rob, Alpha_bnd, Beta_bnd)]
        f_social_ped=Myfunctions_Explicit.f_soc_ped(P_ped, P_rob, Alpha_ped, Beta_ped, Lambda_ped, V_rob, V_ped)
        F_ped=Myfunctions_Explicit.f_des_ped(np.array([-1.33, 0]), V_ped, tau_ped, m_ped)+[0.0, Myfunctions_Explicit.f_bnd(Width, P_ped, Alpha_bnd, Beta_bnd)]+f_social_ped
        V_ped+=F_ped/m_ped*dt
        P_ped+=V_ped*dt
        f_social=Myfunctions_Explicit.f_soc(P_ped, P_rob, Alpha, Beta, V_ped, V_rob, Lambda)
        S=np.linalg.norm(f_boundary)+np.linalg.norm(f_social)
        F=Myfunctions_Explicit.f_des_rob(V_max_vec, V_rob, tau, m, S, sigma)+f_boundary+f_social
        z_dot=Myfunctions_Explicit.Robot_Dynamics(Theta, F, Mn, Cn, Bn, z)
        z+=z_dot*dt
        v=z[0]
        w=z[1]
        Theta+=w*dt
        u_rob=np.array([math.cos(Theta), math.sin(Theta)])
        V_rob=v*u_rob
        Theta_d=w   
        P_rob+=V_rob*dt
        Speed=np.linalg.norm(V_rob)
        Speed_History=np.vstack((Speed_History,Speed))
        P_rob_History=np.vstack((P_rob_History,P_rob))
        P_ped_History=np.vstack((P_ped_History,P_ped))
        PTTC_index=[Myfunctions_Explicit.PTTC(P_ped, P_rob, V_ped, V_rob, PTTC_max), t_i]
        PTTC_History=np.vstack((PTTC_History,PTTC_index))
    T_p=min(PTTC_History[:,0])
    Discomfort=33.9*np.exp(-6.5*T_p)
    ND_metric=Discomfort/5
    V_sum_vicinity=0.0
    T_vicinity=0.0
    P_rel_x=(P_rob_History[:,0]-P_ped_History[:,0])
    for i1 in range(n):
        if (P_rel_x[i1]>=-5) and (P_rel_x[i1]<=2):
            V_sum_vicinity+=Speed_History[i1]
            T_vicinity+=dt
    if T_vicinity==0:
        NV_metric=[1]
    else:
        NV_metric=1-(V_sum_vicinity/(T_vicinity/dt))/V_max
    return NV_metric, ND_metric
    #Parameters for Monte Carlo Simulations

#def MCSimu_NV_ND(Parameters):
    tau, Alpha, Beta, sigma, Lambda=Parameters
    P_ped_x0=15.0
    P_ped_y0=-0.0
    Width=np.array([2.4])
    n_width=len(Width)
    n_y0=7
    Metrics=np.zeros( (n_width*n_y0, 4))
    for i in range(n_width):
        #P_ped_y0_array=np.linspace(-Width[i]/2,Width[i]/2,n_y0+1, endpoint=False)[1:]
        P_ped_y0_array=np.array([-Width[i]/3, -Width[i]/4, -Width[i]/5, 0, Width[i]/5, Width[i]/4, Width[i]/3])
        for j in range(n_y0):
            P_ped_y0=P_ped_y0_array[j]
            P_ped=np.hstack((P_ped_x0,P_ped_y0))
            NV_metric, ND_metric=Myfunctions_Explicit.Trial(P_ped, Width[i], tau, Alpha, Beta, sigma, Lambda)
            Metrics[n_y0*i+j,0]=Width[i]
            Metrics[n_y0*i+j,1]=P_ped_y0
            Metrics[n_y0*i+j,2]=NV_metric[0]
            Metrics[n_y0*i+j,3]=ND_metric
        Cost_function=(Metrics.mean(axis=0)[2]+Metrics.mean(axis=0)[3])/2
    return Cost_function, Metrics

