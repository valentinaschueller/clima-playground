#!/usr/bin/env python

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

def compute_rho(C_AO, C_OI, C_AI, zeta, gamma, a_I):
    """
    Computing the convergence factor (with nu -> 0) and large H_O, H_A
    """
    return (1-a_I)**2*C_AO**2/(np.sqrt(((zeta+C_AO-a_I*(C_AO-C_OI))**2+zeta**2)*((gamma*zeta+C_AO-a_I*(C_AO-C_AI))**2+gamma**2*zeta**2)))

def compute_actual_rho(C_AO, C_OI, C_AI, omega, k_o, alpha_o, alpha_a, k_a, a_I, H_A, H_O):
    sigma_o=np.sqrt(np.abs(omega)/(2*alpha_o))*(1+1j)
    sigma_a=np.sqrt(np.abs(omega)/(2*alpha_a))*(1+1j)
    return (1-a_I)**2*C_AO**2/(k_o*sigma_o*(1/np.tanh(sigma_o*H_O))+(1-a_I)*C_AO+a_I*C_OI)*np.abs(k_a*sigma_a*(1/np.tanh(sigma_a*H_A))+(1-a_I)*C_AO+a_I*C_AI)

def create_video(C_AO, C_OI, C_AI, zeta, gamma, a_I, omega, k_o, alpha_o, alpha_a, k_a,  H_A, H_O, vector="C_AO"):
    global surface
    # Update function for the animation
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    A_I, Zeta = np.meshgrid(a_I, zeta)
    A_I, Omega = np.meshgrid(a_I, omega)

    # Initial surface data
    if vector=="C_AO":
        Rho=compute_actual_rho(C_AO[0], C_OI, C_AI, Omega, k_o, alpha_o, alpha_a, k_a, A_I, H_A, H_O)
        # Rho=compute_rho(C_AO[0], C_OI, C_AI, Zeta, gamma, A_I)
    elif vector=="C_OI":
        Rho=compute_actual_rho(C_AO, C_OI[0], C_AI, Omega, k_o, alpha_o, alpha_a, k_a, A_I, H_A, H_O)
        # Rho=compute_rho(C_AO, C_OI[0], C_AI, Zeta, gamma, A_I)
    elif vector=="C_AI":
        Rho=compute_actual_rho(C_AO, C_OI, C_AI[0], Omega, k_o, alpha_o, alpha_a, k_a, A_I, H_A, H_O)
        # Rho=compute_rho(C_AO, C_OI, C_AI[0], Zeta, gamma, A_I)
    print(A_I*Omega)
    print(Omega)
    # Rho=sc.fft.ifft(Rho, axis=0)#? This is a bad idea.
    # Create a 3D surface plot
    # surface = [ax.plot_surface(A_I, Zeta, Rho, cmap='viridis')]
    surface = [ax.plot_surface(A_I, Omega, Rho, cmap='viridis')]
    index=[1]

    # Axis limits
    ax.set_xlim(np.min(a_I), np.max(a_I))
    ax.set_ylim(np.min(omega), np.max(omega))
    # ax.set_ylim(np.min(zeta), np.max(zeta))
    ax.set_zlim(0, 1)
    def update(frame):
        for surf in surface:
            surf.remove()  # Remove the previous surface
        if vector=="C_AO":
            Rho=compute_actual_rho(C_AO[index], C_OI, C_AI, Omega, k_o, alpha_o, alpha_a, k_a, A_I, H_A, H_O)
            # Rho=compute_rho(C_AO[index], C_OI, C_AI, Zeta, gamma, A_I)
        elif vector=="C_OI":
            Rho=compute_actual_rho(C_AO, C_OI[index], C_AI, Omega, k_o, alpha_o, alpha_a, k_a, A_I, H_A, H_O)
            # Rho=compute_rho(C_AO, C_OI[index], C_AI, Zeta, gamma, A_I)
        elif vector=="C_AI":
            Rho=compute_actual_rho(C_AO, C_OI, C_AI[index], Omega, k_o, alpha_o, alpha_a, k_a, A_I, H_A, H_O)
            # Rho=compute_rho(C_AO, C_OI, C_AI[index], Zeta, gamma, A_I)
        # Rho=sc.fft.ifft(Rho, axis=0)# This is a bad idea.
        index[0]+=1
        surf = ax.plot_surface(A_I, Omega, Rho, cmap='viridis')  # Plot new surface
        # surf = ax.plot_surface(A_I, Zeta, Rho, cmap='viridis')  # Plot new surface
        surface[0] = surf  # Update surface reference
        ax.set_title(f"index: {index[0]-1}", fontsize=16)
        if index[0] == np.size(locals()[vector]):  # Stop the video
            ani.event_source.stop()
        return surface

    # Create animation
    ani = FuncAnimation(fig, update, frames=1000, interval=500, blit=False)

    # Display the animation
    plt.show()

def plot_data(C_AO, zeta, gamma, ice=False, C_OI=0, C_AI=0, a_I=0, to_plot=1):
    """
    Plot data
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if not ice:
        C_AO, Zeta = np.meshgrid(C_AO, zeta)
        Rho=compute_rho(C_AO, C_OI, C_AI, Zeta, gamma, a_I)
        surf = ax.plot_surface(C_AO, Zeta, Rho, cmap='viridis')
        label=f'C_AO'
    else:
        if to_plot == 1:
            C_AO, Zeta = np.meshgrid(C_AO, zeta)
            Rho=compute_rho(C_AO, C_OI, C_AI, Zeta, gamma, a_I)
            surf = ax.plot_surface(C_AO, Zeta, Rho, cmap='viridis')
        if to_plot == 2:
            C_OI, Zeta = np.meshgrid(C_OI, zeta)
            Rho=compute_rho(C_AO, C_OI, C_AI, Zeta, gamma, a_I)
            surf = ax.plot_surface(C_OI, Zeta, Rho, cmap='viridis')
        if to_plot==3:
            C_AI, Zeta = np.meshgrid(C_AI, zeta)
            Rho=compute_rho(C_AO, C_OI, C_AI, Zeta, gamma, a_I)
            surf = ax.plot_surface(C_AI, Zeta, Rho, cmap='viridis')
        label=f'lambda_{to_plot}'
        if to_plot==4:
            A_I, Zeta = np.meshgrid(a_I, zeta)
            Rho=compute_rho(C_AO, C_OI, C_AI, Zeta, gamma, A_I)
            surf = ax.plot_surface(A_I, Zeta, Rho, cmap='viridis')
        label=r"$a^I$"
    # Plot the surface
    fig.colorbar(surf)  # Add color bar for reference

    ax.set_title('Convergence Factor as a function of C_AO and zeta')
    ax.set_xlabel(label)
    ax.set_ylabel('zeta')

    plt.show()

def define_realistic_vals():
    """
    This yields an interesting plot where a^I>0 yields a larger convergence factor than a^I=0
    """
    dt=2
    T=10
    w_min=np.pi/T
    w_max=np.pi/dt
    # w_min=0.01 # 1=> nice plot.
    # w_max=10 # 10=> nice plot

    k_o=0.58
    c_o=4182
    rho_o=1000
    alpha_o=k_o/(rho_o*c_o)
    H_O=100
    
    k_a=0.024
    c_a=1005
    rho_a=1.225
    alpha_a=k_a/(rho_a*c_a)
    H_A=100

    C_AO=10**-5*rho_a*c_a*np.linspace(0, 30, 100)
    C_AI=10**-5*rho_a*c_a*np.linspace(0, 30, 100) # Need help here
    C_OI=5*10**-5*rho_o*c_o*np.linspace(0, 30, 100) # Need help here
    omega=np.linspace(0, 1, 100)*(w_max-w_min)+w_min
    zeta=np.sqrt(np.abs(omega)*k_o**2/(2*alpha_o))
    gamma=np.sqrt(alpha_o*k_a**2/(alpha_a*k_o**2))
    a_I=np.linspace(0, 1, 100)#0.1=> nice plot
    return C_AO, C_OI, C_AI, zeta, gamma, a_I, omega, k_o, alpha_o, alpha_a, k_a,  H_A, H_O

def define_vals_for_plots():
    """
    To create plots for the convergence factor cases
    """
    C_AO=1000
    C_OI=0.001
    C_AI=0.001
    zeta=0.001
    gamma=1
    a_I=0.5
    lambda_1 = (1 - a_I) * C_AO
    lambda_2 = a_I * C_OI
    lambda_3 = a_I * C_AI
    return C_AO, C_OI, C_AI, zeta, gamma, a_I

def plot_conv_facs(var, var_name):
    variables={}
    C_AO, C_OI, C_AI, zeta, gamma, a_I=define_vals_for_plots()
    var_names=("C_AO", "C_OI", "C_AI", "\zeta", "gamma", "a_I")
    variables=dict(zip(var_names, (C_AO, C_OI, C_AI, zeta, gamma, a_I)))
    variables[var_name]=var
    rho=compute_rho(*variables.values())

    if var_name[-3]=="_":
        var_name=var_name[:-3]+var_name[-3]+"{"+var_name[-2:]+"}"
    plt.figure()
    plt.plot(var, rho)
    plt.xlabel(r"${}$".format(var_name))
    plt.ylabel(r"$\rho$")
    if var_name=="a_I":
        plt.xticks(ticks=[0,1], labels=[0,1])
        plt.yticks(ticks=[0], labels=[0])
    else:
        plt.xticks(ticks=[], labels=[])
    if var_name=="C_{AO}":
        plt.plot(var, np.ones_like(var), color="grey", linestyle="--")
        plt.yticks(ticks=[0, 1], labels=[0, 1])
        plt.xticks(ticks=[0], labels=[0])
    if var_name=="\zeta":
        plt.plot(var, np.zeros_like(var), color="grey", linestyle="--")
        plt.yticks(ticks=[0], labels=[0])
        plt.xticks(ticks=[0], labels=[0])
    plt.show()

# The most interesting case is to plot rho wrt zeta and a^I, these are vectors here.
C_AO, C_OI, C_AI, zeta, gamma, a_I, omega, k_o, alpha_o, alpha_a, k_a, H_A, H_O=define_realistic_vals()
# C_AO=C_AO[int(np.size(C_AO)/2)]
C_OI=C_OI[int(np.size(C_AO)/2)]# If large, we get peak closer to zero. If small we get peak further away.
C_AI=C_AI[int(np.size(C_AI)/2)]
create_video(C_AO, C_OI, C_AI, zeta, gamma, a_I, omega,  k_o, alpha_o, alpha_a, k_a,  H_A, H_O, vector="C_AO")
# plot_conv_facs(np.linspace(0,1, 1000), "a_I")

