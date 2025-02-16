#!/usr/bin/env python

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import ScalarFormatter

def compute_actual_rho(eta_AO, eta_OI, eta_AI, a_i, omega, k_a, k_o, alpha_a, alpha_o, H_A, H_O, z_0numA, z_0numO):
    sigma_o=np.sqrt(np.abs(omega)/(2*alpha_o))*(1+1j)
    sigma_a=np.sqrt(np.abs(omega)/(2*alpha_a))*(1+1j)
    return np.abs((1-a_i)**2*eta_AO**2/((k_o*sigma_o*(1/np.tanh(sigma_o*(H_O-z_0numO)))+(1-a_i)*eta_AO+a_i*eta_OI)*(k_a*sigma_a*(1/np.tanh(sigma_a*(H_A-z_0numA)))+(1-a_i)*eta_AO+a_i*eta_AI)))

def plot_data(rho, param, param_name, log=False, label=None, color=None, marker=None, linestyle=None):
    """
    Plot data
    """
    if log:
        rho=np.log(rho)
        x_label=r"$log(\hat{\rho})$"
    else:
        x_label=r"$\hat{\rho}$"
    
    plt.plot(param, rho, label=label, color=color, marker=marker, linestyle=linestyle)
    plt.ylabel(x_label)
    plt.xlabel(param_name)

def psi(xi, atm, stable, heat):
    if atm and stable and heat:
        return -2 / 3 * (xi - (5 / 0.35)) * np.exp(-0.35 * xi) - (1 + (2 * xi / 3))**1.5 - (10 / 1.05) + 1
    elif atm and stable and not heat:
        return -2 / 3 * (xi - (5 / 0.35)) * np.exp(-0.35 * xi) - xi - (10 / 1.05)
    elif atm and not stable and heat:
        return 2 * np.log((1 + (1 - 16 * xi)**(1 / 2)) / 2)
    elif atm and not stable and not heat:
        x = (1 - 16 * xi)**(1 / 4)
        return np.pi / 2 - 2 * np.atan(x) + np.log((1 + x)**2 * (1 + x**2) / 8)
    elif not atm and stable:
        return 1 + 5 * xi
    elif not atm and not stable:
        if heat:
            x = (1 - 25 * xi)**(1 / 3)
        else:
            x = (1 - 14 * xi)**(1 / 3)
        return np.sqrt(3) * (np.atan(np.sqrt(3)) - np.atan(1 / np.sqrt(3)) * (2 * x + 1)) + (3 / 2) * np.log((x**2 + x + 1) / 3)

def define_realistic_vals(a_i, L_AI=None, C_AO=None, u_A=None, u_O=None):
    """
    This yields an interesting plot where a^I>0 yields a larger convergence factor than a^I=0
    """
    k_o=0.58
    c_o=4182
    rho_o=1000
    alpha_o=k_o/(rho_o*c_o)
    H_O=51
    
    k_a=0.024
    c_a=1005
    rho_a=1.225
    alpha_a=k_a/(rho_a*c_a)
    H_A=210

    u_A=5 if u_A is None else u_A
    u_O=1 if u_O is None else u_O

    # Monin-Obukhov parameters
    kappa = 0.4
    z_0numA = 10 # numerical zero height atmosphere [m]
    z_0numO = 1 # numerical zero height atmosphere [m]
    z_ruAO = 2 * 10**-4
    z_ruAI = np.max([10**-3, 0.93 * 10**-3 * (1 - a_i) + 6.05 * 10**-3 * np.exp(-17 * (a_i - 0.5)**2)])
    z_rTAO = 2 * 10**-4
    z_rTAI = 10**-3
    L_AO = 50
    pre_L_AI=L_AI
    L_AI= 50 if L_AI is None else L_AI
    L_OA = 100
    nu_O = 1e-6
    nu_A = 1.5 * 1e-5
    mu = nu_O / nu_A
    lambda_u = np.sqrt(rho_a / rho_o)
    lambda_T = np.sqrt(rho_a / rho_o) * c_a / c_o
    stable_atm=True if L_AO>0 else False
    stable_oce=True if L_OA>0 else False
    C_AO = kappa**2 / ((np.log(z_0numA / z_ruAO) - psi(z_0numA / L_AO, True, stable_atm, False) + lambda_u * (np.log(lambda_u * z_0numO / (z_ruAO * mu)) - psi(z_0numO / L_OA, False, stable_oce, False))) * (np.log(z_0numA / z_rTAO) - psi(z_0numA / L_AO, True, stable_atm, True) + lambda_T * (np.log(lambda_T * z_0numO / (z_rTAO * mu)) - psi(z_0numO / L_OA, False, stable_oce, True)))) if C_AO is None else C_AO
    C_AI = kappa**2 / (np.log(z_0numA / z_ruAI) - psi(z_0numA / L_AI, True, stable_atm, False)) * (np.log(z_0numA / z_rTAI) - psi(z_0numA / L_AI, True, stable_atm, True))
    C_OI = 5*1e-3 # From Nemo
    eta_AO=C_AO*np.abs(u_A-u_O)*rho_a*c_a
    eta_AI=C_AI*np.abs(u_A)*rho_a*c_a
    eta_OI=C_OI*np.abs(u_O)*rho_o*c_o
    # C_AO=10**-5*rho_a*c_a*np.linspace(0, 30, 100)
    # C_AI=10**-5*rho_a*c_a*np.linspace(0, 30, 100) # Need help here
    # C_OI=5*10**-5*rho_o*c_o*np.linspace(0, 30, 100) # Need help here
    if not pre_L_AI is None:
        return eta_AO, eta_OI, eta_AI, k_a, k_o, alpha_a, alpha_o, H_A, H_O, z_0numA, z_0numO, C_AI
    else:
        return eta_AO, eta_OI, eta_AI, k_a, k_o, alpha_a, alpha_o, H_A, H_O, z_0numA, z_0numO

# Plot as function of a_i
# dt=2
# T=1000
# w_min=np.pi/T
# w_max=np.pi/dt
# omegas=np.linspace(0, 1, 100)*(w_max-w_min)+w_min
# a_is=np.linspace(0, 1, 100) #0.1=> nice plot

# rhos=np.zeros((len(a_is), len(omegas)))
# for i, a_i in enumerate(a_is):
#     for j, omega in enumerate(omegas):
#         eta_AO, eta_OI, eta_AI, k_a, k_o, alpha_a, alpha_o, H_A, H_O, z_0numA, z_0numO=define_realistic_vals(a_i)
#         rho=compute_actual_rho(eta_AO, eta_OI, eta_AI, a_i, omega, k_a, k_o, alpha_a, alpha_o, H_A, H_O, z_0numA, z_0numO)
#         rhos[i,j]=rho

# rho_max_over_omega=np.max(rhos, axis=1)
# print(rho_max_over_omega[0])
# plt.figure()
# plot_data(rho_max_over_omega, a_is, r"$a^I$", log=True)
# plt.show()

a_is=[0.1, 0.3, 0.5]
# Plot as function of other parameters for differetn a_i
T=1000
w_min=np.pi/T
# Ts=np.logspace(1, 6, 100)#[10,100,1000, 10000, 100000, 1000000]
# w_mins=[np.pi/T for T in Ts]

# H_As=np.linspace(10.001, 1500, 1000)#[10.001, 15, 50, 100, 200, 300, 500, 700, 1500]

# H_Os=np.linspace(1.001, 200, 1000)# [1.001, 5, 10, 50, 100, 125, 150, 175, 200]

# L_AIs = np.linspace(10, 200, 1000)#[10, 50, 100, 125, 150, 175, 200] # Some special treatment in the loop is required.
# C_AIs=np.zeros((len(a_is), len(L_AIs)))

# C_AOs = np.linspace(0.0007, 0.0013, 100)#[0.0007, 0.0008, 0.0009, 0.001, 0.0011, 0.0012, 0.0013]

# u_As = np.linspace(0,5,100)

# u_Os = np.linspace(0,5,100)

rhos=np.zeros((len(a_is), len(C_AOs)))
for i, a_i in enumerate(a_is):
    # eta_AO, eta_OI, eta_AI, k_a, k_o, alpha_a, alpha_o, H_A, H_O, z_0numA, z_0numO=define_realistic_vals(a_i)
    for j, C_AO in enumerate(C_AOs):
        # eta_AO, eta_OI, eta_AI, k_a, k_o, alpha_a, alpha_o, H_A, H_O, z_0numA, z_0numO, C_AI=define_realistic_vals(a_i, L_AI=L_AI)
        eta_AO, eta_OI, eta_AI, k_a, k_o, alpha_a, alpha_o, H_A, H_O, z_0numA, z_0numO=define_realistic_vals(a_i, L_AI=None, C_AO=C_AO, u_A=None, u_O=None)
        rho=compute_actual_rho(eta_AO, eta_OI, eta_AI, a_i, w_min, k_a, k_o, alpha_a, alpha_o, H_A, H_O, z_0numA, z_0numO)
        rhos[i,j]=rho
        # C_AIs[i,j]=C_AI


plt.figure()
colors=["b", "r", "g"]
for i, a_i in enumerate(a_is):
    a_i_rho=rhos[i,:]
    plot_data(a_i_rho, C_AOs, r"$C^{AO}$", log=False, label=r"$a^I$"+f"={a_i}", color=colors[i], linestyle="-")

# plt.xscale('log')
plt.gca().yaxis.set_major_formatter(ScalarFormatter())
plt.gca().yaxis.get_major_formatter().set_powerlimits((-6, -4))  # Adjust power limits
plt.legend(loc='upper right')
plt.show()
