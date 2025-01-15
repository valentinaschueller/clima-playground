#!/usr/bin/env python

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

def create_video(data_atm, t, z):
    global surface
    # Update function for the animation
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    T, Z = np.meshgrid(t, z)

    surface = [ax.plot_surface(T, Z, data_atm[:,:,0], cmap='viridis')]
    index=[0]

    # Axis limits
    ax.set_xlim(np.min(t), np.max(t))
    ax.set_ylim(np.min(z), np.max(z))
    # ax.set_ylim(np.min(zeta), np.max(zeta))
    ax.set_zlim(np.min(data_atm), np.max(data_atm))
    def update(frame):
        if index[0] == np.shape(data_atm)[2]:  # Stop the video
            print(f"Stopped at iter {index[0]}")
            ani.event_source.stop()
            return surface
        for surf in surface:
            surf.remove()  # Remove the previous surface
        surf = ax.plot_surface(T, Z, data_atm[:,:,index[0]], cmap='viridis')  # Plot new surface
        # surf = ax.plot_surface(A_I, Zeta, Rho, cmap='viridis')  # Plot new surface
        surface[0] = surf  # Update surface reference
        ax.set_title(f"iter: {index[0]+1}", fontsize=16)
        index[0]+=1
        return surface

    # Create animation
    ani = FuncAnimation(fig, update, frames=1, interval=1000, blit=False)

    # Display the animation
    plt.show()

def main():
    # Get the list of files in the directories
    atm_files = os.listdir("atm_temps/")
    oce_files = os.listdir("oce_temps/")

    # Initialize lists to store the data frames
    atm_data_frames = []
    oce_data_frames = []
    times=[]
    atm_z_range=[]
    oce_z_range=[]

    # Read the ATM temperature files and store the data frames
    for i, atm_file in enumerate(atm_files):
        file_path = os.path.join("atm_temps", atm_file)
        df = pd.read_csv(file_path)  # Assuming no header in the CSV file
        if i==0:
            times=list(map(float, df.columns.tolist()[1:]))
            atm_z_range=df.index.tolist()
        temps=df.values[:,1:]
        atm_data_frames.append(temps)

    # Read the OCE temperature files and store the data frames
    for oce_file in oce_files:
        file_path = os.path.join("oce_temps", oce_file)
        df = pd.read_csv(file_path)  # Assuming no header in the CSV file
        if i==0:
            oce_z_range=df.index.tolist()
        temps=df.values[:,1:]
        oce_data_frames.append(temps)

    # You can check the data if needed:
    print(atm_data_frames[0])  # Print the first ATM data frame
    print(oce_data_frames[0])  # Print the first OCE data frame


    # Convert the data frames to numpy arrays and concatenate along the third axis
    data_atm = np.stack([iter_data for iter_data in atm_data_frames], axis=2)
    data_oce = np.stack([iter_data for iter_data in oce_data_frames], axis=2)

    # Print the shape of the resulting 3D arrays
    print(data_atm.shape)
    print(data_oce.shape)
    create_video(data_atm,times,atm_z_range)

if __name__ == '__main__':
    main()
