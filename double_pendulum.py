import numpy as np
import pandas as pd
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class Pendulum():

    def __init__(self, length_1, length_2, mass_1, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, time, steps):
        self.n = steps
        self.g = 9.8
        self.L1 = length_1
        self.L2 = length_2
        self.m1 = mass_2
        self.m2 = mass_2
        self.time = np.linspace(0, time, steps+1)
        self.K1 = KineticEnergy1(length_1, mass_1)
        self.K2 = KineticEnergy2(length_1, length_2, mass_2)
        self.P1 = PotentialEnergy1(length_1, mass_1)
        self.P2 = PotentialEnergy2(length_1, length_2, mass_2)
    
        self.angular_position_1 = np.zeros((len(self.time)))
        self.angular_velocity_1 = np.zeros((len(self.time)))
        self.angular_position_2 = np.zeros((len(self.time)))
        self.angular_velocity_2 = np.zeros((len(self.time)))
        self.kinetic_energy = np.zeros((len(self.time)))
        self.potential_energy = np.zeros((len(self.time)))
        self.total_energy = np.zeros((len(self.time)))
        self.angular_position_1[0] = initial_angular_position_1
        self.angular_velocity_1[0] = initial_angular_velocity_1
        self.angular_position_2[0] = initial_angular_position_2
        self.angular_velocity_2[0] = initial_angular_velocity_2
        self.kinetic_energy[0] = self.K1(initial_angular_velocity_1) + self.K2(initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2)
        self.potential_energy[0] = self.P1(initial_angular_position_1) + self.P2(initial_angular_position_1, initial_angular_position_2)
        self.total_energy[0] = self.kinetic_energy[0] + self.potential_energy[0]

        self.numerical_analysis()

    def linear_pendulum(self, parameters, t):
        theta_1 = parameters[0]
        omega_1 = parameters[1]
        theta_2 = parameters[2]
        omega_2 = parameters[3]
        dwdt_1 = (-self.g*(2*self.m1+self.m2)*np.sin(theta_1)-self.m2*self.g*np.sin(theta_1-2*theta_2)-2*np.sin(theta_1-theta_2)*self.m2*(omega_2**2*self.L2+omega_1**2*self.L1*np.cos(theta_1-theta_2))) / (self.L1*(2*self.m1+self.m2-self.m2*np.cos(2*theta_1-2*theta_2)))
        dwdt_2 = (2*np.sin(theta_1-theta_2)*(omega_1**2*self.L1*(self.m1+self.m2)+self.g*(self.m1+self.m2)*np.cos(theta_1)+omega_2**2*self.L2*self.m2*np.cos(theta_1-theta_2))) / (self.L2*(2*self.m1+self.m2-self.m2*np.cos(2*theta_1-2*theta_2)))

        return [omega_1, dwdt_1, omega_2, dwdt_2]
    
    
    def numerical_analysis(self):
        for i in range(1, self.n+1):
            parameters = odeint(self.linear_pendulum, [self.angular_position_1[i-1], self.angular_velocity_1[i-1], self.angular_position_2[i-1], self.angular_velocity_2[i-1]], self.time[i-1: i+1])
            theta1 = parameters[1, 0]
            omega1 = parameters[1, 1]
            theta2 = parameters[1, 2]
            omega2 = parameters[1, 3]
        
            self.angular_position_1[i] = self.normalize_angle(theta1)
            self.angular_velocity_1[i] = omega1
            self.angular_position_2[i] = self.normalize_angle(theta2)
            self.angular_velocity_2[i] = omega2
            self.kinetic_energy[i] = self.K1(omega1) + self.K2(theta1, omega1, theta2, omega2)
            self.potential_energy[i] = self.P1(theta1) + self.P2(theta1, theta2)
            self.total_energy[i] = self.kinetic_energy[i] + self.potential_energy[i]
    
    def normalize_angle(self, angle):
        twopi = np.pi*2

        while angle > np.pi:
           angle -= twopi
        
        while angle < -np.pi:
           angle += twopi

        return angle        

class KineticEnergy1():

    def __init__(self, L1, m1):
        self.L1 = L1
        self.m1 = m1
    
    def __call__(self, omega1):
        return self.m1/2*self.L1**2*omega1**2

class KineticEnergy2():

    def __init__(self, L1, L2, m2):
        self.L1 = L1
        self.L2 = L2
        self.m2 = m2
    
    def __call__(self, theta1, omega1, theta2, omega2):
        return self.m2/2*(self.L1**2*omega1**2+self.L2**2*omega2**2+2*self.L1*self.L2*omega1*omega2*np.cos(theta1-theta2))

class PotentialEnergy1():

    def __init__(self, L1, m1):
        self.g = 9.8
        self.L1 = L1
        self.m1 = m1

    def __call__(self, theta1):
        return self.m1*self.g*self.L1*(1-np.cos(theta1))

class PotentialEnergy2():

    def __init__(self, L1, L2, m2):
        self.g = 9.8
        self.L1 = L1
        self.L2 = L2
        self.m2 = m2
    
    def __call__(self, theta1, theta2):
        return self.m2*self.g*(self.L1*(1-np.cos(theta1))+self.L2*(1-np.cos(theta2)))


def plotting(pendulum, plot_energy=False):
    plt.plot(pendulum.time, pendulum.angular_position_1, 'r-', label='angular_position_1')
    plt.plot(pendulum.time, pendulum.angular_position_2, 'g-', label='angulat_position_2')
    plt.legend()
    plt.show()

    if plot_energy:
        plt.plot(pendulum.time, pendulum.kinetic_energy, 'b-', label='kinetic')
        plt.plot(pendulum.time, pendulum.potential_energy, 'g-', label='potential')
        plt.plot(pendulum.time, pendulum.total_energy, 'r-,', label='total energy')
        plt.legend()
        plt.show() 

def print_data(pendulum):
    pendulum_data = pd.DataFrame({
    'time': [time for i, time in enumerate(pendulum.time) if i%10 == 0],
    'Angular position': [angular_position for i, angular_position in enumerate(pendulum.angular_position) if i%10 == 0],
    'Angular velocity': [angular_velocity for i, angular_velocity in enumerate(pendulum.angular_velocity) if i%10 == 0]
    })

    pendulum_data.index.name = 'Step'
    print(pendulum_data)

pendulum = Pendulum(1, 1, 1, 1, np.pi/2, 0, np.pi/4, 0, 100, 1000)
plotting(pendulum, plot_energy=True)

