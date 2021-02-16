import numpy as np
import pandas as pd
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class Pendulum():

    def __init__(self, length, mass, initial_angular_position, initial_angular_velocity, steps):
        self.L = length
        self.m = mass
        self.angular_position = np.zeros(steps+1)
        self.angular_velocity = np.zeros(steps+1)
        self.angular_position[0] = initial_angular_position
        self.angular_velocity[0] = initial_angular_velocity

class System():

    def __init__(self, length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps, time):
        self.g = 9.8
        self.p1 = Pendulum(length_1, mass_1, initial_angular_position_1, initial_angular_velocity_1, steps)
        self.p2 = Pendulum(length_2, mass_2, initial_angular_position_2, initial_angular_velocity_2, steps)
        self.n = steps
        self.time = np.linspace(0, time, steps+1)
    
    def model(self, initial_conditions, t):
        theta_1 = initial_conditions[0]
        omega_1 = initial_conditions[1]
        theta_2 = initial_conditions[2]
        omega_2 = initial_conditions[3]
        dwdt_1 = (-self.g*(2*self.p1.m+self.p2.m)*np.sin(theta_1)-self.p2.m*self.g*np.sin(theta_1-2*theta_2)-2*np.sin(theta_1-theta_2)*self.p2.m*(omega_2**2*self.p2.L+omega_1**2*self.p1.L*np.cos(theta_1-theta_2))) / (self.p1.L*(2*self.p1.m+self.p2.m-self.p2.m*np.cos(2*theta_1-2*theta_2)))
        dwdt_2 = (2*np.sin(theta_1-theta_2)*(omega_1**2*self.p1.L*(self.p1.m+self.p2.m)+self.g*(self.p1.m+self.p2.m)*np.cos(theta_1)+omega_2**2*self.p2.L*self.p2.m*np.cos(theta_1-theta_2))) / (self.p2.L*(2*self.p1.m+self.p2.m-self.p2.m*np.cos(2*theta_1-2*theta_2)))

        return [omega_1, dwdt_1, omega_2, dwdt_2]    
    
    def numerical_analysis(self):
        for i in range(1, self.n+1):
            output = odeint(
                self.model, 
                [self.p1.angular_position[i-1], self.p1.angular_velocity[i-1], self.p2.angular_position[i-1], self.p2.angular_velocity[i-1]], 
                self.time[i-1: i+1]
                )
            theta_1 = output[1, 0]
            omega_1 = output[1, 1]
            theta_2 = output[1, 2]
            omega_2 = output[1, 3]
        
            self.p1.angular_position[i] = self.normalize_angle(theta_1)
            self.p1.angular_velocity[i] = omega_1
            self.p2.angular_position[i] = self.normalize_angle(theta_2)
            self.p2.angular_velocity[i] = omega_2
    
    @staticmethod
    def normalize_angle(angle):
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


