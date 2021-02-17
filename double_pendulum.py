import numpy as np
from numpy import pi, cos, sin
import pandas as pd
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class System():

    def __init__(self, length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps, time):
        self.g = 9.8
        self.p1 = Pendulum(length_1, mass_1, initial_angular_position_1, initial_angular_velocity_1, steps)
        self.p2 = Pendulum(length_2, mass_2, initial_angular_position_2, initial_angular_velocity_2, steps)
        self.p1.K = KineticEnergy1(length_1, mass_1)
        self.p1.U = PotentialEnergy1(length_1, mass_1)
        self.p2.K = KineticEnergy2(length_1, length_2, mass_2)
        self.p2.U = PotentialEnergy2(length_1, length_2, mass_2)
        self.n = steps
        self.time = np.linspace(0, time, steps+1)
        self.kinetic_energy = np.zeros(steps+1)
        self.potential_energy = np.zeros(steps+1)
        self.total_energy = np.zeros(steps+1)

        self.kinetic_energy[0] = self.p1.K(initial_angular_velocity_1) + self.p2.K(initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2)
        self.potential_energy[0] = self.p1.U(initial_angular_position_1) + self.p2.U(initial_angular_position_1, initial_angular_position_2)
        self.total_energy[0] = self.kinetic_energy[0] + self.potential_energy[0]

    def model(self, initial_conditions, t):
        theta_1 = initial_conditions[0]
        omega_1 = initial_conditions[1]
        theta_2 = initial_conditions[2]
        omega_2 = initial_conditions[3]
        dwdt_1 = (-self.g*(2*self.p1.m+self.p2.m)*sin(theta_1)-self.p2.m*self.g*sin(theta_1-2*theta_2)-2*sin(theta_1-theta_2)*self.p2.m*(omega_2**2*self.p2.L+omega_1**2*self.p1.L*cos(theta_1-theta_2))) / (self.p1.L*(2*self.p1.m+self.p2.m-self.p2.m*cos(2*theta_1-2*theta_2)))
        dwdt_2 = (2*sin(theta_1-theta_2)*(omega_1**2*self.p1.L*(self.p1.m+self.p2.m)+self.g*(self.p1.m+self.p2.m)*cos(theta_1)+omega_2**2*self.p2.L*self.p2.m*cos(theta_1-theta_2))) / (self.p2.L*(2*self.p1.m+self.p2.m-self.p2.m*cos(2*theta_1-2*theta_2)))

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
            self.kinetic_energy[i] =  self.p1.K(omega_1) + self.p2.K(theta_1, omega_1, theta_2, omega_2)
            self.potential_energy[i] = self.p1.U(theta_1) + self.p2.U(theta_1, theta_2)
            self.total_energy[i] = self.kinetic_energy[i] + self.potential_energy[i] 
    
    @staticmethod
    def normalize_angle(angle):
        while angle > pi:
           angle -= 2*pi
        
        while angle < -pi:
           angle += 2*pi

        return angle    

class Pendulum():

    def __init__(self, length, mass, initial_angular_position, initial_angular_velocity, steps):
        self.L = length
        self.m = mass
        self.angular_position = np.zeros(steps+1)
        self.angular_velocity = np.zeros(steps+1)
        self.angular_position[0] = initial_angular_position
        self.angular_velocity[0] = initial_angular_velocity   
      
class Energy():

    def __init__(self, length, mass):
        self.g = 9.8
        self.L = length
        self.m = mass
class KineticEnergy1(Energy):

    def __init__(self, length, mass):
        super().__init__(length, mass)
    
    def __call__(self, omega):
        return self.m/2*self.L**2*omega**2
class KineticEnergy2(Energy):

    def __init__(self, length_1, length_2, mass_2):
        super().__init__(length_1, mass_2)
        self.L2 = length_2
    
    def __call__(self, theta_1, omega_1, theta_2, omega_2):
        return self.m/2*(self.L**2*omega_1**2+self.L2**2*omega_2**2+2*self.L*self.L2*omega_1*omega_2*cos(theta_1-theta_2))

class PotentialEnergy1(Energy):

    def __init__(self, length, mass):
        super().__init__(length, mass)

    def __call__(self, theta):
        return self.m*self.g*self.L*(1-cos(theta))

class PotentialEnergy2(Energy):

    def __init__(self, length_1, length_2, mass_2):
        super().__init__(length_1, mass_2)
        self.L2 = length_2
         
    def __call__(self, theta_1, theta_2):
        return self.m*self.g*(self.L*(1-cos(theta_1))+self.L2*(1-cos(theta_2)))

def plotting(pendulum, plot_energy=False):
    pendulum.numerical_analysis()
    plt.plot(pendulum.time, pendulum.p1.angular_position, 'r-', label='angular_position_1')
    plt.plot(pendulum.time, pendulum.p2.angular_position, 'g-', label='angulat_position_2')
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

#pendulum = System(1, 1, 1, 1, np.pi/2, 0, np.pi/2, 0, 1000, 50)
#plotting(pendulum, True)


