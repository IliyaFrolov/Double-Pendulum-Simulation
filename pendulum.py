import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class Pendulum():

    def __init__(self, length, periods, b, F, w, initial_angular_position, initial_angular_velocity, steps):
        self.n = steps
        self.g = 9.8
        self.L = length
        self.b = b
        self.w = w
        self.F = F
        self.A = initial_angular_position
        self.P = 2*np.pi*np.sqrt(length/self.g)
        self.T = periods*self.P
        self.time = np.linspace(0, self.T, steps+1)
        self.K = KineticEnergy(length)
        self.P = PotentialEnergy(length)
        
        self.analytical = np.zeros((len(self.time)))
        self.energy_change = np.zeros(len(self.time))
        self.angular_position = np.zeros((len(self.time)))
        self.angular_velocity = np.zeros((len(self.time)))
        self.kinetic_energy = np.zeros((len(self.time)))
        self.potential_energy = np.zeros((len(self.time)))
        self.total_energy = np.zeros((len(self.time)))
        self.analytical[0] = initial_angular_position
        self.energy_change[0] = 0 
        self.angular_position[0] = initial_angular_position
        self.angular_velocity[0] = initial_angular_velocity
        self.kinetic_energy[0] = self.K(initial_angular_velocity)
        self.potential_energy[0] = self.P(initial_angular_position)
        self.total_energy[0] = self.K(initial_angular_velocity) + self.P(initial_angular_position)

        self.numerical_analysis()

    def linear_pendulum(self, parameters, t):
        omega = parameters[1]
        dwdt = -self.g/self.L*parameters[0] + self.F*np.cos(self.w*t) - self.b*self.L**2*parameters[1]

        return [omega, dwdt]
    
    def energy_analysis(self, E, t, omega):
        dedt = -self.b*self.L**2*omega**2
        
        return dedt
    
    def numerical_analysis(self):
        for i in range(1, self.n+1):
            parameters = odeint(self.linear_pendulum, [self.angular_position[i-1], self.angular_velocity[i-1]], self.time[i-1: i+1])
            energy = odeint(self.energy_analysis, self.energy_change[i-1], self.time[i-1: i+1], args=(self.angular_velocity[i-1],))
            
            self.analytical[i] = self.A*np.exp((-self.b/2)*self.time[i])*np.cos(np.sqrt(self.g/self.L - self.b**2/4)*self.time[i])
            self.energy_change[i] = energy[1, 0]
            self.angular_position[i] = parameters[1, 0]
            self.angular_velocity[i] = parameters[1, 1]
            self.kinetic_energy[i] = self.K(parameters[1, 1])
            self.potential_energy[i] = self.P(parameters[1, 0])
            self.total_energy[i] = self.K(parameters[1, 1]) + self.P(parameters[1, 0])

class KineticEnergy():

    def __init__(self, L):
        self.m = 1
        self.L = L
    
    def __call__(self, omega):
        return 0.5*self.m*(self.L*omega)**2

class PotentialEnergy():

    def __init__(self, L):
        self.m = 1
        self.g = 9.8
        self.L = L
    
    def __call__(self, angular_position):
        return self.m*self.g*self.L*(angular_position**2)/2


def plotting(pendulum, plot_energy=False):
    plt.plot(pendulum.time, pendulum.angular_position, 'r:', label='angular_position')
    plt.plot(pendulum.time, pendulum.angular_velocity, 'g--', label='Anangular_velocity velocity')
    #plt.plot(pendulum.time, pendulum.analytical, 'b-.', label='Analytical solution')
    plt.legend()
    plt.show()

    if plot_energy:
        plt.plot(pendulum.time, pendulum.kinetic_energy, 'b--', label='kinetic')
        plt.plot(pendulum.time, pendulum.potential_energy, 'g:', label='potential')
        plt.plot(pendulum.time, pendulum.total_energy, 'r-', label='total energy')
        #plt.plot(pendulum.time, pendulum.energy_change, 'y-.', label='analytical energy')
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

pendulum = Pendulum(1, 7, 0, 0, 0, np.pi*2, 0, 1000)

plotting(pendulum, plot_energy=True)
print_data(pendulum)
