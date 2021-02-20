import numpy as np
from numpy import pi, cos, sin
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class System():

    def __init__(self, length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps, time):
        self.g = 9.8
        self.p1 = Pendulum(length_1, mass_1, initial_angular_position_1, initial_angular_velocity_1, steps)
        self.p2 = Pendulum(length_2, mass_2, initial_angular_position_2, initial_angular_velocity_2, steps)
        self.p1.dwdt = Acceleration(length_1, mass_1, length_2, mass_2, 1)
        self.p2.dwdt = Acceleration(length_1, mass_1, length_2, mass_2, 2)
        self.p1.K = KineticEnergy(length_1, mass_1, length_2, mass_2, 1)
        self.p1.U = PotentialEnergy(length_1, mass_1, length_2, mass_2, 1)
        self.p2.K = KineticEnergy(length_1, mass_1, length_2, mass_2, 2)
        self.p2.U = PotentialEnergy(length_1, mass_1, length_2, mass_2, 2)
        self.n = steps
        self.time = np.linspace(0, time, steps+1)
        self.kinetic_energy = np.zeros(steps+1)
        self.potential_energy = np.zeros(steps+1)
        self.total_energy = np.zeros(steps+1)

        self.p1.angular_acceleration[0] = self.p1.dwdt(initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2)
        self.p2.angular_acceleration[0] = self.p2.dwdt(initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2)
        self.kinetic_energy[0] = self.p1.K(initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2) + self.p2.K(initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2)
        self.potential_energy[0] = self.p1.U(initial_angular_position_1, initial_angular_position_2) + self.p2.U(initial_angular_position_1, initial_angular_position_2)
        self.total_energy[0] = self.kinetic_energy[0] + self.potential_energy[0]

    def model(self, initial_conditions, t):
        theta_1 = initial_conditions[0]
        omega_1 = initial_conditions[1]
        theta_2 = initial_conditions[2]
        omega_2 = initial_conditions[3]
        dwdt_1 = self.p1.dwdt(theta_1, omega_1, theta_2, omega_2)
        dwdt_2 = self.p2.dwdt(theta_1, omega_1, theta_2, omega_2)

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
            self.p1.angular_acceleration[i] = self.p1.dwdt(theta_1, omega_1, theta_2, omega_2)
            self.p2.angular_position[i] = self.normalize_angle(theta_2)
            self.p2.angular_velocity[i] = omega_2
            self.p2.angular_acceleration[i] = self.p2.dwdt(theta_1, omega_1, theta_2, omega_2)
            self.kinetic_energy[i] =  self.p1.K(theta_1, omega_1, theta_2, omega_2) + self.p2.K(theta_1, omega_1, theta_2, omega_2)
            self.potential_energy[i] = self.p1.U(theta_1, theta_2) + self.p2.U(theta_1, theta_2)
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
        self.angular_acceleration = np.zeros(steps+1)
        self.angular_position[0] = initial_angular_position
        self.angular_velocity[0] = initial_angular_velocity
      
class EquationTerms():

    def __init__(self, length_1, mass_1, length_2, mass_2, pendulum_bob):
        self.g = 9.8
        self.L1 = length_1
        self.m1 = mass_1
        self.L2 = length_2
        self.m2 = mass_2
        self.pendulum_bob = pendulum_bob

class Acceleration(EquationTerms):

    def __init__(self, length_1, mass_1, length_2, mass_2, pendulum_bob):
        super().__init__(length_1, mass_1, length_2, mass_2, pendulum_bob)
    
    def __call__(self, theta_1, omega_1, theta_2, omega_2):
        if self.pendulum_bob == 1:
            return (-self.g*(2*self.m1+self.m2)*sin(theta_1)-self.m2*self.g*sin(theta_1-2*theta_2)-2*sin(theta_1-theta_2)*self.m2*(omega_2**2*self.L2+omega_1**2*self.L1*cos(theta_1-theta_2))) / (self.L1*(2*self.m1+self.m2-self.m2*cos(2*theta_1-2*theta_2)))

        elif self.pendulum_bob == 2:
            return (2*sin(theta_1-theta_2)*(omega_1**2*self.L1*(self.m1+self.m2)+self.g*(self.m1+self.m2)*cos(theta_1)+omega_2**2*self.L2*self.m2*cos(theta_1-theta_2))) / (self.L2*(2*self.m1+self.m2-self.m2*cos(2*theta_1-2*theta_2)))
        
        else:
            raise Exception('Parameter "pendulum_bob" must be either 1 or 2.')

class KineticEnergy(EquationTerms):

    def __init__(self, length_1, mass_1, length_2, mass_2, pendulum_bob):
        super().__init__(length_1, mass_1, length_2, mass_2, pendulum_bob)
    
    def __call__(self, theta_1, omega_1, theta_2, omega_2):
        if self.pendulum_bob == 1:
            return self.m1/2*self.L1**2*omega_1**2
        
        elif self.pendulum_bob == 2:
            return self.m1/2*(self.L1**2*omega_1**2+self.L2**2*omega_2**2+2*self.L1*self.L2*omega_1*omega_2*cos(theta_1-theta_2))
        
        else:
            raise Exception('Parameter "pendulum_bob" must be either 1 or 2.')

class PotentialEnergy(EquationTerms):

    def __init__(self, length_1, mass_1, length_2, mass_2, pendulum_bob):
        super().__init__(length_1, mass_1, length_2, mass_2, pendulum_bob)
         
    def __call__(self, theta_1, theta_2):
        if self.pendulum_bob == 1:
            return self.m1*self.g*self.L1*(1-cos(theta_1))  
        
        elif self.pendulum_bob == 2:
            return self.m1*self.g*(self.L1*(1-cos(theta_1))+self.L2*(1-cos(theta_2)))
        
        else:
            raise Exception('Parameter "pendulum_bob" must be either 1 or 2.')

def plotting(pendulum, plot_energy=False):
    pendulum.numerical_analysis()
    plt.plot(pendulum.time, pendulum.p1.angular_position, 'r-', label='angular_position_1')
    plt.plot(pendulum.time, pendulum.p2.angular_position, 'g-', label='angular_position_2')
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

#pendulum = System(1, 1, 1, 1, np.pi/2, 0, np.pi, 0, 1000, 50)
#plotting(pendulum, True)


