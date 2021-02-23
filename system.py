from equations import *
from pendulum import Pendulum, np
from scipy.integrate import odeint

class System():

    def __init__(self, length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps, time):
        self.g = 9.8
        self.p1 = Pendulum(length_1, mass_1, initial_angular_position_1, initial_angular_velocity_1, steps)
        self.p2 = Pendulum(length_2, mass_2, initial_angular_position_2, initial_angular_velocity_2, steps)
        self.p1.convert = Polar_to_Cartesian(length_1)
        self.p2.convert = Polar_to_Cartesian(length_2)
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
            self.p1.x_position[i] = self.p1.convert(theta_1, 'x')
            self.p1.y_position[i] = self.p1.convert(theta_1, 'y')
            self.p1.angular_velocity[i] = omega_1
            self.p1.angular_acceleration[i] = self.p1.dwdt(theta_1, omega_1, theta_2, omega_2)
            self.p2.angular_position[i] = self.normalize_angle(theta_2)
            self.p2.x_position[i] = self.p1.convert(theta_2, 'x') + self.p1.x_position[i]
            self.p2.y_position[i] = self.p1.convert(theta_2, 'y') + self.p1.y_position[i]
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