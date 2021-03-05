import numpy as np
from equations import *

class Pendulum():

    def __init__(self, length, mass, other_length, other_mass, initial_angular_position, initial_angular_velocity, steps, pendulum):
        self.dwdt = Acceleration(length, mass, other_length, other_mass, pendulum)
        self.K = KineticEnergy(length, mass, other_length, other_mass, pendulum)
        self.U = PotentialEnergy(length, mass, other_length, other_mass, pendulum)
        self.convert = Polar_to_Cartesian(length)
        self.L = length
        self.m = mass
        self.angular_position = np.zeros(steps+1)
        self.x_position = np.zeros(steps+1)
        self.y_position = np.zeros(steps+1)
        self.angular_velocity = np.zeros(steps+1)
        self.angular_acceleration = np.zeros(steps+1)
        self.angular_position[0] = initial_angular_position
        self.angular_velocity[0] = initial_angular_velocity