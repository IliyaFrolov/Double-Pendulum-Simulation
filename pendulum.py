from equations import *

class Pendulum():
    '''
    A class to represent a Pendulum Bob.
    ...
    Attributes
    ----------
    '''

    def __init__(self, length, mass, other_length, other_mass, initial_angular_position, initial_angular_velocity, other_initial_angular_position, other_initial_angular_velocity, steps, pendulum):
        self.dwdt = Acceleration(length, mass, other_length, other_mass, pendulum)
        self.K = KineticEnergy(length, mass, other_length, other_mass, pendulum)
        self.U = PotentialEnergy(length, mass, other_length, other_mass, pendulum)
        self.convert = Polar_to_Cartesian(length)

        self.L = length
        self.m = mass
        self.theta = initial_angular_position
        self.omega = initial_angular_velocity
        self.pendulum = pendulum
        self.angular_position = np.zeros(steps+1)
        self.x_position = np.zeros(steps+1)
        self.y_position = np.zeros(steps+1)
        self.angular_velocity = np.zeros(steps+1)
        self.angular_acceleration = np.zeros(steps+1)
        
        self.angular_position[0] = initial_angular_position
        self.angular_velocity[0] = initial_angular_velocity
        self.angular_acceleration[0] = self.dwdt(initial_angular_position, initial_angular_velocity, other_initial_angular_position, other_initial_angular_velocity)

    def __repr__(self):
        return f'Pendulum: {self.pendulum} ("1" for top bob "2" for bottom bob), Rod Length: {self.L}, Bob Mass: {self.m}, Initial Angular Position: {self.theta}, Initial Angular Velocity: {self.omega}'