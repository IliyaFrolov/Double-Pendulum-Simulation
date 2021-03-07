from equations import *

class Pendulum():
    '''
    A class to represent a Pendulum Bob.
    ...
    Attributes
    ----------
    dwdt : object
        Instance of the Acceleration class representing the equations of motion.
    K : object
        Instance of the KineticEnergy class representing the kinetic energy equation for a pendulum bob.
    U : object
        Instance of the PotentialEnergy class representing the potential energy equation for a pendulum bob.
    convert : object
        Instance of the Polar_to_Cartesian class representing the equation to convert from polar to cartesian coordinates.
    pendulum_bob : int
        A number either 1 or 2, identifying the top and bottom pendulum respectively.
    L : int
        Rod length of the pendulum bob.
    m : int
        Mass of the pendulum bob.
    theta : int
        Initial angular displacement of the pendulum bob.
    omega : int
        Initial angular velocity of the pendulum bob.
    angular_position : numpy array
        Stores the angular displacement of the pendulum bob at each step.
    x_position : numpy array
        Stores the linear displacement of the pendulum bob in the x direction at each step.
    y_position : numpy array
        Stores the linear displacement of the pendulum bob in the y direction at each step.
    angular_velocity : numpy array
        Stores the angular velocity of the pendulum bob at each step.
    angular_acceleration : numpy array
        Stores the angular acceleration of the pendulum bob at each step.
    
    Methods
    ---------
    __repr__()
        A dunder method used to print a summary of information about the pendulum bob.
    '''

    def __init__(self, pendulum_bob, length, mass, other_length, other_mass, initial_angular_position, initial_angular_velocity, other_initial_angular_position, other_initial_angular_velocity, steps):
        '''
        Constructs all the necessary attributes for the Pendulum object and sets the initial values for the angular displacement, velocity and acceleration of the pendulum bob.
        
        Parameters
        ----------
        pendulum_bob : int
            A number either 1 or 2, identifying the top and bottom pendulum respectively.
        length : int
            Rod length of the pendulum bob.
        mass : int
            Mass of the pendulum bob.
        other_length : int
            Rod length of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "length" is the rod length of the top pendulum bob and "other_length" is the rod length of the bottom pendulum bob, and vice versa.
        other_mass : int
            Mass of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "mass" is the mass of the top pendulum bob and "other_length" is the mass of the bottom pendulum bob, and vice versa.
        initial_angular_position : int
            Intitial angular displacement of the pendulum bob from its equilibrium point.
        initial_angular_velocity : int
            Initial angular velocity of the pendulum bob (usually 0 if not at equilibrium).
        other_initial_angular_position : int
            Intitial angular displacement of the other pendulum bob from its equilibrium point, i.e. if parameter "pendulum_bob" is 1, "initial_angular_position" is the initial angular displacement of the top pendulum bob and "other_intial_angular_position" is the initial_angular_displacement of the bottom pendulum bob, and vice versa.
        other_initial_angular_velocity_2 : int
             Intitial angular velocity of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "initial_angular_velocity" is the initial angular displacement of the top pendulum bob and "other_intial_angular_velocity" is the initial_angular_velocity of the bottom pendulum bob, and vice versa.
        steps: int
            Number of steps in the simulation.
        '''

        self.dwdt = Acceleration(pendulum_bob, length, mass, other_length, other_mass)
        self.K = KineticEnergy(pendulum_bob, length, mass, other_length, other_mass)
        self.U = PotentialEnergy(pendulum_bob, length, mass, other_length, other_mass)
        self.convert = Polar_to_Cartesian(length)

        self.pendulum_bob = pendulum_bob
        self.L = length
        self.m = mass
        self.theta = initial_angular_position
        self.omega = initial_angular_velocity
        self.angular_position = np.zeros(steps+1)
        self.x_position = np.zeros(steps+1)
        self.y_position = np.zeros(steps+1)
        self.angular_velocity = np.zeros(steps+1)
        self.angular_acceleration = np.zeros(steps+1)
        
        self.angular_position[0] = initial_angular_position
        self.angular_velocity[0] = initial_angular_velocity
        self.angular_acceleration[0] = self.dwdt(initial_angular_position, initial_angular_velocity, other_initial_angular_position, other_initial_angular_velocity)

    def __repr__(self):
        '''
        A dunder method used to print a summary of information about the pendulum bob.

        Parameters
        ---------
        None

        Returns
        ---------
        str
            A string used to present the summary of information about the pendulum rod including the pendulum bob number, rod length, mass, initial angular displacement and velocity.
        '''
        
        return f'Pendulum: {self.pendulum_bob} ("1" for top bob "2" for bottom bob), Rod Length: {self.L}, Bob Mass: {self.m}, Initial Angular Position: {self.theta}, Initial Angular Velocity: {self.omega}'