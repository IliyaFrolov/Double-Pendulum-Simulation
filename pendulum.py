from equations import *
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os

class Pendulum():
    '''
    A class to represent a Pendulum Bob and run a simulation of a single pendulum.
    ...
    Attributes
    ----------
    dwdt : object
        Instance of the Acceleration class representing the equations of motion for the double pendulum.
    K : object
        Instance of the KineticEnergy class representing the kinetic energy equations for a double pendulum.
    U : object
        Instance of the PotentialEnergy class representing the potential energy equations for a double pendulum.
    convert : object
        Instance of the Polar_to_Cartesian class representing the equation to convert from polar to cartesian coordinates.
    pendulum_bob : int
        A number either 1 or 2, identifying the top and bottom pendulum respectively.
    L : int
        Rod length of the pendulum bob.
    m : int
        Mass of the pendulum bob.
    n : int
        Number of steps in the simple pendulum simulation.
    time : numpy array
        Stores the time of the simple pendulum simulation for each step.
    angular_position : numpy array
        Stores the angular displacement of the pendulum bob for each step.
    x_position : numpy array
        Stores the linear displacement of the pendulum bob in the x direction for each step.
    y_position : numpy array
        Stores the linear displacement of the pendulum bob in the y direction for each step.
    angular_velocity : numpy array
        Stores the angular velocity of the pendulum bob for each step.
    angular_acceleration : numpy array
        Stores the angular acceleration of the pendulum bob for each step.
    
    Methods
    ----------
    init_simple_pendulum(length, mass, initial_angular_position, initial_anglar_velocity, steps, time)
        Used as an alternative to instantiate a Pendulum object for running the simple pendulum simulation.  
    __repr__()
        Shows a summary of information about the pendulum bob when printing a Pendulum object.
    model(t, initial_conditions)
        Function is used as an input into solve_ivp to get the angular position and velocity of each pendulum bob at the next step.
    make_simple_pendulum
        Uses solve_ivp to obtain and plot the solutions for the simple pendulum at each step.
    '''

    def __init__(self, pendulum_bob, length, mass, other_length, other_mass, initial_angular_position, initial_angular_velocity, other_initial_angular_position, other_initial_angular_velocity, steps, time=0):
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
            Number of steps in the simple pendulum simulation.
        time : int
            Total time the simple pendulum simulation runs for. Is 0 by default if running double pendulum simulation instead.
        '''

        self.dwdt = Acceleration(pendulum_bob, length, mass, other_length, other_mass)
        self.K = KineticEnergy(pendulum_bob, length, mass, other_length, other_mass)
        self.U = PotentialEnergy(pendulum_bob, length, mass, other_length, other_mass)
        self.convert = Polar_to_Cartesian(length)

        self.pendulum_bob = pendulum_bob
        self.L = length
        self.m = mass
        self.n = steps
        self.time = np.linspace(0, time, steps+1)
        self.x_position = np.zeros(steps+1)
        self.y_position = np.zeros(steps+1)
        self.angular_position = np.zeros(steps+1)
        self.angular_velocity = np.zeros(steps+1)
        self.angular_acceleration = np.zeros(steps+1)
       
        self.angular_position[0] = initial_angular_position
        self.angular_velocity[0] = initial_angular_velocity
        self.angular_acceleration[0] = self.dwdt(initial_angular_position, initial_angular_velocity, other_initial_angular_position, other_initial_angular_velocity)  

    @classmethod
    def init_simple_pendulum(cls, length, mass, initial_angular_position, initial_angular_velocity, steps, time):
        '''
        Used as an alternative to instantiate a Pendulum object for running the simple pendulum simulation.

        Parameters
        ----------
        length : int
            Rod length of the pendulum bob.
        mass : int
            Mass of the pendulum bob.
        initial_angular_position : int
            Intitial angular displacement of the pendulum bob from its equilibrium point.
        initial_angular_velocity : int
            Initial angular velocity of the pendulum bob (usually 0 if not at equilibrium).
        steps: int
            Number of steps in the simple pendulum simulation.
        time : int
            Total time the simple pendulum simulation runs for. Is 0 by default if running double pendulum simulation instead.
        
        Returns
        ----------
        object
            Returns a Pendulum object to run a simple pendulum simulation.
        '''
        if length == 0 or mass == 0:
            raise Exception('Parameters length or mass cannot be equal to 0')

        return cls(3, length, mass, 0, 0, initial_angular_position, initial_angular_velocity, 0, 0, steps, time)

    def __repr__(self):
        '''
        Shows a summary of information about the pendulum bob when printing a Pendulum object.

        Parameters
        ----------
        None

        Returns
        ----------
        str
            A string used to present the summary of information about the pendulum rod including the pendulum bob number, rod length, mass, initial angular displacement and velocity.
        '''
        
        return f'Pendulum: {self.pendulum_bob} ("1" for top bob "2" for bottom bob), Rod Length: {self.L}, Bob Mass: {self.m}, Initial Angular Position: {self.angular_position[0]}, Initial Angular Velocity: {self.angular_velocity[0]}, Initial Angular Acceleration: {self.angular_acceleration[0]}'

    def model(self, t, initial_conditions):
        '''
        Function is a parameter of solve_ivp. Return of the function is used to calculate the angular displacement and velocity of the simple pendulum.
        ...
        Parameters
        ----------
        t : numpy array
            Time interval between current and next step.
        initial_conditions : list
            List containing the initial angular displacement and velocity of the simple pendulum at the current step.

        Returns
        ----------
        list
            List containing the angular velocity and acceleration of the simple pendulum.
        '''

        theta = initial_conditions[0]
        omega = initial_conditions[1]
        dwdt = -g/self.L*theta 

        return [omega, dwdt]  
    
    def make_simple_pendulum(self, method='Radau', file_name=None):
        '''
        Uses solve_ivp within a for loop to obtain and plot the solutions for the simple pendulum at each step.

        Parameters
        ----------
        method : string
            Used as a parameter in solve_ivp to select the approximation method to be used. Is set to the 5th order Runge-Kutta approximation method by default.
        file_name : str
            File name of the produced plot to save. If no name is specified, plot is not saved.
        
        Returns
        ----------
        None
        '''

        for i in range(1, self.n+1):
            output = solve_ivp(self.model, self.time[i-1: i+1],[self.angular_position[i-1], self.angular_velocity[i-1]], t_eval=np.linspace(self.time[i-1], self.time[i], 2), method=method).y
            theta = output[0][1]
            omega = output[1][1]

            self.angular_position[i] = theta
            self.angular_velocity[i] = omega
            self.angular_acceleration[i] = -g/self.L*theta 

        kinetic_energy = self.K(self.angular_position, self.angular_velocity, 0, 0)
        potential_energy = self.U(self.angular_position, 0)
        total_energy = kinetic_energy + potential_energy
        
        fig_1 = plt.figure()
        ax = fig_1.add_subplot(111)
        ax.minorticks_on()
        ax.grid(which='major', linestyle='-', linewidth='0.5', color='black')
        ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
        ax.set(title='Angular displacements of the Simple Pendulum over time.', ylabel='Angular displacement (radians)', xlabel='Time (s)')
        ax.plot(self.time, self.angular_position, 'r-', label='Angular position')
        ax.plot(self.time, self.angular_velocity, 'b-', label='Angular velocity')
        ax.plot(self.time, self.angular_acceleration, 'g-', label='Angular acceleration')
        ax.legend(loc='upper left')
        plt.show()

        fig_2 = plt.figure()
        ax = fig_2.add_subplot(111)
        ax.set(title='Energy of the Simple Pendulum over time.', ylabel='Energy (J)', xlabel='Time (s)')
        ax.plot(self.time, kinetic_energy, 'b-', label='Kinetic Energy')
        ax.plot(self.time, potential_energy, 'g-', label='Potential Energy')
        ax.plot(self.time, total_energy, 'r-', label='Total Energy')
        ax.legend(loc='upper left')
        plt.show()

        if file_name:
            fig_1.savefig(rf'{os.getcwd()}\{file_name}.png')
            fig_2.savefig(rf'{os.getcwd()}\{file_name}.png')

