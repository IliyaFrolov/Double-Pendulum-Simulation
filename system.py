import pandas as pd
from equations import *
from pendulum import Pendulum
from scipy.integrate import solve_ivp

class System():
    '''
    A class to represent the Double Pendulum and run the simulation.
    ...
    Attributes
    ----------
    p1 : object
        Instance of the Pendulum class representing the top pendulum.
    p2 : object
        Instance of the Pendulum class representing the bottom pendulum.
    g : int
        Gravitational acceleration constant.
    n : int
        Number of steps in the simulation.
    flip_time : int
        The time it takes for either of the pendulum bobs to flip. Is None by default if there's no flip.
    has_flipped : boolean
        Is set to True if either of the pendulum bobs has flipped, else is False.
    time : numpy array
        Stores the time of the simulation for each step.
    kinetic_energy : numpy array
        Stores the kinetic energy of the system for each step.
    potential_energy : numpy array
        Stores the potential energy of the system for each step.
    total_energy : numpy array
        Stores the total energy of the system for each step.
    
    Methods
    ----------
    model(t, initial_conditions)
        Function is used as an input into solve_ivp to get the angular position and velocity of each pendulum bob at the next step.
    run_simulation(method)
        Uses solve_ivp within a for loop to compute and store properties such as kinetic energy, potential energy, total energy, position and velocity (angular and linear) for each pendulum bob at each step. 
    check_flip(theta_1, theta_2, i)
        Checks if either of the pendulum bobs has flipped, if so, stores the time of flip.
    make_data(method)
        Calls self.run_simulation and returns all the computed data as a Pandas Dataframe.
    normalise_angle(angle)
        Normalises an angle between -pi and pi. 
    '''
    
    def __init__(self, length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps, time):
        '''
        Constructs all the necessary attributes for the System object and sets the initial values for the kinetic energy, potential energy and total energy.
        ...
        Parameters
        ----------
        length_1 : int
            Rod length of the top pendulum bob.
        mass_1 : int
            Mass of the top pendulum bob.
        length_2 : int
            Rod length of the bottom pendulum bob.
        mass_2 : int
            Mass of the bottom pendulum bob.
        initial_angular_position_1 : int
            Intitial angular displacement of the top pendulum bob from its equilibrium point.
        initial_angular_velocity_1 : int
            Initial angular velocity of the top pendulum bob (usually 0 if not at equilibrium).
        initial_angular_position_2 : int
            Intitial angular displacement of the bottom pendulum bob from its equilibrium point.
        initial_angular_velocity_2 : int
            Initial angular velocity of the bottom pendulum bob (usually 0 if not at equilibrium).
        steps: int
            Number of steps in the simulation.
        time : int
            Total time the simulation runs for.
        '''

        if length_1 == 0 or mass_1 == 0 or length_2 == 0 or mass_2 == 0:
            raise Exception('Parameters length_1, mass_2, length_2 or mass_2 cannot be equal to 0')

        self.p1 = Pendulum(1, length_1, mass_1, length_2, mass_2, self.normalise_angle(initial_angular_position_1), initial_angular_velocity_1, self.normalise_angle(initial_angular_position_2), initial_angular_velocity_2, steps)
        self.p2 = Pendulum(2, length_2, mass_2, length_1, mass_1, self.normalise_angle(initial_angular_position_2), initial_angular_velocity_2, self.normalise_angle(initial_angular_position_1), initial_angular_velocity_1, steps)

        self.g = g
        self.n = steps
        self.flip_time = None
        self.has_flipped = False
        self.time = np.linspace(0, time, steps+1)
        self.kinetic_energy = np.zeros(steps+1)
        self.potential_energy = np.zeros(steps+1)
        self.total_energy = np.zeros(steps+1)

        self.kinetic_energy[0] = self.p1.K(initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2) + self.p2.K(initial_angular_position_2, initial_angular_velocity_2, initial_angular_position_1, initial_angular_velocity_1)
        self.potential_energy[0] = self.p1.U(initial_angular_position_1, initial_angular_position_2) + self.p2.U(initial_angular_position_2, initial_angular_position_1)
        self.total_energy[0] = self.kinetic_energy[0] + self.potential_energy[0]

    def model(self, t, initial_conditions):
        '''
        Function is a parameter of solve_ivp. Return of function is used to calculate the angular displacement and velocity of each pendulum bob at the next step.
        ...
        Parameters
        ----------
        t : numpy array slice
            Time interval between current and next step.
        initial_conditions : list
            List containing the initial angular displacement and velocity of both pendulum bobs at the current step.
        
        Returns
        ----------
        list
            List containing the angular velocity and acceleration of both pendulum bobs.
        '''

        theta_1 = initial_conditions[0]
        omega_1 = initial_conditions[1]
        theta_2 = initial_conditions[2]
        omega_2 = initial_conditions[3]
        dwdt_1 = self.p1.dwdt(theta_1, omega_1, theta_2, omega_2)
        dwdt_2 = self.p2.dwdt(theta_2, omega_2, theta_1, omega_1)

        return [omega_1, dwdt_1, omega_2, dwdt_2]    
    
    def run_simulation(self, method, is_phase=False):
        '''
        Uses solve_ivp within a for loop to compute and store the position and velocity (angular and linear) for each pendulum bob at each step. 

        Parameters
        ----------
        method : string
            Used as a parameter in solve_ivp to select the approximation method to be used.
        is_phase: Boolean
            Is True if the simulation is finding the time of flip, False otherwise.
        
        Returns
        ----------
        None
        '''

        for i in range(1, self.n+1):
            output = solve_ivp(
                self.model, 
                self.time[i-1: i+1],
                [self.p1.angular_position[i-1], self.p1.angular_velocity[i-1], self.p2.angular_position[i-1], self.p2.angular_velocity[i-1]], 
                t_eval=np.linspace(self.time[i-1], self.time[i], 2), method=method
                ).y
            
            theta_1 = output[0][1]
            omega_1 = output[1][1]
            theta_2 = output[2][1]
            omega_2 = output[3][1]

            self.check_flip(theta_1, theta_2, i)

            if self.has_flipped and is_phase: # Ends the function if the flip time is found.
                return

            self.p1.angular_position[i] = self.normalise_angle(theta_1)
            self.p1.x_position[i] = self.p1.convert(theta_1, 'x')
            self.p1.y_position[i] = self.p1.convert(theta_1, 'y')
            self.p1.angular_velocity[i] = omega_1
            self.p1.angular_acceleration[i] = self.p1.dwdt(theta_1, omega_1, theta_2, omega_2)

            self.p2.angular_position[i] = self.normalise_angle(theta_2)
            self.p2.x_position[i] = self.p1.x_position[i] + self.p2.convert(theta_2, 'x')
            self.p2.y_position[i] = self.p1.y_position[i] + self.p2.convert(theta_2, 'y') 
            self.p2.angular_velocity[i] = omega_2
            self.p2.angular_acceleration[i] = self.p2.dwdt(theta_2, omega_2, theta_1, omega_1)
            
            self.kinetic_energy[i] =  self.p1.K(theta_1, omega_1, theta_2, omega_2) + self.p2.K(theta_2, omega_2, theta_1, omega_1)
            self.potential_energy[i] = self.p1.U(theta_1, theta_2) + self.p2.U(theta_2, theta_1)
            self.total_energy[i] = self.kinetic_energy[i] + self.potential_energy[i] 
        
    def check_flip(self, theta_1, theta_2, i):
        '''
        Checks if either of the pendulum bobs has flipped, if so, stores the time of flip.

        Parameters
        ----------
        theta_1 : int
            Angular displacement of the top pendulum bob.
        theta_2 : int
            Angular displacement of the bottom pendulum bob.
        i : int
            Loop variable 
        
        Returns
        ----------
        None
        '''

        if (theta_1 > pi or theta_1 < -pi or theta_2 > pi or theta_2 < -pi) and not self.has_flipped:
            self.flip_time = self.time[i]
            self.has_flipped = True

    def make_data(self, method='Radau'):
        '''
        Calls self.run_simulation and returns all the computed data as a Pandas Dataframe.

        Parameters
        ----------
        method : string
            Used as a parameter in solve_ivp to select the approximation method to be used. Is set to the 5th order Runge-Kutta approximation method by default.
        
        Returns
        ----------
        Dataframe
            A Pandas Dataframe containing the results of the simulation.
        '''

        self.run_simulation(method)

        return pd.DataFrame({
        'Time': self.time,
        'Angular position 1': self.p1.angular_position,
        'x position 1': self.p1.x_position,
        'y position 1': self.p1.y_position,
        'Angular velocity 1': self.p1.angular_velocity,
        'Angular acceleration 1': self.p1.angular_acceleration,
        'Angular position 2': self.p2.angular_position,
        'x position 2': self.p2.x_position,
        'y position 2': self.p2.y_position,
        'Angular velocity 2': self.p2.angular_velocity,
        'Angular acceleration 2': self.p2.angular_acceleration,
        'Kinetic energy': self.kinetic_energy,
        'Potential energy': self.potential_energy,
        'Total energy': self.total_energy
        })
    
    @staticmethod
    def normalise_angle(angle):
        '''
        Normalises an angle between -pi and pi.

        Parameters
        ----------
        angle: int
            Angle to be normalised.
        
        Returns
        ----------
        int
            Normalised angle.
        '''

        while angle > pi:
           angle -= 2*pi
        
        while angle < -pi:
           angle += 2*pi

        return angle  


    
