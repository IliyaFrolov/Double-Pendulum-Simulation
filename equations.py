import numpy as np
from numpy import sin, cos, pi
from scipy.constants import g

class Polar_to_Cartesian():
    '''
    A class representing the equation used to convert from polar to cartesian coordinates.

    Attributes
    ----------
    L : int
        Rod length of the pendulum.
    
    Methods
    ----------
    __call__(theta, coordinate)
        Used to compute and return result for the conversion from polar to cartesian coordiantes for x or y.
    '''

    def __init__(self, length):
        '''
        Construcuts the necessary attribute for the Polar_to_Cartesian object.

        Parameters
        ----------
        length : int
            Rod length of the pendulum.
        '''

        self.L = length

    def __call__(self, theta, coordinate):
        '''
        Computes and returns the result for the conversion from polar to cartesian coordiantes for x or y.

        Parameters
        ----------
        theta : int
            Angular displacement to be converted to cartesian coordinates.
        coordinate : str
            Either "x" or "y" indicating the desired coordinate to be returned from the conversion.
        
        Returns
        ----------
        int
            Returns the computed conversion of the angular displacement for the desired coordinate.
        '''

        if coordinate == 'x':
            return self.L*cos(theta-pi/2)

        elif coordinate == 'y':
            return self.L*sin(theta-pi/2)
        
        else:
            raise Exception('Parameter "coordinate" must be either "x" or "y".')

class EquationTerms():
    '''
    A parent class used to set and store the values of shared terms in the relevant double pendulum equations.

    Subclasses
    ----------
    Acceleration
    KineticEnergy
    PotentialEnergy

    Attributes
    ----------
    pendulum_bob : int
        A number either 1, 2 or 3 identifying the top pendulun bob, bottom pendulum bob or a simple pendulum respectively.
    L1 : int
        Rod length of the top pendulum bob.
    m1 : int
        Mass of the top pendulum bob.
    L2 : int
        Rod length of the bottom pendulum bob.
    m2 : int
        Mass of the bottom pendulum bob.
    '''

    def __init__(self, pendulum_bob, length, mass, other_length, other_mass):
        ''' 
        Constructs all the necessary attributes to be inherited by the subclasses.

        Parameters:
        ----------
        pendulum_bob : int
            A number either 1, 2 or 3 identifying the top pendulun bob, bottom pendulum bob or a simple pendulum respectively.
        length : int
            Rod length of the pendulum bob.
        mass : int
            Mass of the pendulum bob.
        other_length : int
            Rod length of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "length" is the rod length of the top pendulum bob and "other_length" is the rod length of the bottom pendulum bob, and vice versa.
        other_mass : int
            Mass of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "mass" is the mass of the top pendulum bob and "other_length" is the mass of the bottom pendulum bob, and vice versa.
        '''

        if pendulum_bob != 1 and pendulum_bob != 2 and pendulum_bob != 3:
            raise Exception('Parameter "pendulum_bob" must be either 1, 2 or 3.')

        self.pendulum_bob = pendulum_bob
        self.L1 = length if pendulum_bob == 1 or pendulum_bob == 3 else other_length # This assignment technique ensures that the rod lengths and pendulum bob masses are always assigned to the correct pendulum bob, i.e. self.L1 will always correspond to the rod length of the top pendulum bob.
        self.m1 = mass if pendulum_bob == 1 or pendulum_bob == 3 else other_mass
        self.L2 = other_length if pendulum_bob == 1 or pendulum_bob == 3 else length
        self.m2 = other_mass if pendulum_bob == 1 or pendulum_bob == 3 else mass

class Acceleration(EquationTerms):
    '''
    A class representing the equations of motion for the double pendulum.

    Attributes
    ----------
        All inherited from the parent class "EquationTerms".
    
    Methods
    ----------
    __call__(theta, omega, other_theta, other_omega)
        Computes and returns the result of the angular acceleration for a pendulum bob.
    '''

    def __init__(self, pendulum_bob, length, mass, other_length, other_mass):
        ''' 
        Constructs all the necessary attributes for the Acceleration object.

        Parameters:
        ----------
        pendulum_bob : int
            A number either 1, 2 or 3 identifying the top pendulun bob, bottom pendulum bob or a simple pendulum respectively.
        length : int
            Rod length of the pendulum bob.
        mass : int
            Mass of the pendulum bob.
        other_length : int
            Rod length of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "length" is the rod length of the top pendulum bob and "other_length" is the rod length of the bottom pendulum bob, and vice versa.
        other_mass : int
            Mass of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "mass" is the mass of the top pendulum bob and "other_length" is the mass of the bottom pendulum bob, and vice versa.
        '''

        super().__init__(pendulum_bob, length, mass, other_length, other_mass)
    
    def __call__(self, theta, omega, other_theta, other_omega):
        '''
        Compute and returns the result of the angular acceleration for the corresponding pendulum bob.

        Parameters
        ----------
        theta : int
            Angular displacement of the pendulum bob.
        omega : int
            Angular velocity of the pendulum bob.
        other_theta : int
            Angular displacement of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "theta" is the angular displacement of the top pendulum bob and "other_theta" is the angular displacement of the bottom pendulum bob, and vice versa.
        other_omega : int
           Angular velocity of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "omega" is the angular velocity of the top pendulum bob and "other_omega" is the angular velocity of the bottom pendulum bob, and vice versa.
        
        Returns
        ----------
        int
            Returns the computed values of the angular acceleration from the equations of motion.
        '''

        if self.pendulum_bob == 1:
            return (-g*(2*self.m1+self.m2)*sin(theta)-self.m2*g*sin(theta-2*other_theta)-2*sin(theta-other_theta)*self.m2*(other_omega**2*self.L2+omega**2*self.L1*cos(theta-other_theta))) / (self.L1*(2*self.m1+self.m2-self.m2*cos(2*theta-2*other_theta)))

        elif self.pendulum_bob == 2:
            return (2*sin(other_theta-theta)*(other_omega**2*self.L1*(self.m1+self.m2)+g*(self.m1+self.m2)*cos(other_theta)+omega**2*self.L2*self.m2*cos(other_theta-theta))) / (self.L2*(2*self.m1+self.m2-self.m2*cos(2*other_theta-2*theta)))
        
        else:
            return -g/self.L1*theta 

class KineticEnergy(EquationTerms):
    '''
    A class representing the kinetic energy for the double pendulum.

    Attributes
    ----------
        All inherited from the parent class "EquationTerms"
    
    Methods
    ----------
    __call__(theta, omega, other_theta, other_omega)
        Computes and returns the result for the kinetic energy of a pendulum bob.
    '''

    def __init__(self,pendulum_bob, length, mass, other_length, other_mass):
        ''' 
        Constructs all the necessary attributes for the KineticEnergy object.

        Parameters:
        ----------
        pendulum_bob : int
            A number either 1, 2 or 3 identifying the top pendulun bob, bottom pendulum bob or a simple pendulum respectively.
        length : int
            Rod length of the pendulum bob.
        mass : int
            Mass of the pendulum bob.
        other_length : int
            Rod length of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "length" is the rod length of the top pendulum bob and "other_length" is the rod length of the bottom pendulum bob, and vice versa.
        other_mass : int
            Mass of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "mass" is the mass of the top pendulum bob and "other_length" is the mass of the bottom pendulum bob, and vice versa.
        '''

        super().__init__(pendulum_bob, length, mass, other_length, other_mass)
    
    def __call__(self, theta, omega, other_theta, other_omega):
        '''
        Computes and returns the result for the kinetic energy of the corresponding pendulum bob.

        Parameters
        ----------
        theta : int
            Angular displacement of the pendulum bob.
        omega : int
            Angular velocity of the pendulum bob.
        other_theta : int
            Angular displacement of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "theta" is the angular displacement of the top pendulum bob and "other_theta" is the angular displacement of the bottom pendulum bob, and vice versa.
        other_omega : int
           Angular velocity of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "omega" is the angular velocity of the top pendulum bob and "other_omega" is the angular velocity of the bottom pendulum bob, and vice versa.
        
        Returns
        ----------
        int
            Returns the computed values of the kinetic energy.
        '''
        if self.pendulum_bob == 1 or self.pendulum_bob == 3: # The kinetic energy for the top pendulum bob and for a simple pendulum is the same.
            return self.m1/2*self.L1**2*omega**2
        
        else:
            return self.m2/2*(self.L1**2*other_omega**2+self.L2**2*omega**2+2*self.L1*self.L2*other_omega*omega*cos(other_theta-theta))

class PotentialEnergy(EquationTerms):
    '''
    A class representing the potential energy for the double pendulum.

    Attributes
    ----------
        All inherited from the parent class "EquationTerms"
    
    Methods
    ----------
    __call__(theta, omega, other_theta, other_omega)
        Computes and returns the result for the potential energy of a pendulum bob.
    '''
    
    def __init__(self, pendulum_bob, length, mass, other_length, other_mass):
        ''' 
        Constructs all the necessary attributes for the PotentialEnergy object.

        Parameters:
        ----------
        pendulum_bob : int
            A number either 1, 2 or 3 identifying the top pendulun bob, bottom pendulum bob or a simple pendulum respectively.
        length : int
            Rod length of the pendulum bob.
        mass : int
            Mass of the pendulum bob.
        other_length : int
            Rod length of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "length" is the rod length of the top pendulum bob and "other_length" is the rod length of the bottom pendulum bob, and vice versa.
        other_mass : int
            Mass of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "mass" is the mass of the top pendulum bob and "other_length" is the mass of the bottom pendulum bob, and vice versa.
        '''

        super().__init__(pendulum_bob, length, mass, other_length, other_mass)
         
    def __call__(self, theta, other_theta):
        '''
        Computes and returns the result of the potential energy of the corresponding pendulum bob.

        Parameters
        ----------
        theta : int
            Angular displacement of the pendulum bob.
        other_theta : int
            Angular displacement of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "theta" is the angular displacement of the top pendulum bob and "other_theta" is the angular displacement of the bottom pendulum bob, and vice versa.
        
        Returns
        ----------
        int
            Returns the computed value of the potential energy.
        '''
        if self.pendulum_bob == 1:
            return self.m1*g*self.L1*(1-cos(theta))  
        
        elif self.pendulum_bob == 2:
            return self.m2*g*(self.L1*(1-cos(other_theta))+self.L2*(1-cos(theta)))
        
        else:
            return self.m1*g*self.L1*(theta**2)/2  
        
    