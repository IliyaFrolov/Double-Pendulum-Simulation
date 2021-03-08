import numpy as np
from numpy import sin, cos, pi

class Polar_to_Cartesian():
    '''
    A class representing the equation used to convert from polar to cartesian coordinates.

    Attributes
    ----------
    L : int
        Rod length of the pendulum.
    
    Methods
    ---------
    __call__(theta, coordinate)
        Dunder method used to compute and return result for the conversion from polar to cartesian coordiantes for x or y.
    '''

    def __init__(self, length):
        '''
        Construcuts the necessary attribute for the Polar_to_Cartesian object.

        Parameters
        ---------
        length : int
            Rod length of the pendulum
        '''

        self.L = length

    def __call__(self, theta, coordinate):
        '''
        Dunder method used to compute and return result for the conversion from polar to cartesian coordiantes for x or y.

        Parameters
        ---------
        theta : int, required
            Angular displacement to be converted to cartesian coordinates.
        coordinate : str, required
            Either "x" or "y" indicating the desired coordinate to be returned from the conversion.
        
        Returns
        ---------
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
    --------
    Acceleration
    KineticEnergy
    PotentialEnergy

    Attributes
    ---------
    pendulum_bob : int
        A number either 1 or 2, identifying the top and bottom pendulum respectively.
    g : int
        Gravitational acceleration constant.
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
        ---------
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
        '''

        if pendulum_bob != 1 and pendulum_bob != 2:
            raise Exception('Parameter "pendulum_bob" must be either 1 or 2.')

        self.pendulum_bob = pendulum_bob
        self.g = 9.8
        self.L1 = length if pendulum_bob == 1 else other_length
        self.m1 = mass if pendulum_bob == 1 else other_mass
        self.L2 = other_length if pendulum_bob == 1 else length
        self.m2 = other_mass if pendulum_bob == 1 else mass

class Acceleration(EquationTerms):
    '''
    A class representing the equations of motion for the double pendulum.

    Attributes
    ----------
        All inherited from the parent class "EquationTerms"
    
    Methods
    ----------
    __call__(theta, omega, other_theta, other_omega)
        Dunder method used to compute and return the result of the angular acceleration for a pendulum bob.
    '''

    def __init__(self, pendulum_bob, length, mass, other_length, other_mass):
        ''' 
        Constructs all the necessary attributes for the Acceleration object.

        Parameters:
        ---------
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
        '''

        super().__init__(pendulum_bob, length, mass, other_length, other_mass)
    
    def __call__(self, theta, omega, other_theta, other_omega):
        '''
        Dunder method used to compute and return the result of the angular acceleration for a pendulum bob.

        Parameters
        ---------
        theta : int, required
            Angular displacement of the pendulum bob.
        omega : int, required
            Angular velocity of the pendulum bob.
        other_theta : int, required
            Angular displacement of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "theta" is the angular displacement of the top pendulum bob and "other_theta" is the angular displacement of the bottom pendulum bob, and vice versa.
        other_omega : int, required
           Angular velocity of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "omega" is the angular velocity of the top pendulum bob and "other_omega" is the angular velocity of the bottom pendulum bob, and vice versa.
        
        Returns
        ---------
        int
            Returns the computed values of the angular acceleration from the equations of motion.
        '''
        if self.pendulum_bob == 1:
            return (-self.g*(2*self.m1+self.m2)*sin(theta)-self.m2*self.g*sin(theta-2*other_theta)-2*sin(theta-other_theta)*self.m2*(other_omega**2*self.L2+omega**2*self.L1*cos(theta-other_theta))) / (self.L1*(2*self.m1+self.m2-self.m2*cos(2*theta-2*other_theta)))

        else:
            return (2*sin(other_theta-theta)*(other_omega**2*self.L1*(self.m1+self.m2)+self.g*(self.m1+self.m2)*cos(other_theta)+omega**2*self.L2*self.m2*cos(other_theta-theta))) / (self.L2*(2*self.m1+self.m2-self.m2*cos(2*other_theta-2*theta)))
        

class KineticEnergy(EquationTerms):
    '''
    A class representing the kinetic energy for the double pendulum.

    Attributes
    ----------
        All inherited from the parent class "EquationTerms"
    
    Methods
    ----------
    __call__(theta, omega, other_theta, other_omega)
        Dunder method used to compute and return the result for the kinetic energy of a pendulum bob.
    '''

    def __init__(self,pendulum_bob, length, mass, other_length, other_mass):
        ''' 
        Constructs all the necessary attributes for the KineticEnergy object.

        Parameters:
        ---------
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
        '''

        super().__init__(pendulum_bob, length, mass, other_length, other_mass)
    
    def __call__(self, theta, omega, other_theta, other_omega):
        '''
        Dunder method used to compute and return the result for the kinetic energy of a pendulum bob.

        Parameters
        ---------
        theta : int, required
            Angular displacement of the pendulum bob.
        omega : int, required
            Angular velocity of the pendulum bob.
        other_theta : int, required
            Angular displacement of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "theta" is the angular displacement of the top pendulum bob and "other_theta" is the angular displacement of the bottom pendulum bob, and vice versa.
        other_omega : int, required
           Angular velocity of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "omega" is the angular velocity of the top pendulum bob and "other_omega" is the angular velocity of the bottom pendulum bob, and vice versa.
        
        Returns
        ---------
        int
            Returns the computed values of the kinetic energy from the equations for the kinetic energy of the double pendulum.
        '''
        if self.pendulum_bob == 1:
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
        Dunder method used to compute and return the result for the potential energy of a pendulum bob.
    '''
    
    def __init__(self, pendulum_bob, length, mass, other_length, other_mass):
        ''' 
        Constructs all the necessary attributes for the PotentialEnergy object.

        Parameters:
        ---------
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
        '''

        super().__init__(pendulum_bob, length, mass, other_length, other_mass)
         
    def __call__(self, theta, other_theta):
        '''
        Dunder method used to compute and return the result of the potential energy for a pendulum bob.

        Parameters
        ---------
        theta : int, required
            Angular displacement of the pendulum bob.
        other_theta : int, required
            Angular displacement of the other pendulum bob, i.e. if parameter "pendulum_bob" is 1, "theta" is the angular displacement of the top pendulum bob and "other_theta" is the angular displacement of the bottom pendulum bob, and vice versa.
        
        Returns
        ---------
        int
            Returns the computed values of the potential energy from the equations for the potential energy of the double pendulum.
        '''
        if self.pendulum_bob == 1:
            return self.m1*self.g*self.L1*(1-cos(theta))  
        
        else:
            return self.m2*self.g*(self.L1*(1-cos(other_theta))+self.L2*(1-cos(theta)))
        
    