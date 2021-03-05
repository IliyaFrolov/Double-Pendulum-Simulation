from numpy import sin, cos, pi

class Polar_to_Cartesian():

    def __init__(self, length):
        self.L = length

    def __call__(self, theta, coordinate):
        if coordinate == 'x':
            return self.L*cos(theta-pi/2)

        elif coordinate == 'y':
            return self.L*sin(theta-pi/2)
        
        else:
            raise Exception('Parameter "coordinate" must be either "x" or "y".')

class EquationTerms():

    def __init__(self, length, mass, other_length, other_mass, pendulum_bob):
        if pendulum_bob != 1 and pendulum_bob != 2:
            raise Exception('Parameter "pendulum_bob" must be either 1 or 2.')

        self.g = 9.8
        self.L1 = length if pendulum_bob == 1 else other_length
        self.m1 = mass if pendulum_bob == 1 else other_mass
        self.L2 = length if pendulum_bob == 2 else other_length
        self.m2 = mass if pendulum_bob == 2 else other_mass
        self.pendulum_bob = pendulum_bob

class Acceleration(EquationTerms):

    def __init__(self, length, mass, other_length, other_mass, pendulum_bob):
        super().__init__(length, mass, other_length, other_mass, pendulum_bob)
    
    def __call__(self, theta_1, omega_1, theta_2, omega_2):
        if self.pendulum_bob == 1:
            return (-self.g*(2*self.m1+self.m2)*sin(theta_1)-self.m2*self.g*sin(theta_1-2*theta_2)-2*sin(theta_1-theta_2)*self.m2*(omega_2**2*self.L2+omega_1**2*self.L1*cos(theta_1-theta_2))) / (self.L1*(2*self.m1+self.m2-self.m2*cos(2*theta_1-2*theta_2)))

        else:
            return (2*sin(theta_1-theta_2)*(omega_1**2*self.L1*(self.m1+self.m2)+self.g*(self.m1+self.m2)*cos(theta_1)+omega_2**2*self.L2*self.m2*cos(theta_1-theta_2))) / (self.L2*(2*self.m1+self.m2-self.m2*cos(2*theta_1-2*theta_2)))
        

class KineticEnergy(EquationTerms):

    def __init__(self, length, mass, other_length, other_mass, pendulum_bob):
        super().__init__(length, mass, other_length, other_mass, pendulum_bob)
    
    def __call__(self, theta_1, omega_1, theta_2, omega_2):
        if self.pendulum_bob == 1:
            return self.m1/2*self.L1**2*omega_1**2
        
        else:
            return self.m2/2*(self.L1**2*omega_1**2+self.L2**2*omega_2**2+2*self.L1*self.L2*omega_1*omega_2*cos(theta_1-theta_2))

class PotentialEnergy(EquationTerms):

    def __init__(self, length, mass, other_length, other_mass, pendulum_bob):
        super().__init__(length, mass, other_length, other_mass, pendulum_bob)
         
    def __call__(self, theta_1, theta_2):
        if self.pendulum_bob == 1:
            return self.m1*self.g*self.L1*(1-cos(theta_1))  
        
        else:
            return self.m2*self.g*(self.L1*(1-cos(theta_1))+self.L2*(1-cos(theta_2)))
        
    