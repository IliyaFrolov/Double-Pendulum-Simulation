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
        self.L2 = other_length if pendulum_bob == 1 else length
        self.m2 = other_mass if pendulum_bob == 1 else mass
        self.pendulum_bob = pendulum_bob

class Acceleration(EquationTerms):

    def __init__(self, length, mass, other_length, other_mass, pendulum_bob):
        super().__init__(length, mass, other_length, other_mass, pendulum_bob)
    
    def __call__(self, theta, omega, other_theta, other_omega):
        if self.pendulum_bob == 1:
            return (-self.g*(2*self.m1+self.m2)*sin(theta)-self.m2*self.g*sin(theta-2*other_theta)-2*sin(theta-other_theta)*self.m2*(other_omega**2*self.L2+omega**2*self.L1*cos(theta-other_theta))) / (self.L1*(2*self.m1+self.m2-self.m2*cos(2*theta-2*other_theta)))

        else:
            return (2*sin(other_theta-theta)*(other_omega**2*self.L1*(self.m1+self.m2)+self.g*(self.m1+self.m2)*cos(other_theta)+omega**2*self.L2*self.m2*cos(other_theta-theta))) / (self.L2*(2*self.m1+self.m2-self.m2*cos(2*other_theta-2*theta)))
        

class KineticEnergy(EquationTerms):

    def __init__(self, length, mass, other_length, other_mass, pendulum_bob):
        super().__init__(length, mass, other_length, other_mass, pendulum_bob)
    
    def __call__(self, theta, omega, other_theta, other_omega):
        if self.pendulum_bob == 1:
            return self.m1/2*self.L1**2*omega**2
        
        else:
            return self.m2/2*(self.L1**2*omega**2+self.L2**2*other_omega**2+2*self.L1*self.L2*omega*other_omega*cos(theta-other_theta))

class PotentialEnergy(EquationTerms):

    def __init__(self, length, mass, other_length, other_mass, pendulum_bob):
        super().__init__(length, mass, other_length, other_mass, pendulum_bob)
         
    def __call__(self, theta, other_theta):
        if self.pendulum_bob == 1:
            return self.m1*self.g*self.L1*(1-cos(theta))  
        
        else:
            return self.m2*self.g*(self.L1*(1-cos(theta))+self.L2*(1-cos(other_theta)))
        
    