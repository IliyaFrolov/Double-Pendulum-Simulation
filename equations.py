from numpy import sin, cos, pi

class EquationTerms():

    def __init__(self, length_1, mass_1, length_2, mass_2, pendulum_bob):
        self.g = 9.8
        self.L1 = length_1
        self.m1 = mass_1
        self.L2 = length_2
        self.m2 = mass_2
        self.pendulum_bob = pendulum_bob

        if pendulum_bob != 1 and pendulum_bob != 2:
            raise Exception('Parameter "pendulum_bob" must be either 1 or 2.')

class Acceleration(EquationTerms):

    def __init__(self, length_1, mass_1, length_2, mass_2, pendulum_bob):
        super().__init__(length_1, mass_1, length_2, mass_2, pendulum_bob)
    
    def __call__(self, theta_1, omega_1, theta_2, omega_2):
        if self.pendulum_bob == 1:
            return (-self.g*(2*self.m1+self.m2)*sin(theta_1)-self.m2*self.g*sin(theta_1-2*theta_2)-2*sin(theta_1-theta_2)*self.m2*(omega_2**2*self.L2+omega_1**2*self.L1*cos(theta_1-theta_2))) / (self.L1*(2*self.m1+self.m2-self.m2*cos(2*theta_1-2*theta_2)))

        elif self.pendulum_bob == 2:
            return (2*sin(theta_1-theta_2)*(omega_1**2*self.L1*(self.m1+self.m2)+self.g*(self.m1+self.m2)*cos(theta_1)+omega_2**2*self.L2*self.m2*cos(theta_1-theta_2))) / (self.L2*(2*self.m1+self.m2-self.m2*cos(2*theta_1-2*theta_2)))
        

class KineticEnergy(EquationTerms):

    def __init__(self, length_1, mass_1, length_2, mass_2, pendulum_bob):
        super().__init__(length_1, mass_1, length_2, mass_2, pendulum_bob)
    
    def __call__(self, theta_1, omega_1, theta_2, omega_2):
        if self.pendulum_bob == 1:
            return self.m1/2*self.L1**2*omega_1**2
        
        elif self.pendulum_bob == 2:
            return self.m1/2*(self.L1**2*omega_1**2+self.L2**2*omega_2**2+2*self.L1*self.L2*omega_1*omega_2*cos(theta_1-theta_2))

class PotentialEnergy(EquationTerms):

    def __init__(self, length_1, mass_1, length_2, mass_2, pendulum_bob):
        super().__init__(length_1, mass_1, length_2, mass_2, pendulum_bob)
         
    def __call__(self, theta_1, theta_2):
        if self.pendulum_bob == 1:
            return self.m1*self.g*self.L1*(1-cos(theta_1))  
        
        elif self.pendulum_bob == 2:
            return self.m1*self.g*(self.L1*(1-cos(theta_1))+self.L2*(1-cos(theta_2)))
        
    