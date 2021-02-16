from double_pendulum import Pendulum, System, Energy, KineticEnergy1, KineticEnergy2, PotentialEnergy1, PotentialEnergy2
import numpy as np
from numpy import pi, sin, cos
from math import isclose

def test_system_init():
    system = System(1, 5, 2, 4, pi/2, 0, pi/4, 0, 1000, 100)

    assert system.g == 9.8
    assert system.n == 1000
    assert len(system.time) == 1001
    assert isinstance(system.p1, Pendulum)
    assert isinstance(system.p2, Pendulum)

    assert system.p1.L == 1
    assert system.p1.m == 5 
    assert len(system.p1.angular_position) == 1001
    assert len(system.p1.angular_velocity) == 1001
    assert system.p1.angular_position[0] == pi/2
    assert system.p1.angular_velocity[0] == 0

    assert system.p2.L == 2
    assert system.p2.m == 4
    assert len(system.p2.angular_position) == 1001
    assert len(system.p2.angular_velocity) == 1001
    assert system.p2.angular_position[0] == pi/4
    assert system.p2.angular_velocity[0] == 0

def test_normalise_angle():
    assert isclose(System.normalize_angle(pi/2), pi/2)
    assert isclose(System.normalize_angle(pi), pi)
    assert isclose(System.normalize_angle(-1.5*pi), pi/2)
    assert isclose(System.normalize_angle(2*pi), 0)
    assert isclose(System.normalize_angle(-pi), -pi)
    assert isclose(System.normalize_angle(-3*pi), -pi)
    assert isclose(System.normalize_angle(13*pi), -pi)

def test_kineticenergy1():
    kinetic = KineticEnergy1(2, 3)

    assert kinetic.L == 2
    assert kinetic.m == 3
    assert kinetic(5) == 150

def test_kineticenergy2():
    length, length_2, mass, theta_1, omega_1, theta_2, omega_2 = 2, 3, 4, 5, 6, 7, 8
    kinetic = KineticEnergy2(length, length_2, mass)

    assert kinetic.L == 2
    assert kinetic.L2 == 3
    assert kinetic.m == 4
    assert kinetic(theta_1, omega_1, theta_2, omega_2) == mass/2*(length**2*omega_1**2+length_2**2*omega_2**2+2*length*length_2*omega_1*omega_2*cos(theta_1-theta_2))

'''
def test_numerical_analysis():
    system = System(1, 1, 1, 1, 1, 1, 1, 1, 2, 0.5)
    system.numerical_analysis()

def test_model():
    length_1 = 1
    mass_1 = 5
    length_2 = 2
    mass_2 = 4
    theta_1 = pi/4
    omega_1 = pi
    theta_2 = pi/2
    omega_2 = pi/8
    system = System(length_1, mass_1, length_2, mass_2, theta_1, omega_1, theta_2, omega_2, 1000, 100)
    output =  system.model([theta_1, omega_1, theta_2, omega_2], [1, 2]) 

    assert round(output[1]) == -3
'''



    
    
   