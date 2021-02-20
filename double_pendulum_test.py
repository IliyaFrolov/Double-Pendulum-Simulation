from double_pendulum import Pendulum, System, EquationTerms, KineticEnergy, PotentialEnergy, Acceleration
import numpy as np
import pytest
from numpy import pi, sin, cos
from math import isclose

def test_system_init():
    length_1, mass_1, length_2, mass_2 = 1, 5, 2, 4
    inital_angluar_position_1 = pi/2
    inital_angluar_velocity_1 = 0
    inital_angluar_position_2 = pi/4
    inital_angluar_velocity_2 = 0
    steps, time =  1000, 100
    system = System(length_1, mass_1, length_2, mass_2, inital_angluar_position_1, inital_angluar_velocity_1, inital_angluar_position_2, inital_angluar_velocity_2, steps, time)

    assert system.g == 9.8
    assert system.n == steps
    assert len(system.time) == steps+1
    assert isinstance(system.p1, Pendulum)
    assert isinstance(system.p2, Pendulum)

    assert system.p1.L == length_1
    assert system.p1.m == mass_1 
    assert len(system.p1.angular_position) == steps+1
    assert len(system.p1.angular_velocity) == steps+1
    assert len(system.p1.angular_acceleration) == steps+1
    assert system.p1.angular_position[0] == inital_angluar_position_1
    assert system.p1.angular_velocity[0] == inital_angluar_velocity_1
    assert system.p1.angular_acceleration[0] == system.p1.dwdt(inital_angluar_position_1, inital_angluar_velocity_1, inital_angluar_position_2, inital_angluar_velocity_2)


    assert system.p2.L == length_2
    assert system.p2.m == mass_2
    assert len(system.p2.angular_position) == steps+1
    assert len(system.p2.angular_velocity) == steps+1
    assert len(system.p2.angular_acceleration) == steps+1
    assert system.p2.angular_position[0] == inital_angluar_position_2
    assert system.p2.angular_velocity[0] == inital_angluar_velocity_2
    assert system.p2.angular_acceleration[0] == system.p2.dwdt(inital_angluar_position_1, inital_angluar_velocity_1, inital_angluar_position_2, inital_angluar_velocity_2)

@pytest.mark.parametrize('angle, expected', [
    (pi/2, pi/2), (pi, pi),
    (-1.5*pi, pi/2), (2*pi, 0),
    (-pi, -pi), (-3*pi, -pi), 
    (13*pi, -pi) 
    ])
def test_normalise_angle(angle, expected):
    assert isclose(System.normalize_angle(angle), expected)

@pytest.mark.parametrize('length_1, mass_1, length_2, mass_2, theta_1, omega_1, theta_2, omega_2', [
    (2, 3, 4, 5, pi, pi/2, pi/4, 2*pi), 
    (4, 4, 3, 1, 5*pi, 0, pi/2, 6),
    (0, 1, 10, 4, 2, 3*pi, 0, 2)
])
class TestEnergy():

    def test_kineticenergy(self, length_1, mass_1, length_2, mass_2, theta_1, omega_1, theta_2, omega_2):
        kinetic_1 = KineticEnergy(length_1, mass_1, length_2, mass_2, 1)
        kinetic_2 = KineticEnergy(length_1, mass_1, length_2, mass_2, 2)

        assert kinetic_1.L1 == length_1
        assert kinetic_1.m1 == mass_1
        assert kinetic_1.L2 == length_2
        assert kinetic_1.m2 == mass_2
        assert kinetic_1(theta_1, omega_1, theta_2, omega_2) == mass_1/2*length_1**2*omega_1**2
        assert kinetic_2(theta_1, omega_1, theta_2, omega_2) == mass_1/2*(length_1**2*omega_1**2+length_2**2*omega_2**2+2*length_1*length_2*omega_1*omega_2*cos(theta_1-theta_2))

    def test_potentialenergy(self, length_1, mass_1, length_2, mass_2, theta_1, omega_1, theta_2, omega_2):
        potential_1 = PotentialEnergy(length_1, mass_1, length_2, mass_2, 1)
        potential_2 = PotentialEnergy(length_1, mass_1, length_2, mass_2, 2)

        assert potential_1.L1 == length_1
        assert potential_1.m1 == mass_1
        assert potential_1.L2 == length_2
        assert potential_1.m2 == mass_2
        assert potential_1(theta_1, theta_2) == mass_1*9.8*length_1*(1-cos(theta_1))
        assert potential_2(theta_1, theta_2) == mass_1*9.8*(length_1*(1-cos(theta_1))+length_2*(1-cos(theta_2)))


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



    
    
   