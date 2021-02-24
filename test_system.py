from system import System, Pendulum, pi
import pytest
from math import isclose

def test_system_init():
    length_1, mass_1, length_2, mass_2 = 1, 5, 2, 4
    initial_angular_position_1 = pi/2
    initial_angular_velocity_1 = 0
    initial_angular_position_2 = pi/4
    initial_angular_velocity_2 = 0
    steps, time =  1000, 100
    system = System(length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps, time)

    assert system.g == 9.8
    assert system.n == steps
    assert len(system.time) == steps+1
    assert isinstance(system.p1, Pendulum)
    assert isinstance(system.p2, Pendulum)
    assert system.p1.angular_acceleration[0] == system.p1.dwdt(initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2)
    assert system.p2.angular_acceleration[0] == system.p2.dwdt(initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2)

@pytest.mark.parametrize('angle, expected', [
    (pi/2, pi/2), (pi, pi),
    (-1.5*pi, pi/2), (2*pi, 0),
    (-pi, -pi), (-3*pi, -pi), 
    (13*pi, -pi) 
    ])
def test_normalise_angle(angle, expected):
    assert isclose(System.normalize_angle(angle), expected)


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



    
    
   