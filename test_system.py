from system import System, Pendulum
from equations import *
import pytest
from math import isclose

@pytest.fixture
def system():
    length_1, mass_1, length_2, mass_2 = 1, 5, 2, 4
    initial_angular_position_1 = pi/2
    initial_angular_velocity_1 = 0
    initial_angular_position_2 = pi/4
    initial_angular_velocity_2 = 0
    steps, time =  1000, 100
    return System(length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps, time)

def test_system_init(system):
    assert isinstance(system.p1, Pendulum)
    assert isinstance(system.p2, Pendulum)
    assert system.g == g
    assert system.n == 1000
    assert system.flip_time == None
    assert system.has_flipped == False
    assert len(system.time) == 1001
    assert len(system.kinetic_energy) == 1001
    assert len(system.potential_energy) == 1001
    assert len(system.total_energy) == 1001
    assert system.kinetic_energy[0] == system.p1.m/2*system.p1.L**2*system.p1.angular_velocity[0]**2 + system.p2.m/2*(system.p1.L**2*system.p1.angular_velocity[0]**2+system.p2.L**2*system.p2.angular_velocity[0]**2+2*system.p1.L*system.p2.L*system.p1.angular_velocity[0]*system.p2.angular_velocity[0]*cos(system.p1.angular_position[0]-system.p2.angular_position[0]))
    assert system.potential_energy[0] ==  system.p1.m*g*system.p1.L*(1-cos(system.p1.angular_position[0])) + system.p2.m*g*(system.p1.L*(1-cos(system.p1.angular_position[0]))+system.p2.L*(1-cos(system.p2.angular_position[0])))  
    assert system.total_energy[0] == system.kinetic_energy[0] + system.potential_energy[0]

    with pytest.raises(Exception):
        System(1, 0, 1, 1, pi/2, 0, pi/2, 0, 10, 10)
        System(1, 1, 0, 0, pi/2, 0, pi/2, 0, 10, 10)

@pytest.mark.parametrize('t, theta_1, omega_1, theta_2, omega_2', [
    ([2.5, 3.5], pi/2, pi, pi/4, -pi),
    ([4, 5], 5*pi, -pi, -3*pi/2, pi/2),
    ([0, 4], 6*pi, -pi/4, -2*pi, 8*pi/3)
])
def test_model(system, t, theta_1, omega_1, theta_2, omega_2):
    dwdt_1 = system.p1.dwdt(theta_1, omega_1, theta_2, omega_2)
    dwdt_2 = system.p2.dwdt(theta_2, omega_2, theta_1, omega_1)
    assert system.model(t, [theta_1, omega_1, theta_2, omega_2]) == [omega_1, dwdt_1, omega_2, dwdt_2]

def test_check_flip(system):
    theta_1, theta_2, i = pi/2, pi/4 , 2
    system.has_flipped = False
    system.check_flip(theta_1, theta_2, i)
    assert system.has_flipped == False

    theta_1, theta_2, i = pi, -pi , 2
    system.has_flipped = False
    system.check_flip(theta_1, theta_2, i)
    assert system.has_flipped == False

    theta_1, theta_2, i = 2*pi, pi/4 , 2
    system.has_flipped = False
    system.check_flip(theta_1, theta_2, i)
    assert system.has_flipped == True


@pytest.mark.parametrize('angle, expected', [
    (pi/2, pi/2), (pi, pi),
    (-1.5*pi, pi/2), (2*pi, 0),
    (-pi, -pi), (-3*pi, -pi), 
    (13*pi, -pi) 
    ])
def test_normalise_angle(angle, expected):
    assert isclose(System.normalise_angle(angle), expected)







    
    
   