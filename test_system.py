from system import System, pi, cos, sin, Pendulum
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
    assert system.g == 9.8
    assert system.n == 1000
    assert system.flip_time == None
    assert system.has_flipped == False
    assert len(system.time) == 1001
    assert len(system.kinetic_energy) == 1001
    assert len(system.potential_energy) == 1001
    assert len(system.total_energy) == 1001
    assert system.kinetic_energy[0] == system.p1.m/2*system.p1.L**2*system.p1.angular_velocity[0]**2 + system.p2.m/2*(system.p1.L**2*system.p1.angular_velocity[0]**2+system.p2.L**2*system.p2.angular_velocity[0]**2+2*system.p1.L*system.p2.L*system.p1.angular_velocity[0]*system.p2.angular_velocity[0]*cos(system.p1.angular_position[0]-system.p2.angular_position[0]))
    assert system.potential_energy[0] ==  system.p1.m*9.8*system.p1.L*(1-cos(system.p1.angular_position[0])) + system.p2.m*9.8*(system.p1.L*(1-cos(system.p1.angular_position[0]))+system.p2.L*(1-cos(system.p2.angular_position[0])))  
    assert system.total_energy[0] == system.kinetic_energy[0] + system.potential_energy[0]

    with pytest.raises(Exception):
        System(1, 0, 1, 1, pi/2, 0, pi/2, 0, 10, 10)
        System(1, 1, 0, 0, pi/2, 0, pi/2, 0, 10, 10)

def test_model(system):
    assert system.model([2.5, 3.5], [pi/2, pi, pi/4, -pi]) == [pi, system.p1.dwdt(pi/2, pi, pi/4, -pi), -pi, system.p2.dwdt(pi/4, -pi, pi/2, pi)]

@pytest.mark.parametrize('angle, expected', [
    (pi/2, pi/2), (pi, pi),
    (-1.5*pi, pi/2), (2*pi, 0),
    (-pi, -pi), (-3*pi, -pi), 
    (13*pi, -pi) 
    ])
def test_normalise_angle(angle, expected):
    assert isclose(System.normalise_angle(angle), expected)







    
    
   