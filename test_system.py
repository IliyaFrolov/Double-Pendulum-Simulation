from system import System, Pendulum, pi
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
    assert system.g == 9.8
    assert system.n == 1000
    assert len(system.time) == 1001
    assert isinstance(system.p1, Pendulum)
    assert isinstance(system.p2, Pendulum)
    assert system.p1.angular_acceleration[0] == system.p1.dwdt(pi/2, 0, pi/4, 0)
    assert system.p2.angular_acceleration[0] == system.p2.dwdt(pi/2, 0, pi/4, 0)

def test_model(system):
    assert system.model([2.5, 3.5], [pi/2, pi, pi/4, -pi]) == [pi, system.p1.dwdt(pi/2, pi, pi/4, -pi), -pi, system.p2.dwdt(pi/2, pi, pi/4, -pi)]

@pytest.mark.parametrize('angle, expected', [
    (pi/2, pi/2), (pi, pi),
    (-1.5*pi, pi/2), (2*pi, 0),
    (-pi, -pi), (-3*pi, -pi), 
    (13*pi, -pi) 
    ])
def test_normalise_angle(angle, expected):
    assert isclose(System.normalize_angle(angle), expected)







    
    
   