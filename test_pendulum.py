from system import Pendulum, pi, cos, sin
from equations import *
import pytest

@pytest.fixture
def pendulum():
    length, mass = 1, 2
    initial_angular_position = pi/4
    initial_angular_velocity = pi/2
    steps, time = 100, 10
    return Pendulum.init_simple_pendulum(length, mass, initial_angular_position, initial_angular_velocity, steps, time)

def test_pendulum():
    length_1, mass_1, length_2, mass_2 = 1, 5, 3, 2
    initial_angular_position_1 = pi/2
    initial_angular_velocity_1 = 0
    initial_angular_position_2 = pi/4
    initial_angular_velocity_2 = 0
    steps =  1000
    pendulum_1 = Pendulum(1, length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps)
    pendulum_2 = Pendulum(2, length_2, mass_2, length_1, mass_1, initial_angular_position_2, initial_angular_velocity_2, initial_angular_position_1, initial_angular_velocity_1, steps)

    assert isinstance(pendulum_1.dwdt, Acceleration)
    assert isinstance(pendulum_1.K, KineticEnergy)
    assert isinstance(pendulum_1.U, PotentialEnergy)
    assert isinstance(pendulum_2.convert, Polar_to_Cartesian)
    assert pendulum_1.pendulum_bob == 1
    assert pendulum_1.L == 1
    assert pendulum_1.m == 5
    assert pendulum_1.theta == pi/2
    assert pendulum_1.omega ==  0
    assert pendulum_1.n == 1000
    assert len(pendulum_2.time) == 1000
    assert len(pendulum_1.angular_position) == 1001
    assert len(pendulum_1.x_position) == 1001
    assert len(pendulum_1.y_position) == 1001
    assert len(pendulum_1.angular_velocity) == 1001
    assert len(pendulum_1.angular_acceleration) == 1001
    assert pendulum_1.angular_position[0] == pi/2
    assert pendulum_1.angular_velocity[0] == 0
    assert pendulum_1.angular_acceleration[0] == (-9.8*(2*mass_1+mass_2)*sin(initial_angular_position_1)-mass_2*9.8*sin(initial_angular_position_1-2*initial_angular_position_2)-2*sin(initial_angular_position_1-initial_angular_position_2)*mass_2*(initial_angular_velocity_2**2*length_2+initial_angular_velocity_1**2*length_1*cos(initial_angular_position_1-initial_angular_position_2))) / (length_1*(2*mass_1+mass_2-mass_2*cos(2*initial_angular_position_1-2*initial_angular_position_2)))
    
    assert pendulum_2.pendulum_bob == 2
    assert pendulum_2.L == 3
    assert pendulum_2.m == 2
    assert pendulum_2.angular_position[0] == pi/4
    assert pendulum_2.angular_velocity[0] == 0
    assert pendulum_2.angular_acceleration[0] == (2*sin(initial_angular_position_1-initial_angular_position_2)*(initial_angular_velocity_1**2*length_1*(mass_1+mass_2)+9.8*(mass_1+mass_2)*cos(initial_angular_position_1)+initial_angular_velocity_2**2*length_2*mass_2*cos(initial_angular_position_1-initial_angular_position_2))) / (length_2*(2*mass_1+mass_2-mass_2*cos(2*initial_angular_position_1-2*initial_angular_position_2)))

def test_init_simple_pendulum(pendulum):
    assert pendulum.pendulum_bob == 1
    assert pendulum.L == 1
    assert pendulum.m == 2
    assert pendulum.theta == pi/4
    assert pendulum.omega == pi/2
    assert pendulum.n == 100
    assert len(pendulum.time) == 100

@pytest.mark.parametrize('t, theta, omega', [
    ([1, 2], pi/2, pi),
    ([5, 7], 3*pi, pi/4),
    ([0, 4.5], 2*pi, 0)
])
def test_model(pendulum, t, theta, omega):
    dwdt = -9.8/pendulum.L*sin(theta) 
    assert pendulum.model(t, [theta, omega]) == [omega, dwdt] 