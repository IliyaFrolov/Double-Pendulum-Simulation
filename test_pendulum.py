from system import Pendulum, pi

def test_pendulum():
    length = 1
    mass = 5
    initial_angular_position = pi/2
    initial_angular_velocity = 0
    steps =  1000
    pendulum = Pendulum(length, mass, initial_angular_position, initial_angular_velocity, steps)

    assert pendulum.L == length
    assert pendulum.m == mass 
    assert len(pendulum.angular_position) == steps+1
    assert len(pendulum.x_position) == steps+1
    assert len(pendulum.y_position) == steps+1
    assert len(pendulum.angular_velocity) == steps+1
    assert len(pendulum.angular_acceleration) == steps+1
    assert pendulum.angular_position[0] == initial_angular_position
    assert pendulum.angular_velocity[0] == initial_angular_velocity
    