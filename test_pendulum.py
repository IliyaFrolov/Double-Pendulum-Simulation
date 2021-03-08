from system import Pendulum, pi, cos, sin

def test_pendulum():
    length_1 = 1
    mass_1 = 5
    length_2 = 3
    mass_2 = 2
    initial_angular_position_1 = pi/2
    initial_angular_velocity_1 = 0
    initial_angular_position_2 = pi/4
    initial_angular_velocity_2 = 0
    steps =  1000
    pendulum_1 = Pendulum(1, length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps)
    pendulum_2 = Pendulum(2, length_2, mass_2, length_1, mass_1, initial_angular_position_2, initial_angular_velocity_2, initial_angular_position_1, initial_angular_velocity_1, steps)

    assert pendulum_1.L == 1
    assert pendulum_1.m == 5
    assert len(pendulum_1.angular_position) == 1001
    assert len(pendulum_1.x_position) == 1001
    assert len(pendulum_1.y_position) == 1001
    assert len(pendulum_1.angular_velocity) == 1001
    assert len(pendulum_1.angular_acceleration) == 1001
    assert pendulum_1.angular_position[0] == pi/2
    assert pendulum_1.angular_velocity[0] == 0
    assert pendulum_1.angular_acceleration[0] == (-9.8*(2*mass_1+mass_2)*sin(initial_angular_position_1)-mass_2*9.8*sin(initial_angular_position_1-2*initial_angular_position_2)-2*sin(initial_angular_position_1-initial_angular_position_2)*mass_2*(initial_angular_velocity_2**2*length_2+initial_angular_velocity_1**2*length_1*cos(initial_angular_position_1-initial_angular_position_2))) / (length_1*(2*mass_1+mass_2-mass_2*cos(2*initial_angular_position_1-2*initial_angular_position_2)))
    
    assert pendulum_2.L == 3
    assert pendulum_2.m == 2
    assert len(pendulum_2.angular_position) == 1001
    assert len(pendulum_2.x_position) == 1001
    assert len(pendulum_2.y_position) == 1001
    assert len(pendulum_2.angular_velocity) == 1001
    assert len(pendulum_2.angular_acceleration) == 1001
    assert pendulum_2.angular_position[0] == pi/4
    assert pendulum_2.angular_velocity[0] == 0
    assert pendulum_2.angular_acceleration[0] == (2*sin(initial_angular_position_1-initial_angular_position_2)*(initial_angular_velocity_1**2*length_1*(mass_1+mass_2)+9.8*(mass_1+mass_2)*cos(initial_angular_position_1)+initial_angular_velocity_2**2*length_2*mass_2*cos(initial_angular_position_1-initial_angular_position_2))) / (length_2*(2*mass_1+mass_2-mass_2*cos(2*initial_angular_position_1-2*initial_angular_position_2)))