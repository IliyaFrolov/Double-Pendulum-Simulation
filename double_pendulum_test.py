from double_pendulum import Pendulum, System
import numpy as np

def test_system_init():
    system = System(1, 5, 2, 4, np.pi/2, 0, np.pi/4, 0, 1000, 100)

    assert system.g == 9.8
    assert len(system.time) == 1001
    assert isinstance(system.pendulum_1, Pendulum)
    assert isinstance(system.pendulum_2, Pendulum)

    assert system.pendulum_1.L == 1
    assert system.pendulum_1.m == 5 
    assert len(system.pendulum_1.angular_position) == 1001
    assert len(system.pendulum_1.angular_velocity) == 1001
    assert system.pendulum_1.angular_position[0] == np.pi/2
    assert system.pendulum_1.angular_velocity[0] == 0

    assert system.pendulum_2.L == 2
    assert system.pendulum_2.m == 4
    assert len(system.pendulum_2.angular_position) == 1001
    assert len(system.pendulum_2.angular_velocity) == 1001
    assert system.pendulum_2.angular_position[0] == np.pi/4
    assert system.pendulum_2.angular_velocity[0] == 0




    
    
   