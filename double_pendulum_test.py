from double_pendulum import Pendulum
import numpy as np

def test_pendulum_init():
    pendulum = Pendulum(1, 1, np.pi/2, 0, 100, 1000)
    
    assert hasattr(pendulum, 'L') == True
    assert hasattr(pendulum, 'm') == True
    assert hasattr(pendulum, 'g') == True
    assert hasattr(pendulum, 't') == True
    assert hasattr(pendulum, 'n') == True
   