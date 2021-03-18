from main import System, save_data, fetch_data, pi, make_animation, make_plot, find_phase_space
from pendulum import Pendulum

pendulum = System(1, 1, 1, 1, pi/2, 0, pi/2, 0, 500, 100)
pend_data_1 = pendulum.make_data(method='RK45')
pend_data_2 = pendulum.make_data(method='Radau')
pend_data_3 = pendulum.make_data(method='BDF')
make_plot(pend_data_1)
make_plot(pend_data_2)
make_plot(pend_data_3)

#pendulum = Pendulum.init_simple_pendulum(1, 1, pi/2, 0, 100, 10)
#pendulum.make_simple_pendulum(method='RK23')




