from main import System, save_data, fetch_data, pi, make_animation, make_plot_double, make_plot_single, find_phase_space
from pendulum import Pendulum

pendulum = System(1, 1, 1, 1, pi/6, 0, pi/2, 0, 5000, 500)
#pend_data_1 = pendulum.make_data(method='RK45')
#pend_data_2 = pendulum.make_data(method='Radau')
#pend_data_3 = pendulum.make_data(method='BDF')
#make_plot(pend_data_1)
#make_plot(pend_data_2)
#make_plot(pend_data_3)

pendulums = [Pendulum.init_simple_pendulum(1, 1, pi/6, 0, 10, 10), Pendulum.init_simple_pendulum(1, 1, pi/6, 0, 100, 10), Pendulum.init_simple_pendulum(1, 1, pi/6, 0, 1000, 10)]

for pendulum, i in zip(pendulums, [10, 100, 1000]):
    pend_data_1 = pendulum.make_data(method='RK45')
    pend_data_2 = pendulum.make_data(method='Radau')
    make_plot_single(pendulum, pend_data_1, file_name=f'single RK45 S={i} T=10')
    make_plot_single(pendulum, pend_data_2, file_name=f'single Radau S={i} T=10')




