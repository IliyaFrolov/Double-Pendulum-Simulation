from main import System, save_data, fetch_data, pi, make_animation, make_plot, find_phase_space
from pendulum import Pendulum

#pendulum = System(1, 2, 3, 4, 3*pi, 0, 2*pi, 0, 1000, 50)
#pend_data = pendulum.make_data(method='Radau')
#make_plot(pend_data, plot_energy=True)
#make_animation(pendulum, pend_data)

#find_phase_space(1, 1, 1, 1, 100, 10, angle_step_size=10)

pendulum = Pendulum.init_simple_pendulum(1, 1, pi/2, 0, 1000, 10)
pendulum.make_simple_pendulum(method='Radau')




