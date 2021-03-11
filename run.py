from main import System, save_data, fetch_data, pi, make_animation, make_plot, find_phase_space
from pendulum import Pendulum

pendulum = System(1, 2, 3, 4, 3*pi, 0, 2*pi, 0, 1000, 50)
pend_data = pendulum.make_data(method='Radau')
make_plot(pend_data, plot_energy=True)
make_animation(pendulum, pend_data)

#find_phase_space(1, 1, 1, 1, 5000, 500, angle_step_size=200, file_name="1111 S5000 T500 P50")
#find_phase_space(1, 1, 1, 1, 500, 300, angle_step_size=50)

#pendulum = Pendulum.init_simple_pendulum(1, 2, pi/2, 0, 1000, 50)
#pendulum.make_simple_pendulum()




