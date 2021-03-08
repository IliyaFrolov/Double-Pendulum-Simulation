from main import System, save_data, fetch_data, pi, make_animation, make_plot, find_phase_space

pendulum = System(1, 0.1, 1, 1, pi/2, 0, pi/2, 0, 1500, 300)
pend_data = pendulum.make_data()
make_animation(pendulum, pend_data)
#make_plot(pend_data, plot_energy=True)

#find_phase_space(1, 1, 1, 1, 5000, 500, phasespace_step_size=50, file_name="1111 S5000 T500 P50")
#find_phase_space(1, 1, 1, 1, 500, 300, phasespace_step_size=50)






