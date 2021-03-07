from main import System, save_data, fetch_data, pi, make_animation, make_plot, find_phase_space

pendulum = System(1, 1, 1, 1, 2.74, 0, 0.128, 0, 500, 10)
pend_data = pendulum.make_data()
#make_animation(pend_data)
make_plot(pend_data, plot_energy=True)

#find_phase_space(1, 1, 1, 1, 50, 100, phasespace_step_size=50)







