from system import System, save_data, fetch_data, pi, make_animation, make_plot, find_phase_space

#pendulum = System(1, 1, 1, 1, pi, 0, pi, 0, 1000, 100)
#pend_data = pendulum.make_data()
#make_plot(pend_data)

find_phase_space(1, 2, 1, 2, 100, 500, phasespace_step_size=50)







