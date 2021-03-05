from system import System, save_data, fetch_data, pi, make_animation, make_plot, find_phase_space

pendulum = System(1, 2, 3, 4, pi/2, 0, pi/4, 0, 500, 10)
print(pendulum.p1)
print(pendulum.p2)
#pend_data = pendulum.make_data()
#make_animation(pend_data)
#make_plot(pend_data, plot_energy=True)

#find_phase_space(1, 1, 1, 1, 500, 50, phasespace_step_size=50, save=False)







