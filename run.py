from system import System, save_data, fetch_data, pi, make_animation, make_plot

pendulum = System(2, 3, 4, 1, pi/2, 0, pi/2, 0, 500, 50)
pend_data = pendulum.make_data(method='DOP853')
make_plot(pend_data, plot_energy=True)








