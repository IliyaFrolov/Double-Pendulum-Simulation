from system import System, save_data, fetch_data, pi, make_animation, plot_data

pendulum = System(1, 1, 1, 1, pi/2, 0, pi/4, 0, 1000, 50)
#pendulum_data = fetch_data('pendulum data')
data = pendulum.make_data()
plot_data(data, True)




