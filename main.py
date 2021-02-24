from system import System, save_data, fetch_data, pi, make_animation

pendulum = System(1, 1, 1, 1, pi/2, 0, pi/4, 0, 1000, 50)
pendulum_data = fetch_data('pendulum data')
make_animation(pendulum_data)


