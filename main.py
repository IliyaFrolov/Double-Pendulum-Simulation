from system import System, save_data, fetch_data, pi, make_animation, make_plot

pendulum1 = System(1, 1, 1, 1, pi/2, 0, pi/4, 0, 50, 50)
make_plot(pendulum1.make_data(), plot_energy=True, save=True, file_name='Runge-Kutta 5th order - 50 steps')
make_plot(pendulum1.make_data(method='RK23'), plot_energy=True, save=True, file_name='Runge-Kutta 3rd order - 50 steps')
make_plot(pendulum1.make_data(method='DOP853'), plot_energy=True, save=True, file_name='Runge-Kutta 8th order - 50 steps')
make_plot(pendulum1.make_data(method='Radau'), plot_energy=True, save=True, file_name='Implicit Runge-kutta 5th order - 50 steps')
make_plot(pendulum1.make_data(method='BDF'), plot_energy=True, save=True, file_name='Implicit multi-step variable order 1-5 - 50 steps')
make_plot(pendulum1.make_data(method='LSODA'),plot_energy=True, save=True, file_name='multi-step with stiffness detection - 50 steps')

pendulum2 = System(1, 1, 1, 1, pi/2, 0, pi/4, 0, 500, 50)
make_plot(pendulum2.make_data(), plot_energy=True, save=True, file_name='Runge-Kutta 5th order - 500 steps')
make_plot(pendulum2.make_data(method='RK23'), plot_energy=True, save=True, file_name='Runge-Kutta 3rd order - 500 steps')
make_plot(pendulum2.make_data(method='DOP853'), plot_energy=True, save=True, file_name='Runge-Kutta 8th order - 500 steps')
make_plot(pendulum2.make_data(method='Radau'), plot_energy=True, save=True, file_name='Implicit Runge-kutta 5th order - 500 steps')
make_plot(pendulum2.make_data(method='BDF'), plot_energy=True, save=True, file_name='Implicit multi-step variable order 1-5 - 500 steps')
make_plot(pendulum2.make_data(method='LSODA'), plot_energy=True, save=True, file_name='multi-step with stiffness detection - 500 steps')









