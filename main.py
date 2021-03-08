import matplotlib.pyplot as plt
import os
from matplotlib import animation 
from system import System, np, pd, pi

def find_phase_space(length_1, mass_1, length_2, mass_2, steps, time, phasespace_step_size=10, file_name=None):
    initial_theta_1 = np.linspace(-pi, pi, phasespace_step_size)
    initial_theta_2 = np.linspace(-pi, pi, phasespace_step_size)
    angle_1, angle_2 = np.meshgrid(initial_theta_1, initial_theta_2)
    time_to_flip = np.zeros((phasespace_step_size, phasespace_step_size))
    counter = phasespace_step_size**2
    units = np.sqrt(length_1/9.8)

    for x, theta_1 in enumerate(initial_theta_1):
        for y, theta_2 in enumerate(initial_theta_2):
            pendulum = System(length_1, mass_1, length_2, mass_2, theta_1, 0, theta_2, 0, steps, time)
            pendulum.run_simulation(method='DOP853', is_phase=True)
            time_to_flip[x][y] = pendulum.flip_time 
            counter -= 1
            print(f'{counter} iterations remaining')

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set(title='Phase space plot of the time it takes for the Double Pendulum to flip', ylabel='Initial displacement of top Pendulum (radians)', xlabel='Initial displacement of bottom Pendulum (radians)')
    cp = ax.contourf(angle_1, angle_2, time_to_flip, levels=[0, 5*units, 10*units, 50*units, 100*units, 500*units, 1000*units], cmap='jet' ,extend='max')
    fig.colorbar(cp) 
    plt.show()

    if file_name:
        plt.savefig(rf'{os.getcwd()}\graphs\{file_name}.png')

def make_plot(pendulum_data, plot_energy=False, save=False):
    fig1 = plt.figure()
    gs = fig1.add_gridspec(2, 2, hspace=0.5)
    position_ax = fig1.add_subplot(gs[1, :])
    position_ax.set(title='Angular displacements of the Double Pendulum over time.', ylabel='Angular displacement (radians)', xlabel='Time (s)')
    position_ax.plot(pendulum_data['Time'], pendulum_data['Angular position 1'], 'r-', label='Angular position 1')
    position_ax.plot(pendulum_data['Time'], pendulum_data['Angular position 2'], 'g-', label='Angular position 2')
    position_ax.legend(loc='upper left')

    velocity_ax = fig1.add_subplot(gs[0, 0])
    velocity_ax.set(title='Angular velocities of the Double Pendulum over time.', ylabel='Angular velocity (radians/s)', xlabel='Time (s)')
    velocity_ax.plot(pendulum_data['Time'], pendulum_data['Angular velocity 1'], 'r-', label='Angular velocity 1')
    velocity_ax.plot(pendulum_data['Time'], pendulum_data['Angular velocity 2'], 'g-', label='Angular velocity 2')
    velocity_ax.legend(loc='upper left')

    acceleration_ax = fig1.add_subplot(gs[0, 1])
    acceleration_ax.set(title='Angular acceleration of the Double Pendulum over time.', ylabel='Angular acceleration (radians/s^2)', xlabel='Time (s)')
    acceleration_ax.plot(pendulum_data['Time'], pendulum_data['Angular acceleration 1'], 'r-', label='Angular acceleration 1')
    acceleration_ax.plot(pendulum_data['Time'], pendulum_data['Angular acceleration 2'], 'g-', label='Angular acceleration 2')
    acceleration_ax.legend(loc='upper left')
    plt.show()

    if save:
        file_name = input('Enter file name: ')
        fig1.savefig(rf'{os.getcwd()}\graphs\{file_name}.png')

    if plot_energy:
        fig2 = plt.figure()
        energy_ax = fig2.add_subplot(111)
        energy_ax.set(title='Energy of the Double Pendulum over time.', ylabel='Energy (J)', xlabel='Time (s)')
        energy_ax.plot(pendulum_data['Time'], pendulum_data['Kinetic energy'], 'b-', label='Kinetic energy')
        energy_ax.plot(pendulum_data['Time'], pendulum_data['Potential energy'], 'g-', label='Potential energy')
        energy_ax.plot(pendulum_data['Time'], pendulum_data['Total energy'], 'r-,', label='Total energy')
        plt.legend(loc='upper left')
        plt.show()

        if save:
            file_name = input('Enter file name: ')
            fig2.savefig(rf'{os.getcwd()}\graphs\{file_name}.png')

def make_animation(pendulum, pendulum_data, file_name=None):
    fig = plt.figure()
    total_length = pendulum.p1.L + pendulum.p2.L
    ax = plt.axes(xlim=(-total_length-0.5, total_length+0.5), ylim=(-total_length-0.5, total_length+0.5))
    ax.set(title='Animation of the Double Pendulum', ylabel='Angular displacement (radians)', xlabel='Angular displacement (radians)')
    line, = ax.plot([], [], 'o-')
    frames = pendulum.n

    def init():
        line.set_data([], [])
        return line,

    def animate(i):
        x1 = pendulum_data['x position 1'][i]
        y1 = pendulum_data['y position 1'][i]
        x2 = pendulum_data['x position 2'][i]
        y2 = pendulum_data['y position 2'][i]
        line.set_data([0, x1, x2], [0, y1, y2])
    
        return line,

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=frames, interval=10, blit=True)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()

    if file_name:
        anim.save(rf'{os.getcwd()}\animations\{file_name}.gif')

def save_data(pendulum_data, file_name):
    pendulum_data.to_pickle(rf'{os.getcwd()}\saved_data\{file_name}')

def fetch_data(file_name):
    return pd.read_pickle(rf'{os.getcwd()}\saved_data\{file_name}')