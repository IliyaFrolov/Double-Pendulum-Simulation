import pandas as pd
from equations import *
from pendulum import Pendulum, np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import animation 

class System():
    '''
    A class to represent the Double Pendulum.
    ...
    Attributes
    ----------
    p1 : class
        Instance of the Pendulum class representing the top pendulum.
    p2 : class
        Instance of the Pendulum class representing the bottom pendulum.
    g : int
        Gravitational acceleration constant.
    n : int
        Number of steps in the simulation.
    flip_time : int
        The time it takes for one of the pendulum bobs to flip. Is set to None by default if there's no flip.
    has_flipped : boolean
        Is set to True if one of the pendulum bobs has flipped, is False otherwise.
    time : numpy array
        Stores the time of the simulation at each step.
    kinetic_energy : numpy array
        Stores the kinetic energy of the system at each step.
    potential_energy : numpy array
        Stores the potential energy of the system at each step.
    total_energy : numpy array
        Stores the total energy of the system at each step.
    
    Methods
    ----------
    model(t, initial_conditions)
        Function is used as an input into solve_ivp to get the angular position and velocity of each pendulum bob at the next step.
    run_simulation(method)
        Uses solve_ivp within a "for loop" to compute and store the position and velocity (angular and linear) for each pendulum bob at each step. 
    check_flip(theta_1, theta_2, i)
        Checks if one of the pendulum bobs has flipped, if so, stores the time of flip.
    make_data(method)
        Calls self.run_simulation and returns all the computed data as a Pandas Dataframe.
    normalise_angle(anlge)
        Normalises the angle between -pi and pi. Returns the angle.
    '''
    def __init__(self, length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps, time):
        self.p1 = Pendulum(length_1, mass_1, length_2, mass_2, initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2, steps, 1)
        self.p2 = Pendulum(length_2, mass_2, length_1, mass_1, initial_angular_position_2, initial_angular_velocity_2, initial_angular_position_1, initial_angular_velocity_1, steps, 2)

        self.g = 9.8
        self.n = steps
        self.flip_time = None
        self.has_flipped = False
        self.time = np.linspace(0, time, steps+1)
        self.kinetic_energy = np.zeros(steps+1)
        self.potential_energy = np.zeros(steps+1)
        self.total_energy = np.zeros(steps+1)

        self.kinetic_energy[0] = self.p1.K(initial_angular_position_1, initial_angular_velocity_1, initial_angular_position_2, initial_angular_velocity_2) + self.p2.K(initial_angular_position_2, initial_angular_velocity_2, initial_angular_position_1, initial_angular_velocity_1)
        self.potential_energy[0] = self.p1.U(initial_angular_position_1, initial_angular_position_2) + self.p2.U(initial_angular_position_2, initial_angular_position_1)
        self.total_energy[0] = self.kinetic_energy[0] + self.potential_energy[0]

    def model(self, t, initial_conditions):
        theta_1 = initial_conditions[0]
        omega_1 = initial_conditions[1]
        theta_2 = initial_conditions[2]
        omega_2 = initial_conditions[3]
        dwdt_1 = self.p1.dwdt(theta_1, omega_1, theta_2, omega_2)
        dwdt_2 = self.p2.dwdt(theta_2, omega_2, theta_1, omega_1)

        return [omega_1, dwdt_1, omega_2, dwdt_2]    
    
    def run_simulation(self, method):
        for i in range(1, self.n+1):
            output = solve_ivp(
                self.model, 
                self.time[i-1: i+1],
                [self.p1.angular_position[i-1], self.p1.angular_velocity[i-1], self.p2.angular_position[i-1], self.p2.angular_velocity[i-1]], 
                t_eval=np.linspace(self.time[i-1], self.time[i], 2), method=method
                ).y
            
            theta_1 = output[0][1]
            omega_1 = output[1][1]
            theta_2 = output[2][1]
            omega_2 = output[3][1]

            self.check_flip(theta_1, theta_2, i)

            self.p1.angular_position[i] = self.normalize_angle(theta_1)
            self.p1.x_position[i] = self.p1.convert(theta_1, 'x')
            self.p1.y_position[i] = self.p1.convert(theta_1, 'y')
            self.p1.angular_velocity[i] = omega_1
            self.p1.angular_acceleration[i] = self.p1.dwdt(theta_1, omega_1, theta_2, omega_2)

            self.p2.angular_position[i] = self.normalize_angle(theta_2)
            self.p2.x_position[i] = self.p1.x_position[i] + self.p2.convert(theta_2, 'x')
            self.p2.y_position[i] = self.p1.y_position[i] + self.p2.convert(theta_2, 'y') 
            self.p2.angular_velocity[i] = omega_2
            self.p2.angular_acceleration[i] = self.p2.dwdt(theta_2, omega_2, theta_1, omega_1)
            
            self.kinetic_energy[i] =  self.p1.K(theta_1, omega_1, theta_2, omega_2) + self.p2.K(theta_2, omega_2, theta_1, omega_1)
            self.potential_energy[i] = self.p1.U(theta_1, theta_2) + self.p2.U(theta_2, theta_1)
            self.total_energy[i] = self.kinetic_energy[i] + self.potential_energy[i] 
        
    def check_flip(self, theta_1, theta_2, i):
        if (theta_1 > pi or theta_1 < -pi or theta_2 > pi or theta_2 < -pi) and not self.has_flipped:
            self.flip_time = self.time[i]
            self.has_flipped = True

    def make_data(self, method='DOP853'):
        self.run_simulation(method)

        return pd.DataFrame({
        'Time': self.time,
        'Angular position 1': self.p1.angular_position,
        'x position 1': self.p1.x_position,
        'y position 1': self.p1.y_position,
        'Angular velocity 1': self.p1.angular_velocity,
        'Angular acceleration 1': self.p1.angular_acceleration,
        'Angular position 2': self.p2.angular_position,
        'x position 2': self.p2.x_position,
        'y position 2': self.p2.y_position,
        'Angular velocity 2': self.p2.angular_velocity,
        'Angular acceleration 2': self.p2.angular_acceleration,
        'Kinetic energy': self.kinetic_energy,
        'Potential energy': self.potential_energy,
        'Total energy': self.total_energy
        })
    
    @staticmethod
    def normalize_angle(angle):
        while angle > pi:
           angle -= 2*pi
        
        while angle < -pi:
           angle += 2*pi

        return angle  

def find_phase_space(length_1, mass_1, length_2, mass_2, steps, time, phasespace_step_size=10, save=True):
    initial_theta_1 = np.linspace(-pi, pi, phasespace_step_size)
    initial_theta_2 = np.linspace(-pi, pi, phasespace_step_size)
    angle_1, angle_2 = np.meshgrid(initial_theta_1, initial_theta_2)
    time_to_flip = np.zeros((phasespace_step_size, phasespace_step_size))

    for x, theta_1 in enumerate(initial_theta_1):
        for y, theta_2 in enumerate(initial_theta_2):
            pendulum = System(length_1, mass_1, length_2, mass_2, theta_1, 0, theta_2, 0, steps, time)
            pendulum.make_data()
            time_to_flip[x][y] = pendulum.flip_time 

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set(title='Phase space plot of the time it takes for the Double Pendulum to flip', ylabel='Initial displacement of top Pendulum (radians)', xlabel='Initial displacement of bottom Pendulum (radians)')
    cp = ax.contourf(angle_1, angle_2, time_to_flip, levels=[0, 10/np.sqrt(9.8), 100/np.sqrt(9.8), 1000/np.sqrt(9.8)], extend='max')
    fig.colorbar(cp) 
    plt.show()

    if save:
        file_name = input('Enter file name: ')
        plt.savefig(rf'C:\Users\iliya\OneDrive\phys_389\graphs\{file_name}.png')

def make_plot(pendulum_data, plot_energy=False, file_name=None):
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

    if plot_energy:
        fig2 = plt.figure()
        energy_ax = fig2.add_subplot(111)
        energy_ax.set(title='Energy of the Double Pendulum over time.', ylabel='Energy (J)', xlabel='Time (s)')
        energy_ax.plot(pendulum_data['Time'], pendulum_data['Kinetic energy'], 'b-', label='Kinetic energy')
        energy_ax.plot(pendulum_data['Time'], pendulum_data['Potential energy'], 'g-', label='Potential energy')
        energy_ax.plot(pendulum_data['Time'], pendulum_data['Total energy'], 'r-,', label='Total energy')
        plt.legend(loc='upper left')
    
    if file_name:
        fig1.savefig(rf'C:\Users\Iliya Frolov\OneDrive\phys_389\graphs\{file_name}.png')
        fig2.savefig(rf'C:\Users\Iliya Frolov\OneDrive\phys_389\graphs\{file_name}.png')
    
    else:
        plt.show()

def make_animation(pendulum_data, file_name=None):
    fig = plt.figure()
    ax = plt.axes(xlim=(-2, 2), ylim=(-2, 2))
    ax.set(title='Animation of the Double Pendulum', ylabel='Angular displacement (radians)', xlabel='Angular displacement (radians)')
    line, = ax.plot([], [], 'o-')

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

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=10000, interval=10, blit=True)
    plt.show()

    if file_name:
        anim.save(rf'C:\Users\Iliya Frolov\OneDrive\phys_389\animations\{file_name}.gif')

def save_data(pendulum_data, file_name):
    pendulum_data.to_pickle(rf'C:\Users\Iliya Frolov\OneDrive\phys_389\modelling\{file_name}')

def fetch_data(file_name):
    return pd.read_pickle(rf'C:\Users\Iliya Frolov\OneDrive\phys_389\modelling\{file_name}')
    
