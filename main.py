import matplotlib.pyplot as plt
import pandas as pd
from system import System, np

def plotting(pendulum, plot_energy=False):
    plt.plot(pendulum['Time'], pendulum['Angular position 1'], 'r-', label='angular_position_1')
    plt.plot(pendulum['Time'], pendulum['Angular position 2'], 'g-', label='angular_position_2')
    plt.legend()
    plt.show()

    if plot_energy:
        plt.plot(pendulum['Time'], pendulum['Kinetic energy'], 'b-', label='kinetic')
        plt.plot(pendulum['Time'], pendulum['Potential energy'], 'g-', label='potential')
        plt.plot(pendulum['Time'], pendulum['Total energy'], 'r-,', label='total energy')
        plt.legend()
        plt.show() 

def pandas_data(pendulum):
    pendulum_data = pd.DataFrame({
    'Time': pendulum.time,
    'Angular position 1': pendulum.p1.angular_position,
    'Angular velocity 1': pendulum.p1.angular_velocity,
    'Angular acceleration 1': pendulum.p1.angular_acceleration,
    'Angular position 2': pendulum.p2.angular_position,
    'Angular velocity 2': pendulum.p2.angular_velocity,
    'Angular acceleration 2': pendulum.p2.angular_acceleration,
    'Kinetic energy': pendulum.kinetic_energy,
    'Potential energy': pendulum.potential_energy,
    'Total energy': pendulum.total_energy
    })

    return pendulum_data

pendulum = System(1, 1, 1, 1, np.pi/2, 0, np.pi, 0, 1000, 50)
pendulum.numerical_analysis()
plotting(pandas_data(pendulum), True)
