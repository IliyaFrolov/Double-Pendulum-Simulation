import matplotlib.pyplot as plt
from matplotlib import animation 
import pandas as pd
from system import System, np, pi

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
    'x position 1': pendulum.p1.x_position,
    'y position 1': pendulum.p1.y_position,
    'Angular velocity 1': pendulum.p1.angular_velocity,
    'Angular acceleration 1': pendulum.p1.angular_acceleration,
    'Angular position 2': pendulum.p2.angular_position,
    'x position 2': pendulum.p2.x_position,
    'y position 2': pendulum.p2.y_position,
    'Angular velocity 2': pendulum.p2.angular_velocity,
    'Angular acceleration 2': pendulum.p2.angular_acceleration,
    'Kinetic energy': pendulum.kinetic_energy,
    'Potential energy': pendulum.potential_energy,
    'Total energy': pendulum.total_energy
    })

    return pendulum_data

pendulum = System(1, 1, 1, 1, pi/2, 0, pi/4, 0, 1000, 50)
pendulum.numerical_analysis()
pendulum_data = pandas_data(pendulum)

fig = plt.figure()
ax = plt.axes(xlim=(-2, 2), ylim=(-2, 2))
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

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=1000, interval=20, blit=True)


anim.save('basic_animation.gif')

plt.show()
    