Instructions:
The simulation is ran by creating objects and calling functions or methods from the "run.py" file. A step by step guide to running a double pendulum  and simple pendulum simulation is outlined below:

Double Pendulum
- Instantiate an object of the 'Sytem' class with the desired parameters.
- Call the 'make_'data' method from the 'System' class to return a Pandas Dataframe containing the simulation results.
- Choose to save the data locally by calling the 'save_data' function. Consequently, you can retreieve the saved data by calling the 'fetch_data' funciton returning the saved Pandas Dataframe.
- Choose to call 'make_plot_doublepend', 'make_plot_error' or 'make_animation' function to generate the desired plots and animations.
- Spereately, call 'make_phase_space' function to generate a phase space plot of the double pendulum.

Single Pendulum
- Instantiate an object of the 'Pendulum' class by calling the class method 'init_simple_pendulum' with the desired parameters. 
- Call the 'make_'data' method from the 'Pendulum' class to return a Pandas Dataframe containing the simulation results.
- Choose to save the data locally by calling the 'save_data' function. Consequently, you can retreieve the saved data by calling the 'fetch_data' funciton returning the saved Pandas Dataframe.
- Call 'make_plot_singlepend' function to generate plots for the simple pendulum.

Files:

"run.py" - The main file where objects are created and functions/methods are called to run simulations, produce plots, generate animations and save data.

"main.py" - This file contains functions to produce plots, animations and save data.

"system.py" - This file contains the main class 'System' that represents the double pendulum and runs the double pendulum simulation, returning the results in a Pandas Dataframe.

"pendulum.py" - This file contains the main class 'Pendulum' that represents a pendulum and runs the simple pendulum simulation, returning the results in a Pandas Dataframe.

"equations.py" - This file contains the equations for the double pendulum such as the equations of motion, kinetic energy, potential energy and a conversion from polar to Cartesian coordinates.
