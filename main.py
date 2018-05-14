"""
Run the N body simulation

George Fryer 2018
"""


from orr import read_file
from nbody import NBody



params = read_file()
NBody(params[2:]).run(time=params[0], dt=params[1], force_video=False)
