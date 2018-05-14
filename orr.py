"""
Reads initial state from file

George Fryer 2018
"""


import body



def read_file(filename="initial.orr"):
    """
    Accepts name of file as string and returns parameters of the simulation
    based on the file contents. Returns an array: [time, dt, b1, b2, ...]

    The file must be in the following format (.orr):

    time
    dt

    obj0_name
    obj0_colour,obj0_marker
    obj0_mass_kg
    obj0_pos_x_km
    obj0_pos_y_km
    obj0_vel_x_km/s
    obj0_vel_y_km/s

    obj1_name
    obj1_colour,obj1_marker
    obj1_mass_kg
    obj1_pos_x_km
    obj1_pos_y_km
    obj1_vel_x_km/s
    obj1_vel_y_km/s

    ...

    A .orr file must end with a blank line. The time and dt parameters will
    be converted to integers, and any fractional part truncated

    """

    params = []
    with open(filename, mode="r") as f:
        if not f:
            raise IOError("Could not open input file")
        for _ in range(2):
            params.append(int(float(f.readline())))
        f.readline()
        if params[0] <= 0 or params[1] <= 0:
            raise ValueError("Time and/or step size parameter non-positive")
        while True:
            name = f.readline()[:-1]
            if name == "":
                break
            colour,marker = f.readline()[:-1].split(",")
            m = float(f.readline())
            r_x = float(f.readline())
            r_y = int(float(f.readline()))
            v_x = float(f.readline())*(60*60*24)
            v_y = float(f.readline())*(60*60*24)
            f.readline()
            params.append(body.Body(name, colour, marker, m, r_x, r_y, v_x, v_y))
    return params
