"""
Provides a class to store data about the state of a body

George Fryer 2018
"""


import numpy



class Body:
    """
    Stores data about an individual body
    """


    def __init__(self, name, colour, marker, m, r_x, r_y, v_x, v_y):
        self._name = name
        self._colour = colour
        self._marker = marker
        self.m = m
        self._r = numpy.array([r_x, r_y])
        self._v = numpy.array([v_x, v_y])


    @property
    def name(self):
        return self._name

    @property
    def colour(self):
        return self._colour

    @property
    def marker(self):
        return self._marker

    @property
    def m(self):
        return self._m

    @m.setter
    def m(self, m):
        if m <= 0:
            raise ValueError("Non-positive mass for body " + self.name)
        self._m = m

    @property
    def r(self):
        return self._r

    @property
    def v(self):
        return self._v
