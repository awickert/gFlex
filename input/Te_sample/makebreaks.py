#! /usr/bin/python

import numpy as np
from matplotlib import pyplot as plt

TeScalar = 80000
# Te = TeScalar*np.ones((100,100))
Te = TeScalar * np.ones(2000)


# Make a discontinuous break: will be unstable, but a check of how to program
# this
def discontinuous(orientation, number, proportion):
    output = np.zeros(Te.shape)
    if orientation == "row":
        output[number, :] = TeScalar
    elif orientation == "column":
        output[:, number] = TeScalar
    return output


def slice2d(rowcol, proportion):
    output = np.zeros(Te.shape)
    for i in rowcol:
        # "max" b/c I am thinking of whole-grid Gaussian in future
        # max(output[i-1:i+1,:])
        output[i - 1 : i + 1, :] = proportion * TeScalar
        output[:, i - 1 : i + 1] = proportion * TeScalar
    return output


# output = np.zeros(Te.shape)
# rowcol = [25,50,75]
# for i in rowcol:
#  output[

# Te -= slice2d([25,50,75],.5)

# for i in 25,50,75:
#  Te-=discontinuous('row',i,.99)
#  Te-=discontinuous('column',i,.99)


# Make a Gaussian function in the middle of a grid with the shape of mine
def gaussian(rowcol, proportion):
    # Only for square grids
    g = np.zeros(Te.shape)
    for i in rowcol:
        a = TeScalar * proportion
        b = i
        c = 8.0
        x = np.arange(0, len(Te))
        gaussian1d = a * np.exp((-((x - b) ** 2)) / (2 * c**2))
        for i in range(len(Te)):
            if len(Te.shape) == 2:
                for j in range(len(Te)):
                    g[i, j] = max(g[i, j], gaussian1d[j])
                    g[j, i] = max(g[j, i], gaussian1d[j])
            elif len(Te.shape) == 1:
                g[i] = max(g[i], gaussian1d[i])
            else:
                print("Error!")
                raise SystemExit
    return g


# Te -= gaussian([25,75],.8)
Te -= gaussian([200, 1000, 1800], 0.99)

if len(Te.shape) == 2:
    plt.imshow(Te)
    plt.colorbar()
elif len(Te.shape) == 1:
    plt.plot(Te)

plt.show()
