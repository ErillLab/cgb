import math


def mean(xs):
    """Returns the average of a given list of numbers."""
    return sum(xs) / float(len(xs))


def log2(x):
    """Returns the logarithm in base 2."""
    return math.log(x, 2)
