import math
import numpy as np


def mean(xs):
    """Returns the average of a given list of numbers."""
    return sum(xs) / float(len(xs))


def log2(x):
    """Returns the logarithm in base 2."""
    return math.log(x, 2)


def weighted_choice(xs, weights, count=1):
    """Given a list of elements, randomly selects one based on weights.
    Args:
        xs (list): elements
        weights (list): weights
    """
    ps = [float(w)/sum(weights) for w in weights]
    return np.random.choice(xs, size=count, p=ps)


def unique(xs, f):
    """Makes the list xs unique by comparing elements with f(x)."""
    d = {f(x): x for x in xs}
    return d.values()
