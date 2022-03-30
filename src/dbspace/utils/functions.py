def unity(x):
    return x

def quadratic(x, c=[0, 0]):
    return (x - c[0]) * (x - c[1])

# This function is a general function that finds the "nearest" object to the pivot
# Used in this module for finding the nearest datetime
def nearest(items, pivot):
    return min(items, key=lambda x: abs(x - pivot))
