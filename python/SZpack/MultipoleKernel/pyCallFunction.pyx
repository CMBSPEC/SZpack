cdef public double callDistribution(func, double p0):
    return func(p0)

cdef public double callTest(func):
    return func(2.0)