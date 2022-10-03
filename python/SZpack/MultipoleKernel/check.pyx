cdef extern from "pyFunctionWrapper.h":
    cdef cppclass PyElectronDist:
        PyElectronDist()
        PyElectronDist(object)

cdef extern from "test.h":
    double call_some_std_func(PyElectronDist) except +


def example(a):
    print(a)
    return a

def test_example():
    cdef PyElectronDist f = PyElectronDist(example)
    return call_some_std_func(f)

def test_call(func):
    cdef PyElectronDist f = PyElectronDist(func)
    return call_some_std_func(f)

#NB: If the function isn't compiled in cython an UnboundLocalError is thrown