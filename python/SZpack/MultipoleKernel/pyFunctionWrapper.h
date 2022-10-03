#include <Python.h>
#include <string>
#include "pyCallFunction.h"

// Solution courtesy of https://stackoverflow.com/a/39052204
class PyElectronDist {
public:
    PyObject* function;
public:
    // constructors and destructors mostly do reference counting
    PyElectronDist(PyObject* o): function(o) {
        Py_XINCREF(o);
    }

    ~PyElectronDist() {
        Py_XDECREF(function);
    }

    PyElectronDist(const PyElectronDist& wrapper): PyElectronDist(wrapper.function) {}

    PyElectronDist(PyElectronDist&& wrapper): function(wrapper.function) {
        wrapper.function = 0;
    }

    PyElectronDist(): PyElectronDist(nullptr) {}


   PyElectronDist& operator=(const PyElectronDist& wrapper) {
        PyElectronDist temp = wrapper;
        return (*this = std::move(temp));
    }

    PyElectronDist& operator=(PyElectronDist&& wrapper) {
        function = wrapper.function;
        wrapper.function = 0;
        return *this;
    }

    double operator()(double eta) {
        if (function) { // nullptr check 
            return callDistribution(function,eta); // This will be checked in python
        }
        return 0.0;
    }
};
