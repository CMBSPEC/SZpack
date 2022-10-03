#include <functional>

inline double call_some_std_func(std::function<double(double)> func) {
    return func(2.0);
}