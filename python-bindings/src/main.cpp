#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

#include "../include/Python_ConstantQTransform.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(Cqt, m) 
{
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: Cqt

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    PYBIND11_NUMPY_DTYPE(Cqt::ScheduleElement, sample, octave, delayOctaveRate);

    py::class_<Cqt::Python_ConstantQTransform<12, 9>>(m, "ConstantQTransform")
        .def(py::init<>())
        .def("init", py::overload_cast<int>(&Cqt::Python_ConstantQTransform<12, 9>::init))
        .def("initFs", &Cqt::Python_ConstantQTransform<12, 9>::initFs)
        .def("inputBlock", &Cqt::Python_ConstantQTransform<12, 9>::Python_inputBlock)
        .def("getCqtSchedule", &Cqt::Python_ConstantQTransform<12, 9>::Python_getCqtSchedule)
        .def("cqt", &Cqt::Python_ConstantQTransform<12, 9>::Python_cqt)
        .def("getOctaveCqtBuffer", &Cqt::Python_ConstantQTransform<12, 9>::getOctaveCqtBuffer);


#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}









