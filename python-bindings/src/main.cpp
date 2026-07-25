#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

#include "../include/Python_ConstantQTransform.h"
#include "../include/Python_ResamplingFilterbank.h"
#include "../include/Python_SlidingCqt.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

static constexpr bool use_windowing{true};

PYBIND11_MODULE(prtcqt, m)
{
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: prtcqt

        .. autosummary::
           :toctree: _generate

    )pbdoc";

    py::class_<Cqt::ScheduleElement>(m, "ScheduleElement")
        .def(py::init<const int, const int, const int>())
        .def("sample", &Cqt::ScheduleElement::sample)
        .def("octave", &Cqt::ScheduleElement::octave)
        .def("delayOctaveRate", &Cqt::ScheduleElement::delayOctaveRate);

    py::class_<Cqt::Python_ConstantQTransform<12, 9>>(m, "Cqt12")
        .def(py::init<>())
        .def("init", &Cqt::Python_ConstantQTransform<12, 9>::init)
        .def("initFs", &Cqt::Python_ConstantQTransform<12, 9>::initFs)
        .def("inputBlock", &Cqt::Python_ConstantQTransform<12, 9>::inputBlock)
        .def("outputBlock", &Cqt::Python_ConstantQTransform<12, 9>::outputBlock)
        .def("getCqtSchedule", &Cqt::Python_ConstantQTransform<12, 9>::getCqtSchedule)
        .def("cqt", &Cqt::Python_ConstantQTransform<12, 9>::cqt)
        .def("icqt", &Cqt::Python_ConstantQTransform<12, 9>::icqt)
        .def("getOctaveCqtBuffer", &Cqt::Python_ConstantQTransform<12, 9>::getOctaveCqtBuffer);

    py::class_<Cqt::Python_ConstantQTransform<24, 9>>(m, "Cqt24")
        .def(py::init<>())
        .def("init", &Cqt::Python_ConstantQTransform<24, 9>::init)
        .def("initFs", &Cqt::Python_ConstantQTransform<24, 9>::initFs)
        .def("inputBlock", &Cqt::Python_ConstantQTransform<24, 9>::inputBlock)
        .def("outputBlock", &Cqt::Python_ConstantQTransform<24, 9>::outputBlock)
        .def("getCqtSchedule", &Cqt::Python_ConstantQTransform<24, 9>::getCqtSchedule)
        .def("cqt", &Cqt::Python_ConstantQTransform<24, 9>::cqt)
        .def("icqt", &Cqt::Python_ConstantQTransform<24, 9>::icqt)
        .def("getOctaveCqtBuffer", &Cqt::Python_ConstantQTransform<24, 9>::getOctaveCqtBuffer);

    py::class_<Cqt::Python_SlidingCqt<24, 9, use_windowing>>(m, "SlidingCqt24")
        .def(py::init<>())
        .def("init", &Cqt::Python_SlidingCqt<24, 9, use_windowing>::init)
        .def("inputBlock", &Cqt::Python_SlidingCqt<24, 9, use_windowing>::inputBlock)
        .def("outputBlock", &Cqt::Python_SlidingCqt<24, 9, use_windowing>::outputBlock)
        .def("getOctaveValues", &Cqt::Python_SlidingCqt<24, 9, use_windowing>::getOctaveValues)
        .def("getOctaveBinFreqs", &Cqt::Python_SlidingCqt<24, 9, use_windowing>::getOctaveBinFreqs);

    py::class_<Cqt::Python_SlidingCqt<12, 9, use_windowing>>(m, "SlidingCqt12")
        .def(py::init<>())
        .def("init", &Cqt::Python_SlidingCqt<12, 9, use_windowing>::init)
        .def("inputBlock", &Cqt::Python_SlidingCqt<12, 9, use_windowing>::inputBlock)
        .def("outputBlock", &Cqt::Python_SlidingCqt<12, 9, use_windowing>::outputBlock)
        .def("getOctaveValues", &Cqt::Python_SlidingCqt<12, 9, use_windowing>::getOctaveValues)
        .def("getOctaveBinFreqs", &Cqt::Python_SlidingCqt<12, 9, use_windowing>::getOctaveBinFreqs);

    py::class_<Cqt::Python_ResamplingFilterbank<9>>(m, "ResamplingFilterbank9")
        .def(py::init<>())
        .def("init", &Cqt::Python_ResamplingFilterbank<9>::init)
        .def("process", &Cqt::Python_ResamplingFilterbank<9>::process)
        .def("getProcessingBlockSize", &Cqt::Python_ResamplingFilterbank<9>::getProcessingBlockSize)
        .def("getLatencySamples", &Cqt::Python_ResamplingFilterbank<9>::getLatencySamples);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
