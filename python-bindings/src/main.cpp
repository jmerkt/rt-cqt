#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/python_constant_q_transform.h"
#include "../include/python_resampling_filterbank.h"
#include "../include/python_sliding_cqt.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

static constexpr bool USE_WINDOWING{true};

PYBIND11_MODULE(prtcqt, module)
{
    module.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: prtcqt

        .. autosummary::
           :toctree: _generate

    )pbdoc";

    py::class_<rt_cqt::ScheduleElement>(module, "ScheduleElement")
        .def(py::init<const int, const int, const int>())
        .def("sample", &rt_cqt::ScheduleElement::sample)
        .def("octave", &rt_cqt::ScheduleElement::octave)
        .def("delay_at_octave_rate", &rt_cqt::ScheduleElement::delay_at_octave_rate)
        .def("synthesis_offset", &rt_cqt::ScheduleElement::synthesis_offset);

    py::class_<rt_cqt::PythonConstantQTransform<12, 9>>(module, "Cqt12")
        .def(py::init<>())
        .def("init", &rt_cqt::PythonConstantQTransform<12, 9>::init)
        .def("init_sample_rate", &rt_cqt::PythonConstantQTransform<12, 9>::init_sample_rate)
        .def("input_block", &rt_cqt::PythonConstantQTransform<12, 9>::input_block)
        .def("output_block", &rt_cqt::PythonConstantQTransform<12, 9>::output_block)
        .def("get_cqt_schedule", &rt_cqt::PythonConstantQTransform<12, 9>::get_cqt_schedule)
        .def("cqt", &rt_cqt::PythonConstantQTransform<12, 9>::cqt)
        .def("icqt", &rt_cqt::PythonConstantQTransform<12, 9>::icqt)
        .def("get_octave_cqt_buffer", &rt_cqt::PythonConstantQTransform<12, 9>::get_octave_cqt_buffer);

    py::class_<rt_cqt::PythonConstantQTransform<24, 9>>(module, "Cqt24")
        .def(py::init<>())
        .def("init", &rt_cqt::PythonConstantQTransform<24, 9>::init)
        .def("init_sample_rate", &rt_cqt::PythonConstantQTransform<24, 9>::init_sample_rate)
        .def("input_block", &rt_cqt::PythonConstantQTransform<24, 9>::input_block)
        .def("output_block", &rt_cqt::PythonConstantQTransform<24, 9>::output_block)
        .def("get_cqt_schedule", &rt_cqt::PythonConstantQTransform<24, 9>::get_cqt_schedule)
        .def("cqt", &rt_cqt::PythonConstantQTransform<24, 9>::cqt)
        .def("icqt", &rt_cqt::PythonConstantQTransform<24, 9>::icqt)
        .def("get_octave_cqt_buffer", &rt_cqt::PythonConstantQTransform<24, 9>::get_octave_cqt_buffer);

    py::class_<rt_cqt::PythonSlidingCqt<24, 9, USE_WINDOWING>>(module, "SlidingCqt24")
        .def(py::init<>())
        .def("init", &rt_cqt::PythonSlidingCqt<24, 9, USE_WINDOWING>::init)
        .def("input_block", &rt_cqt::PythonSlidingCqt<24, 9, USE_WINDOWING>::input_block)
        .def("output_block", &rt_cqt::PythonSlidingCqt<24, 9, USE_WINDOWING>::output_block)
        .def("get_octave_values", &rt_cqt::PythonSlidingCqt<24, 9, USE_WINDOWING>::get_octave_values)
        .def("get_octave_bin_frequencies", &rt_cqt::PythonSlidingCqt<24, 9, USE_WINDOWING>::get_octave_bin_frequencies);

    py::class_<rt_cqt::PythonSlidingCqt<12, 9, USE_WINDOWING>>(module, "SlidingCqt12")
        .def(py::init<>())
        .def("init", &rt_cqt::PythonSlidingCqt<12, 9, USE_WINDOWING>::init)
        .def("input_block", &rt_cqt::PythonSlidingCqt<12, 9, USE_WINDOWING>::input_block)
        .def("output_block", &rt_cqt::PythonSlidingCqt<12, 9, USE_WINDOWING>::output_block)
        .def("get_octave_values", &rt_cqt::PythonSlidingCqt<12, 9, USE_WINDOWING>::get_octave_values)
        .def("get_octave_bin_frequencies", &rt_cqt::PythonSlidingCqt<12, 9, USE_WINDOWING>::get_octave_bin_frequencies);

    py::class_<rt_cqt::PythonResamplingFilterbank<9>>(module, "ResamplingFilterbank9")
        .def(py::init<>())
        .def("init", &rt_cqt::PythonResamplingFilterbank<9>::init)
        .def("process", &rt_cqt::PythonResamplingFilterbank<9>::process)
        .def("get_processing_block_size", &rt_cqt::PythonResamplingFilterbank<9>::get_processing_block_size)
        .def("get_latency_samples", &rt_cqt::PythonResamplingFilterbank<9>::get_latency_samples);

#ifdef VERSION_INFO
    module.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    module.attr("__version__") = "dev";
#endif
}
