/// @file python/bind_frame.cpp
/// @brief CpuFrame binding: numpy <-> cpu::Frame transfer.

#include "cpu_frame_wrapper.h"

void bind_frame(nb::module_& m) {
    nb::class_<CpuFrame>(
        m, "CpuFrame",
        "Internal CPU frame wrapper (numpy <-> cpu::Frame bridge).")
        .def(nb::init<int>(), nb::arg("n_atoms"),
             "Create a CpuFrame with the given number of atoms.")
        .def_prop_ro("n_atoms", &CpuFrame::n_atoms)
        // Populate
        .def("set_positions", &CpuFrame::set_positions, nb::arg("x"),
             nb::arg("y"), nb::arg("z"))
        .def("set_velocities", &CpuFrame::set_velocities, nb::arg("vx"),
             nb::arg("vy"), nb::arg("vz"))
        .def("set_mass", &CpuFrame::set_mass, nb::arg("mass"))
        .def("set_type", &CpuFrame::set_type, nb::arg("type"))
        .def("set_box", &CpuFrame::set_box, nb::arg("box"))
        .def("set_bonds", &CpuFrame::set_bonds, nb::arg("atomi"),
             nb::arg("atomj"), nb::arg("type"))
        .def("set_angles", &CpuFrame::set_angles, nb::arg("atomi"),
             nb::arg("atomj"), nb::arg("atomk"), nb::arg("type"))
        .def("set_dihedrals", &CpuFrame::set_dihedrals, nb::arg("atomi"),
             nb::arg("atomj"), nb::arg("atomk"), nb::arg("atoml"),
             nb::arg("type"))
        // Extract
        .def("get_x", &CpuFrame::get_x)
        .def("get_y", &CpuFrame::get_y)
        .def("get_z", &CpuFrame::get_z)
        .def("get_vx", &CpuFrame::get_vx)
        .def("get_vy", &CpuFrame::get_vy)
        .def("get_vz", &CpuFrame::get_vz)
        .def("get_fx", &CpuFrame::get_fx)
        .def("get_fy", &CpuFrame::get_fy)
        .def("get_fz", &CpuFrame::get_fz)
        .def("get_mass", &CpuFrame::get_mass)
        .def("get_type", &CpuFrame::get_type)
        .def("get_box", &CpuFrame::get_box)
        .def_prop_ro("energy", &CpuFrame::get_energy)
        .def_prop_ro("kinetic", &CpuFrame::get_kinetic);
}
