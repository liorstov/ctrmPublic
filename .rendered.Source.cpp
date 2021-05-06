// cppimport
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "box.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;


int main() {
	box mainBox;
	mainBox.readCoord("asd");
	mainBox.readVelo("asd");
	mainBox.readEnergy("");
	mainBox.createImageSpace();
	mainBox.CalcSurfaceDist();
	mainBox.writeIP();
	mainBox.corrolationOnGeo();
	mainBox.writeSemblence();

	cin;
 	return 0;
}

PYBIND11_MODULE(crp, m) {

	py::class_<box> c(m, "box");
		c.def(py::init<>())
			.def("readCoord", &box::readCoord)
			.def("readVelo", &box::readVelo)
			.def("readEnergy", &box::readEnergy)
			.def("createImageSpace", &box::createImageSpace)
			.def("CalcSurfaceDist", &box::CalcSurfaceDist)
			.def("corrolationOnGeo", &box::corrolationOnGeo)
			.def("getIP", &box::getIP)
			.def("getCoord", &box::getCoord)
	;

#ifdef VERSION_INFO
	m.attr("__version__") = VERSION_INFO;
#else
	m.attr("__version__") = "dev";
#endif
}
/*

*/