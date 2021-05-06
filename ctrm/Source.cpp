#include "box.h"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>
namespace py = pybind11;


PYBIND11_MODULE(crp, m) {

	py::class_<box> c(m, "box");
	c.def(py::init<int,int,int,int>())
		.def("readCoord", &box::readCoord)
		.def("readVelo", &box::readVelo)
		.def("readEnergy", &box::readEnergy)
		.def("createImageSpace", &box::createImageSpace)
		.def("CalcSurfaceDist", &box::CalcSurfaceDist)
		.def("corrolationOnGeo", &box::corrolationOnGeo)
		.def("getIP", &box::getIP)
		.def("writeIP", &box::writeIP)
		.def("writeSemblence", &box::writeSemblence)
		.def("getCoord", &box::getCoord)
		.def("getSample", &box::getSample)
		;

#ifdef VERSION_INFO
	m.attr("__version__") = VERSION_INFO;
#else
	m.attr("__version__") = "dev";
#endif
}


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)




int main() {
	box mainBox(1,2,3,4);
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

