#include "box.h"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>
namespace py = pybind11;
PYBIND11_MODULE(ctrm, m) {
	
	py::class_<box> c(m, "box");
	c.def(py::init<int,int,int,int,int,int,int,int, int, int>(),py::arg("xmax"),py::arg("ymax"),  py::arg("startRad"), py::arg("endRad"), py::arg("dx"), py::arg("dy"), py::arg("nsamp"), py::arg("coordcy0"), py::arg("vrange"), py::arg("dv"))
		.def("readCoord", &box::readCoord)
		.def("readVelo", &box::readVelo)
		.def("setEnergy", &box::setEnergy)
		.def("setCoord", &box::setCoord)
		.def("readEnergy", &box::readEnergy)
		.def("createImageSpace", &box::createImageSpace)
		.def("CalcSurfaceDist", &box::CalcSurfaceDist)
		.def("corrolationOnGeo", &box::corrolationOnGeo)
		.def("getIP", &box::getIP)
		.def("writeIP", &box::writeIP)
		.def("writeSemblence", &box::writeSemblence)
		.def("getCoord", &box::getCoord)
		.def("getSample", &box::getSample)
		.def("getEnergy", &box::getEnergy)
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
	/*box mainBox(1,2,3,4,1,1,1);
	mainBox.readCoord("asd");
	mainBox.readVelo("asd");
	mainBox.readEnergy("");
	mainBox.createImageSpace();
	mainBox.CalcSurfaceDist();
	mainBox.writeIP();
	mainBox.corrolationOnGeo();
	mainBox.writeSemblence();*/

	cin;
 	return 0;
}

