#include "box.h"

//
//
//#define STRINGIFY(x) #x
//#define MACRO_STRINGIFY(x) STRINGIFY(x)
//
//namespace py = pybind11;
//PYBIND11_MODULE(ctrm, m) {
//
//	/*py::class_<box> c(m, "box");
//	c.def(py::init<int, int,int, int, int, int, int,int, int, int, int, int, int, Eigen::MatrixXf, Eigen::MatrixXf, Eigen::MatrixXf, int>(),py::arg("xmax"),py::arg("ymax"), py::arg("xmin"), py::arg("ymin"),  py::arg("startRad"), py::arg("endRad"), py::arg("dx"), py::arg("dy"), py::arg("jbeg"), py::arg("jend"), py::arg("vrange"), py::arg("dv"), py::arg("minDist"), py::arg("veloList"), py::arg("geom"), py::arg("energy"), py::arg("windowSize"))
//		.def("readCoord", &box::readCoord)
//		.def("setEnergy", &box::setEnergy)
//		.def("setCoord", &box::setCoord)
//		.def("readEnergy", &box::readEnergy)
//		.def("createImageSpace", &box::createImageSpace)
//		.def("CalcSurfaceDist", &box::CalcSurfaceDist)
//		.def("corrolationOnGeo", &box::corrolationOnGeo,py::call_guard<py::scoped_ostream_redirect,	py::scoped_estream_redirect>())
//		.def("getIP", &box::getIP)
//		.def("writeIP", &box::writeIP)
//		.def("writeSemblence", &box::writeSemblence)
//		.def("getCoord", &box::getCoord)
//		.def("getSample", &box::getSample)
//		.def("getEnergy", &box::getEnergy)
//		;*/
//
//#ifdef VERSION_INFO
//    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
//#else
//    m.attr("__version__") = "dev";
//#endif
//}

int main(int argc, char *argv[])
{

    int *params = new int[argc - int(4)];
    printf("Program Name Is: %s", argv[0]);
    for (int i = 1; i < argc - 4; i++)
    {
        params[i - 1] = atoi(argv[i]);
    }

    //box mainBox('0', '0', '30', '-40', '0', '35', '1', '1', '300', '350', '100', '20', '9999','1');
    box mainBox(atoi(argv[1]),
                atoi(argv[2]),
                atoi(argv[3]),
                atoi(argv[4]),
                atoi(argv[5]),
                atoi(argv[6]),
                atoi(argv[7]),
                atoi(argv[8]),
                atoi(argv[9]),
                atoi(argv[10]),
                atoi(argv[11]),
                atoi(argv[12]),
                atoi(argv[13]),
                atoi(argv[14]),
                atoi(argv[15]),
                atof(argv[16]));
    mainBox.readCoord(argv[17]);
    mainBox.readVelo(argv[18]);
    // mainBox.readEnergy(argv[18]);
    mainBox.readEnergyFromNpy(argv[19]);
    mainBox.createImageSpace();
    mainBox.CalcSurfaceDist();
    mainBox.CalcTimeDeltaOnly();
    //mainBox.writeIP();
    //mainBox.corrolationOnGeo(0);
    mainBox.writeSemblenceNpy(std::string(argv[19]) + std::string(argv[20]) + ".npy");
    mainBox.writeTimeDeltasNpy(std::string(argv[19]) + std::string(argv[20]) + "_deltas" + ".npy");

    cin;
    return 0;
}
