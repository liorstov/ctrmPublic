#include "box.h"
#include <filesystem>

int main(int argc, char *argv[])
{

    int *params = new int[argc - int(4)];
    printf("Program Name Is: %s\n", argv[0]);
    std::filesystem::path p{std::string(argv[20])};
    

    for (int i = 1; i < argc - 4; i++)
    {
        params[i - 1] = atoi(argv[i]);
    }
    cout << "Asdasdasd"<<std::string(argv[20]).substr(std::string(argv[20]).find_last_of("/\\")+1) + ".npy" << endl;
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
                atof(argv[16]),
                atof(argv[17]));
    mainBox.readCoordnumpy(argv[18]);
    mainBox.readVelo(argv[19]);
    // mainBox.readEnergy(argv[18]);
    mainBox.readEnergyFromNpy(argv[20]);
    mainBox.createImageSpace();
    mainBox.CalcSurfaceDist();
    mainBox.CalcTimeDeltaOnly();
    //mainBox.writeIP();
    //mainBox.corrolationOnGeo(0);
    mainBox.writeSemblenceNpy(std::string(argv[21])+= p.filename().replace_extension("cube"));
    //mainBox.writeTimeDeltasNpy(std::string(argv[21])+= p.filename().replace_extension("deltas"));

    cin;
    return 0;
}
