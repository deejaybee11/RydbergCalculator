#include <iostream>
#include <fstream>

void save_file(double *rdata, double *wfdata, int length, char* filename) {

    remove(filename);
    std::ofstream output(filename);
    output << length << "\n";
    for (int i = 0; i < length; i++) {
        output << rdata[i] << ";" << wfdata[i] << "\n";
    }
    output.close();
}
