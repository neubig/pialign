#include "pialign.h"

using namespace pialign;

int main(int argc, const char** argv) {
    PIAlign model;
    model.loadConfig(argc,argv);
    model.loadCorpora();
    model.initialize();
    model.train();
}
