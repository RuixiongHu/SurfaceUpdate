#include <apf.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <cassert>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <iostream>
#include "lionPrint.h"
#include <parma.h>
#include "spr.h"
#include <apfConvert.h>
//Modify the CMakeList.txt to ensure correct file name etc
//source setup.sh to do cmake
//build with ./build.sh writeVTK
//run with mpirun -n 4 ./build/writeVTK 40x20x5umSi.dmg restart_0.05_.smb restart_0.05_updated


// main surface update function
int main(int argc, char** argv){
  assert(argc==4);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  const char* outFile = argv[3];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  MS_init();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_mesh();
  gmi_register_sim();
  lion_set_verbosity(1);
  gmi_model* model = gmi_load(modelFile);
  apf::Mesh2* mesh = apf::loadMdsMesh(model,meshFile);

  PCU_Barrier();
  apf::writeASCIIVtkFiles(outFile,mesh);
}
