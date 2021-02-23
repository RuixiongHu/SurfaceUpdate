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
//build with ./build.sh SurfaceUpdate
//run with mpirun -n 4 ./build/Fake_Field 40x20x5umSi.dmg restart_0.05_.smb restart_0.05_updated.smb 4

namespace{


	// define a surface class, congtaining three point (triangle)
	class Surface{
  	public:

	  Surface(double *A, double *B, double *C); // constructor

	  // Accessor
	  double getAx() const;
	  double getAy() const;
	  double getAz() const;
	  double getBx() const;
	  double getBy() const;
	  double getBz() const;
	  double getCx() const;
	  double getCy() const;
	  double getCz() const;

	  double* getAptr() ;
	  double* getBptr() ;
	  double* getCptr() ;

 	 // Print
 	 void print();
	
	private:

	  double A_[3];
	  double B_[3];
	  double C_[3];

	};

	Surface::Surface(double *A, double *B, double *C){
  		for (int i=0; i<3; i++){
	    A_[i]=*(A+i);
	    B_[i]=*(B+i);
	    C_[i]=*(C+i);
	  }
	}

	  double Surface::getAx() const{ return A_[0]; }
	  double Surface::getAy() const{ return A_[1]; }
	  double Surface::getAz() const{ return A_[2]; }
	  double Surface::getBx() const{ return B_[0]; }
	  double Surface::getBy() const{ return B_[1]; }
	  double Surface::getBz() const{ return B_[2]; }
	  double Surface::getCx() const{ return C_[0]; }
	  double Surface::getCy() const{ return C_[1]; }
	  double Surface::getCz() const{ return C_[2]; }

	  double* Surface::getAptr() { return A_; }
	  double* Surface::getBptr() { return B_; }
	  double* Surface::getCptr() { return C_; }

	  void Surface::print(){
	    //std::cout<<"This surface contains ("<<A_[0]<<","<<A_[1]<<","<<A_[2]<<"), ("<<B_[0]<<","<<B_[1]<<","<<B_[2]<<"), ("<<C_[0]<<","<<C_[1]<<","<<C_[2]<<").\n";
	  }	 


	//********************************************************
	// Extra Functions

	// This function switch the order of two nodes
	void exchange(int *first, int *second, int size){   
	  int storage[size];
	  for (int i=0; i<size; i++){
	    storage[i]=*(first+i);
	    *(first+i)=*(second+i);
	    *(second+i)=*(storage+i);
	  }
	}

	  // this fuction sort a [0,1][1,2] into p[0,1,2]
	void sort(int *recorder1, int* recorder2, int*result ){
	  int a,b,c,d;
	  a=*recorder1;
	  b=*(recorder1+1);
	  c=*recorder2;
	  d=*(recorder2+1);
	  //std::cout<<"enter sort loop, abcd is "<<a<<b<<c<<d<<"\n";
	  if( (a==c)|(b==c) ){
	    *result=a;
	    *(result+1)=b;
	    *(result+2)=d;
	  }
	  else if( (a==d)|(b==d) ) {
	    *result=a;
	    *(result+1)=b;
	    *(result+2)=c;
	  }
	  else{
	    std::cerr<<"!!!!Can not sort 3 out of 4 since all four vertices are different\n";
	  }
	  //std::cout<<"after sorting, value is "<< *result <<", " << *(result+1) << ", " <<*(result+2)<<"\n";
	}

}  



//*****************************************************************************************************
// main surface update function
int main(int argc, char** argv){
  assert(argc==5);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  const char* outFile = argv[3];
  int partitionFactor = atoi(argv[4]);
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
  apf::Field* Psi_Field = createIPField(mesh,"Psi_1",1,1);
  apf::MeshIterator* it ;   // loop over mesh element 
  apf::MeshEntity* e;  
  std::vector<apf::MeshEntity*> nodes;      // store 2 nodes of an edge
  double qp[3];                           // coordinate of the quadrature point
  double weight[4]={0.25,0.25,0.25,0.25}; // the weighting function for an order 1 tetrahedral 
  apf::Vector3 QPcoordinate;
  double UpdateValue = 1;// fake psi value
  double zero = 0;
    apf::Vector3 v1,v2,v3,v4;           // 2 Vector3 object to store nodes' coordinate for an edge, then 4 of them store element for qp location
      apf::Downward adjacent;         // target container for downward adjacent
      int numAdj;       


    it = mesh->begin(3);
      while( (e = mesh->iterate(it) ) ){
      //std::cout<<"Enter a new mesh for integrrogation\n";
        nodes.clear();
        int NumProjection = 0;
          numAdj = mesh->getDownward(e,0,adjacent);             // get 4 nodes of current mesh element to determine quadrature point
          for (int j = 0; j < numAdj; ++j) {
            nodes.push_back(adjacent[j]);
          }
          mesh -> getPoint_( nodes[0], 0, v1);
          mesh -> getPoint_( nodes[1], 0, v2);
          mesh -> getPoint_( nodes[2], 0, v3);
          mesh -> getPoint_( nodes[3], 0, v4);

        for (int i=0; i<3; i++){                                    // calculate x,y,z of the quadrature point
          qp[i]= weight[0]*v1[i]+weight[1]*v2[i]+weight[2]*v3[i]+weight[3]*v4[i];
          //std::cout<<"Enter new mesh, equation is "<< weight[0]<< "x" <<v1[i]<< "+"<< weight[1]<< "x" <<v2[i]<< "+" <<weight[2]<< "x" <<v3[i]<< "+" <<weight[3]<< "x" <<v4[i]<<"\n";
        }        
        if (qp[2]<20){
            //std::cout<<"found QP, with z=" << qp[2]<<std::endl;
            apf::setComponents( Psi_Field, e, 0, &UpdateValue);
          }
        else{
            apf::setComponents( Psi_Field, e, 0, &zero);          
        }
      }
      mesh->end(it);

		apf::writeASCIIVtkFiles("outVtk",mesh);
		mesh->writeNative(outFile);
	       
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  MS_exit();
  PCU_Comm_Free();
  MPI_Finalize();

	}




