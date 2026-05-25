#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

int main() {
  //std::ifstream fsin("KURAMA_field.map");
  std::ifstream fsin("KURAMA/KuramaFieldMap_E07_20160607_cordMod");

  if(!fsin){
    std::cerr  << ": file open fail" << std::endl;
    exit(-1);
  }

  int Nx, Ny, Nz;
  double X0, Y0, Z0, dX, dY, dZ;

  if( !(fsin >> Nx >> Ny >> Nz >> X0 >> Y0 >> Z0 >> dX >> dY >> dZ ) ){
    std::cerr << ": Invalid format " << std::endl;
    exit(-1);
  }

  double x, y, z, bx, by, bz;
  int npoint=0;

  double pre_x=0, pre_y=0, pre_z=0, pre_bx=0, pre_by=0, pre_bz=0;
  int start = 0;

  while(fsin){
    fsin >> x >> y >> z >> bx >> by >> bz;
    /*
    if (z < pre_z && start != 0) {
      for (double zz=pre_z; zz<=40.; zz +=1.0) {
	//printf("%.6lf   %.6lf   %.6lf   %.6lf   %.6lf   %.6lf \n", pre_x, pre_y, zz, pre_bx, pre_by, pre_bz);
	//std::cout << pre_x << "   " << pre_y << "   " << zz << "   " << pre_bx << "   " << pre_by << "   " << pre_bz << std::endl;	
      }

      for (double zz=-40.; zz<z; zz += 1.0) {
	//printf("%.6lf   %.6lf   %.6lf   %.6lf   %.6lf   %.6lf \n", x, y, zz, bx, by, bz);
	//std::cout << x << "   " << y << "   " << zz << "   " << bx << "   " << by << "   " << bz << std::endl;	
      }
      //printf("%.6lf   %.6lf   %.6lf   %.6lf   %.6lf   %.6lf \n", x, y, z, bx, by, bz);
      //std::cout << x << "   " << y << "   " << z << "   " << bx << "   " << by << "   " << bz << std::endl;	
    } else {
      //printf("%.6lf   %.6lf   %.6lf   %.6lf   %.6lf   %.6lf \n", x, y, z, bx, by, bz);
      //std::cout << x << "   " << y << "   " << z << "   " << bx << "   " << by << "   " << bz << std::endl;	
    }
    */
    pre_x = x;
    pre_y = y;
    pre_z = z;
    pre_bx = bx;
    pre_by = by;
    pre_bz = bz;
 
    if (std::abs(y-(0))<0.1 && std::abs(z)<0.1) {
      printf("%.6lf   %.6lf   %.6lf   %.6lf   %.6lf   %.6lf \n", x, y, z, bx, by, bz);
    }
    
    start = 1;
  }


  return 0;
}
