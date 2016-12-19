#ifndef _VTK_HPP_
#define _VTK_HPP_


#include <memory>
#include <string>

namespace carpio{

namespace VTU{

void unstructured_grid_file_head(std::ofstream& fs);
void unstructured_grid_file_end(std::ofstream& fs);

//void drawtofile_vtu(std::ofstream& fs, ListT<Point3D>&);



}

}

#endif
