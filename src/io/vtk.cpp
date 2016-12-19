#include "vtk.hpp"
#include <iostream>
#include <fstream>
namespace carpio{

namespace VTK{


void vtu_unstructured_grid_file_head(std::ofstream& fs) {
	fs << "<?xml version=\"1.0\"?>";
	fs
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">";
	fs << "<UnstructuredGrid>";
}

void vtu_unstructured_grid_file_end(std::ofstream& fs) {
	fs << "</UnstructuredGrid>";
	fs << "</VTKFile>";
}

}
}
