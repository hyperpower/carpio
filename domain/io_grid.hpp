#ifndef IOGRID_H_
#define IOGRID_H_


#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "boundary.hpp"
#include "grid.hpp"

// vtk include
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkVoxel.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
//

namespace carpio {


struct vtk_RenderingArg {
	int flag;
};

int vtk_show_actor(vtkSmartPointer<vtkActor> actor[], vtkIdType len) {
	//
	//actor->GetProperty()->SetRepresentationToWireframe();
	// Create a renderer, render window and interactor
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<
			vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
			vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	// Add the actors to the scene
	for(vtkIdType i = 0;i<len;++i){
		renderer->AddActor(actor[i]);
	}
	renderer->SetBackground(.1, .2, .3); // Background color dark blue

	// Render and interact
	renderWindow->Render();
	renderWindowInteractor->Start();

	return EXIT_SUCCESS;
}
static const Orientation vtk_VOXEL[8][3] = { //
		//         x    y    z
				{ _M_, _M_, _M_ }, //
				{ _P_, _M_, _M_ }, //
				{ _M_, _P_, _M_ }, //
				{ _P_, _P_, _M_ }, //
				{ _M_, _M_, _P_ }, //
				{ _P_, _M_, _P_ }, //
				{ _M_, _P_, _P_ }, //
				{ _P_, _P_, _P_ }, //
		};

template<typename COO_VALUE, typename VALUE, int DIM>
void _vtkUnstructuredGrid_add_node(
		const Node<COO_VALUE, VALUE, DIM>* pn, //
		vtkSmartPointer<vtkPoints> points,
		vtkSmartPointer<vtkUnstructuredGrid> ugrid) {
	ASSERT(pn!=nullptr);
	typedef COO_VALUE vt;
	typedef Node<COO_VALUE, VALUE, DIM> node_t;
	vtkIdType n = points->GetNumberOfPoints();

	for (size_t i = 0; i < pn->NumVertexes; ++i) {
		vt x = pn->cell->Get(vtk_VOXEL[i][0], _X_);
		vt y = (node_t::Dim >= 2) ? pn->cell->Get(vtk_VOXEL[i][1], _Y_) : 0.0;
		vt z = (node_t::Dim == 3) ? pn->cell->Get(vtk_VOXEL[i][2], _Z_) : 0.0;
		points->InsertNextPoint(x, y, z);
	}

	vtkSmartPointer<vtkVoxel> voxel = vtkSmartPointer<vtkVoxel>::New();
	// append idx

	for (size_t i = 0; i < pn->NumVertexes; ++i) {
		voxel->GetPointIds()->SetId(i, n + i);
	}
	// Insert to Unstrcture Grid
	ugrid->InsertNextCell(voxel->GetCellType(), voxel->GetPointIds());
}

template<typename COO_VALUE, typename VALUE, int DIM>
vtkSmartPointer<vtkActor> vtk_new_actor(const Node<COO_VALUE, VALUE, DIM>* pn) {
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<
			vtkUnstructuredGrid>::New();
	// Add the geometry and topology to the unstructure grid
	_vtkUnstructuredGrid_add_node(pn, points, ugrid);

	ugrid->SetPoints(points);

	// Create a mapper and actor
	vtkSmartPointer<vtkDataSetMapper> mapper =
			vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInput(ugrid);
#else
	mapper->SetInputData(ugrid);
#endif
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	return actor;

}

template<typename COO_VALUE, typename VALUE, int DIM>
vtkSmartPointer<vtkActor> vtk_new_actor(const Grid<COO_VALUE, VALUE, DIM>& grid) {
	typedef Node<COO_VALUE, VALUE, DIM>* pnode;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<
			vtkUnstructuredGrid>::New();
	//
	for (size_t i = 0; i < grid.nodes.size(); ++i) {
		pnode pn = grid.nodes.at_1d(i);
		if (pn != nullptr) {
			_vtkUnstructuredGrid_add_node(pn, points, ugrid);
		}
	}

	ugrid->SetPoints(points);

	// Create a mapper and actor
	vtkSmartPointer<vtkDataSetMapper> mapper =
			vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInput(ugrid);
#else
	mapper->SetInputData(ugrid);
#endif
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetRepresentationToWireframe();
	return actor;
}
/*
 * high level =============================================
 */
template<typename COO_VALUE, typename VALUE, int DIM>
int vtk_show(const Grid<COO_VALUE, VALUE, DIM>& grid) {
	vtkSmartPointer<vtkActor> actor = vtk_new_actor(grid);
	return vtk_show_actor(&actor, 1);
}

template<typename COO_VALUE, typename VALUE, int DIM>
int vtk_show(const Node<COO_VALUE, VALUE, DIM>* pnode) {
	vtkSmartPointer<vtkActor> actor = vtk_new_actor(pnode);
	return vtk_show_actor(&actor, 1);
}


}

#endif
