#ifndef _TEST_HPP_
#define _TEST_HPP_

#include "../domain/domain.hpp"
#include "../carpio_define.hpp"
#include "../domain/io_grid.hpp"


namespace carpio {
    // neighbor find
    void test_neighbor_find() {
        Grid<Float, Float, 3> g(3, 0, 1, //
                                3, 0, 1, //
                                3, 0, 1);
        g.connect_root();

        std::list<vtkSmartPointer<vtkProp> > actors;
        actors.push_back(vtk_new_actor_grid_root(g));
        actors.push_back(vtk_new_actor_axes(0, 0, 0));
        Node<Float, Float, 3> *cn = g(1, 1, 1);
        actors.push_back(vtk_new_actor(cn));

        for (int i = 0; i < 8; i++) {

            Node<Float, Float, 3> *pn = cn->get_root_neighbor(56 + i);
            vtkSmartPointer<vtkActor> actor = vtk_new_actor(pn);
            actor->GetProperty()->SetColor(1 / 8.0 * i, 0.0, 0.0); //(R,G,B)
            actors.push_back(actor);
        }
        vtk_show_actor(actors);

    }

    void test_adaptive_init() {
        Grid<Float, Float, 3> g(1, 0, 1, //
                                1, 0, 1, //
                                1, 0, 1);
        g.connect_root();

        Adaptive<Float, Float, 3> adp(&g, 2, 5);

        adp.adapt();


        std::list<vtkSmartPointer<vtkProp> > actors;
        actors.push_back(vtk_new_actor_grid_leaf(g));
        actors.push_back(vtk_new_actor_axes(0, 0, 0));
        vtk_show_actor(actors);

    }


    void test_adaptive_find() {
        Grid<Float, Float, 3> g(2, 0, 1, //
                                2, 0, 1, //
                                2, 0, 1);
        delete g.nodes(1,1,1);
        g.nodes(1,1,1)= nullptr;
        g.connect_root();

        Adaptive<Float, Float, 3> adp(&g, 2, 5);

        adp.adapt();

        Node<Float, Float, 3> *cn = GetpNodeAt(&g, 0.9, 0.9, 0.9);


        std::list<vtkSmartPointer<vtkProp> > actors;
        actors.push_back(vtk_new_actor_grid_leaf(g));
        actors.push_back(vtk_new_actor(cn));
        st n = 6;
        std::cout<<"ori ==== "<<cn->get_idx()<<endl;
        for (size_t i: IdxRange<size_t>(0, 26)) {
            std::cout<<"-----------"<<endl;
            std::cout<<"idx = "<< i<<endl;
            std::cout<<"dir = "<<DirectionInOrder(i)<<endl;
            Node<Float, Float, 3> *nei = cn->get_neighbor(cn, DirectionInOrder(i));
            if (nei == nullptr) {
                std::cout << "empty nei  " << i << endl;
            } else {
                vtkSmartPointer<vtkActor> actornei = vtk_new_actor(nei);
                actornei->GetProperty()->SetColor(1.0 / 8 * i, 0.0, 0.0); //(R,G,B)
                actors.push_back(actornei);
            }
        }
        //Node<Float, Float, 3> *nei = cn->get_neighbor(cn, DirectionInOrder(22));
        //vtkSmartPointer<vtkActor> actornei = vtk_new_actor(nei);
        //actornei->GetProperty()->SetColor(1.0 / 4, 0.0, 0.0); //(R,G,B)
        //actors.push_back(actornei);
        actors.push_back(vtk_new_actor_axes(0, 0, 0));
        vtk_show_actor(actors);

    }

}


#endif
