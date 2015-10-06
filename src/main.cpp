#include <iostream>
#include "algebra/polynomial.hpp"
#include "domain/node.hpp"
#include "domain/io_grid.hpp"

#include "_test/_test_grid.hpp"


using namespace std;
using namespace carpio;

void show(const Polynomial<int, Float, int>& poly){
    cout<<" coe   term   exp\n";
    Polynomial<int, double, int>::const_iterator iter;
    for(iter = poly.begin(); iter!=poly.end();++iter){
        cout<< iter->second <<"      ";
        cout<< iter->first.first <<"      ";
        cout<< iter->first.second <<"      \n";
    }
}

inline void test_show_node() {
    typedef Node<CooValueType, ValueType, 3> Node;

    Node* pnode = new Node(nullptr, 0, 0, 0, 0, //
                           0.5, 0.5, //
                           0.5, 0.5, //
                           0.5, 0.5);
    //

    vtk_show(pnode);

    delete pnode;

}

int test_ploynomial() {
    Polynomial<int, Float, int> poly;
    poly.insert(typename Polynomial<int, Float, int>::Term(1, 1, 1));
    poly.insert(typename Polynomial<int, Float, int>::Term(1, 2, 1));
    cout << "size     " << poly.size() << "\n";
    show(poly);
    cout << " ============ \n";
    Polynomial<int, Float, int> poly2;
    poly2.insert(typename Polynomial<int, Float, int>::Term(1, 2, 1));
    poly2.insert(typename Polynomial<int, Float, int>::Term(1, 4, 0));
    poly2.insert(typename Polynomial<int, Float, int>::Term(1, 5, 0));
    show(poly2);
    cout << "size 2   " << poly2.size() << "\n";

    poly.minus(poly2);
    poly.concise();

    show(poly);
    test_show_node();
    cout << "Hello, World!" << endl;
    return 0;
}

int main(){
    test_neighbor_find();
    cout<< "======Fin======" << endl;
    return 0;
}