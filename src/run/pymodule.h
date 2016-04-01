#ifndef _PYMODULE_HPP_
#define _PYMODULE_HPP_

#include "../carpio_define.hpp"
#include <Python.h>
namespace carpio {
class PyParse {
protected:
	std::string filename;
	PyObject* pModule;
public:
	PyParse(const std::string& fn) :
			filename(fn) {
		Py_Initialize();        //Initialize
		PyRun_SimpleString("import sys");
		PyRun_SimpleString("sys.path.append('./')");
		pModule = PyImport_ImportModule("para");
		if (!pModule) {
			std::cerr << "get module failed!" << "\n";
			ASSERT(false);
		}
	}
	~PyParse() {
		Py_DECREF(pModule);
		Py_Finalize();
	}
	//
	PyObject* get_function_op(const std::string& name) const {
		PyObject * po = PyObject_GetAttrString(pModule, name.c_str());
		if (!(po && PyCallable_Check(po))) {
			std::cerr << " >! Get function " << name << " is failed!" << "\n";
			ASSERT(false);
		}
		return po;
	}

	PyObject* get_int_op(const std::string& name) const {
		PyObject * po = PyObject_GetAttrString(pModule, name.c_str());
		if (!(po && PyInt_Check(po))) {
			std::cerr << " >! Get function " << name << " is failed!" << endl;
			ASSERT(false);
		}
		return po;
	}
	PyObject* get_op(const std::string& name) const {
		PyObject * po = PyObject_GetAttrString(pModule, name.c_str());
		if (po==nullptr) {
			std::cerr << " >! Get " << name << " is failed!" << endl;
			ASSERT(false);
		}
		return po;
	}
	int get_int(const std::string& name) const {
		PyObject* po = get_int_op(name);
		return int(PyFloat_AsDouble(po));
	}
	int po_to_int(PyObject* t, const string& err_msg) const {
		if (!PyInt_Check(t)) {
			std::cerr << err_msg << " type=" << t->ob_type->tp_name
					<< " should be int" << "\n";
			exit(0);
		}
		return int(PyFloat_AsDouble(t));
	}

	Float po_to_float(PyObject* t, const string& err_msg) const {
		if (!PyFloat_Check(t)) {
			std::cerr << err_msg << " type=" << t->ob_type->tp_name
					<< " should be float" << "\n";
			exit(0);
		}
		return Float(PyFloat_AsDouble(t));
	}
	std::string po_to_string(PyObject* t, const string& err_msg) const {
		if (1 != PyString_Check(t)) {
			cerr << err_msg << endl;
			exit(0);
		}
		return PyString_AsString(t);
	}
}
;
}

#endif
