#ifndef _EVENT_H_
#define _EVENT_H_

#include "../carpio_define.hpp"
#include "../utility/format.h"
#include "gnuplot.h"
#include "../calculation/event.hpp"

namespace carpio {

template<typename COO_VALUE, typename VALUE, int DIM>
class EventOutputImage_: public EventOutput_<COO_VALUE, VALUE, DIM> {
public:
	typedef EventOutput_<COO_VALUE, VALUE, DIM> EventOutput;
	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Domain_<cvt, vt, EventOutput::Dim> Domain;
	typedef Domain_<cvt, vt, EventOutput::Dim>* pDomain;
	typedef Domain_<cvt, vt, EventOutput::Dim>& ref_Domain;
	typedef const Domain_<cvt, vt, EventOutput::Dim>& const_ref_Domain;
protected:
	std::string _prefix;
	std::string _format;

public:
	EventOutputImage_(int is = -1, int ie = -1, int istep = -1, int flag = 0,
			const std::string& prefix,
			const std::string& format) :
			EventOutput(is, ie, istep, flag, nullptr) {
		this->_prefix = prefix;
		this->_format = format;
	}

	int execute(st step, vt t, int fob, pDomain pd = nullptr) {

		return -1;
	}

	~EventOutputImage_() {

	}
protected:
	std::string _get_filename(st step, vt t, int fob){
		std::stringstream sst;
		sst<<this->_prefix<<"_"<<step<<"_"<<t;
		return sst.str();
	}

	void _draw_gnuplot_jpg(std::string& filename){

	}


};



}
