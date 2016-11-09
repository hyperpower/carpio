#ifndef _EVENT_H_
#define _EVENT_H_

#include "../carpio_define.hpp"
#include "../utility/format.h"
#include "../utility/Clock.h"

namespace carpio {
/*
 * The event class
 * This is the base class
 */

template<typename COO_VALUE, typename VALUE, int DIM>
class Event_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Domain_<cvt, vt, Dim> Domain;
	typedef Domain_<cvt, vt, Dim>* pDomain;
	typedef Domain_<cvt, vt, Dim>& ref_Domain;
	typedef const Domain_<cvt, vt, Dim>& const_ref_Domain;

	static const int START = 0x10;
	static const int END = 0x20;
	static const int BEFORE = 0x100;
	static const int AFTER = 0x200;

protected:
	/*
	 *  _flag = 0      disabled flag
	 *  _flag = -1     forward
	 *  _flag = 1      backward      (default)
	 *  _flag = 2      both run, forward and backward
	 *  _flag = -2     both do not run, forward and backward
	 */
	int _flag;
	int _istart, _iend, _istep;

	bool _has_flag(int f) const{
		return (_flag | f) == _flag ? true : false;
	}

public:
	/**
	 *  @brief Event constructor
	 *
	 *  default constructor is running event after each step
	 */

	Event_(int is = -1, int ie = -1, int istep = -1, int flag = 0) {
		_flag = flag;
		_istart = is;
		_iend = ie;
		_istep = istep;
	}

	virtual int execute(st step, vt t, int fob, pDomain pd = nullptr) {
		std::cout << "Event : execute \n";
		return -1;
	}

	virtual int flag() const {
		return 0;
	}

	virtual void set_flag(int i) {
		return;
	}

	virtual ~Event_() {
	}

	/**
	 *  @brief Dose the event need to be executed
	 *
	 *  @param [step] is step right now
	 *  @param [t]    is time right now
	 *  @param [fob]  forward or backward,
	 *                fob = -1 is before the step,
	 *                fob = 1 is after the step
	 */
	virtual bool _do_execute(st step, vt t, int fob) const {
		bool res = ((step - _istart) % _istep == 0);
		if (_iend != -1) {
			res = res && int(step) <= _iend;
		}
		if (_istart != -1) {
			res = res && int(step) >= _istart;
		}
		res = res && (fob == _flag);
		// self value
		//  &&
		return res;
	}

};

/*
 * derived class
 *
 */
template<typename COO_VALUE, typename VALUE, int DIM>
class EventOutput_: public Event_<COO_VALUE, VALUE, DIM> {
public:
	typedef Event_<COO_VALUE, VALUE, DIM> Event;
	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Domain_<cvt, vt, Event::Dim> Domain;
	typedef Domain_<cvt, vt, Event::Dim>* pDomain;
	typedef Domain_<cvt, vt, Event::Dim>& ref_Domain;
	typedef const Domain_<cvt, vt, Event::Dim>& const_ref_Domain;
protected:
	std::FILE* _f;
	int _flag_output;
public:
	EventOutput_(int is, int ie, int istep, int flag, std::FILE* f = nullptr) :
			Event(is, ie, istep, flag) {
		if (f == nullptr) {
			_flag_output = 0;
		}
		_f = f;
	}

	int flag() const {
		return _flag_output;
	}

	std::string name() const {
		return "Output";
	}

	void set_flag(int i) {
		_flag_output = i;
	}

	virtual ~EventOutput_() {

	}

};

/*
 *  derived class
 */
template<typename COO_VALUE, typename VALUE, int DIM>
class EventOutputTime_: public EventOutput_<COO_VALUE, VALUE, DIM> {
public:
	typedef EventOutput_<COO_VALUE, VALUE, DIM> EventOutput;
	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Domain_<cvt, vt, EventOutput::Dim> Domain;
	typedef Domain_<cvt, vt, EventOutput::Dim>* pDomain;
	typedef Domain_<cvt, vt, EventOutput::Dim>& ref_Domain;
	typedef const Domain_<cvt, vt, EventOutput::Dim>& const_ref_Domain;
protected:
	st _count_exe;

	std::list<vt> l_t;     //time
	std::list<st> l_step;  //step
	std::list<tick_t> l_cpu;   //cpu time
	std::list<tick_t> l_wall;  //wall time

public:
	EventOutputTime_(int is = -1, int ie = -1, int istep = -1, int flag = 0,
			std::FILE* f = nullptr) :
			EventOutput(is, ie, istep, flag, f) {
		this->_f = f;
		this->_count_exe = 0;
	}

	int execute(st step, vt t, int fob, pDomain pd = nullptr) {
		_count_exe++;
		_record_time(step, t);
		_output(step, t);
		return -1;
	}

	~EventOutputTime_() {

	}

protected:
	tick_t _d_last(const std::list<tick_t>& list) const {
		tick_t res = 0;
		if (list.size() < 2) {
			return res;
		}
		auto last = list.rbegin();
		tick_t cpu_l = (*last);
		last++;
		tick_t cpu_ll = (*last);
		return (cpu_l - cpu_ll);
	}

	void _output(const st& step, const vt& t) {
		vt d_cpu = Clock::TicksToSecondsD(_d_last(l_cpu));
		vt d_wall = Clock::TicksToSecondsD(_d_last(l_wall));
		if (this->_flag_output == 0) {
			fmt::print(
					"step: {:>8d} t: {:>8.5f} cpu: {:>8.5f} wall: {:>8.5f}\n",
					step, t, d_cpu, d_wall);
		} else {
			fmt::print(this->_f,
					"step: {:>8d} t: {:>8.5f} cpu: {:>8.5f} wall: {:>8.5f}\n",
					step, t, d_cpu, d_wall);
		}
	}

	void _record_time(const st& step, vt& t) {
		l_t.push_back(step);
		l_t.push_back(t);
		l_cpu.push_back(Clock::SystemTime());
		l_wall.push_back(Clock::Tick());
	}

	std::string name() const {
		return "OutputTime";
	}

	bool _do_execute(st step, vt t, int flag) const {
		// flag = -2 all before
		// flag =
		// flag = 2  all end
		if (this->_has_flag(flag)) {
			bool res = ((step - this->_istart) % this->_istep == 0);
			if (this->_iend != -1) {
				res = res && int(step) <= this->_iend;
			}
			if (this->_istart != -1) {
				res = res && int(step) >= this->_istart;
			}
			// self value
			//  &&
			return res;
		} else {
			return false;
		}
	}
};

template<typename COO_VALUE, typename VALUE, int DIM>
class EventFlag_: public Event_<COO_VALUE, VALUE, DIM> {
public:
	typedef Event_<COO_VALUE, VALUE, DIM> Event;
	typedef COO_VALUE cvt;
	typedef VALUE vt;

	typedef Domain_<cvt, vt, Event::Dim> Domain;
	typedef Domain_<cvt, vt, Event::Dim>* pDomain;
	typedef Domain_<cvt, vt, Event::Dim>& ref_Domain;
	typedef const Domain_<cvt, vt, Event::Dim>& const_ref_Domain;
protected:
	int _flag;
	std::string _name; // the name of flag, it is also the key

public:
	EventFlag_(const std::string& name, int f) :
			Event(-1, -1, -1) {
		this->_flag = f;
		this->_name = name;
	}

	~EventFlag_() {

	}

	bool _do_execute(st step, vt t, int fob) const {
		return false;
	}

	int execute(st step, vt t, int fob, pDomain pd = nullptr) {
		return -1;
	}

	int flag() const {
		return _flag;
	}

	void set_flag(int i) {
		_flag = i;
	}

	std::string name() const {
		return _name;
	}
};

}
#endif
