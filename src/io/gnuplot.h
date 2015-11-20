#ifndef _GNUPLOT_H_
#define _GNUPLOT_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>              // for std::ostringstream
#include <stdexcept>
#include <cstdio>
#include <cstdlib>              // for getenv()
#include <list>                 // for std::list

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
//defined for 32 and 64-bit environments
#include <io.h>                // for _access(), _mktemp()
#define GP_MAX_TMP_FILES  27   // 27 temporary files it's Microsoft restriction
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
//all UNIX-like OSs (Linux, *BSD, MacOSX, Solaris, ...)
#include <unistd.h>            // for access(), mkstemp()
#define GP_MAX_TMP_FILES  64
#else
#error unsupported or unknown operating system
#endif

namespace carpio {

class Gnuplot {
protected:
	/*
	 * \brief pointer to the stream that can be used to write to the pipe
	 */
	FILE *gnucmd;
	///\brief name of executed GNUPlot file
	std::string m_sGNUPlotFileName;
	///\brief gnuplot path
	std::string m_sGNUPlotPath;
	///\brief standart terminal, used by showonscreen
	std::string terminal_std;

	bool _valid;
	/*
	 * Opens up a gnuplot session, ready to receive commands
	 */
	void _init();
	/*
	 * Find out if a command lives in m_sGNUPlotPath or in PATH
	 */
	bool _get_program_path();

	/*
	 * Check if file exists
	 */
	bool _file_exists(const std::string &filename, int mode);
public:
	Gnuplot();
	~Gnuplot();

	/*
	 * Sends a command to an active gnuplot session
	 */
	Gnuplot& cmd(const std::string &cmdstr);

	/// turns grid on/off
	inline Gnuplot& set_grid() {
		cmd("set grid");
		return *this;
	}

	/// grid is not set by default
	inline Gnuplot& unset_grid() {
		cmd("unset grid");
		return *this;
	}

	// -----------------------------------------------
	/// set the mulitplot mode
	///
	/// \param ---
	///
	/// \return <-- reference to the gnuplot object
	// -----------------------------------------------
	inline Gnuplot& set_multiplot() {
		cmd("set multiplot");
		return *this;
	}

	// -----------------------------------------------
	/// unsets the mulitplot mode
	///
	/// \param ---
	///
	/// \return <-- reference to the gnuplot object
	// -----------------------------------------------
	inline Gnuplot& unset_multiplot() {
		cmd("unset multiplot");
		return *this;
	}

	// -----------------------------------------------------------------------
	/// \brief sets and clears the title of a gnuplot session
	///
	/// \param title --> the title of the plot [optional, default == ""]
	///
	/// \return <-- reference to the gnuplot object
	// -----------------------------------------------------------------------
	inline Gnuplot& set_title(const std::string &title = "") {
		std::string cmdstr;
		cmdstr = "set title \"";
		cmdstr += title;
		cmdstr += "\"";
		cmd(cmdstr);
		return *this;
	}

	//----------------------------------------------------------------------------------
	///\brief Clears the title of a gnuplot session
	/// The title is not set by default.
	///
	/// \param ---
	///
	/// \return <-- reference to the gnuplot object
	// ---------------------------------------------------------------------------------
	inline Gnuplot& unset_title() {
		this->set_title();
		return *this;
	}

	/// set x axis label
	Gnuplot& set_xlabel(const std::string &label = "x");
	/// set y axis label
	Gnuplot& set_ylabel(const std::string &label = "y");
	/// set z axis label
	Gnuplot& set_zlabel(const std::string &label = "z");

	Gnuplot& set_equal_ratio();

	//------------------------------------------------------------------------------
	//

	// set
	Gnuplot& set(const std::string& str);

	Gnuplot& set_palette_blue_red();
	// set range
	// set the xrange
	Gnuplot& set_xrange(const double iFrom, const double iTo);
	Gnuplot& set_xrange_reverse(const double iFrom, const double iTo);
	Gnuplot& set_yrange(const double iFrom, const double iTo);
	Gnuplot& set_yrange_reverse(const double iFrom, const double iTo);
	Gnuplot& set_zrange(const double iFrom, const double iTo);
	Gnuplot& set_zrange_reverse(const double iFrom, const double iTo);
	Gnuplot& set_cbrange(const double iFrom, const double iTo);

	/*
	 *  plot
	 */
	template<typename CONTAINER>
	Gnuplot& plot_1(     //
			const CONTAINER& x,  //
			const std::string &str = "") {
		std::ostringstream ss;
		ss << "plot \"-\" using 1 " << str << "\n";
		cmd(ss.str());
		ss.str("");
		for (typename CONTAINER::const_iterator it = x.begin(); it != x.end();
				++it) {
			ss << (*it) << "\n";
			cmd(ss.str());
			ss.str("");
		}
		cmd("e\n");
		return *this;
	}
	template<typename X, typename Y>
	Gnuplot& plot_2(   //  type has [] and size()
			const X& x, //
			const Y& y, //
			const std::string &str = "") {  //
		if (x.size() != y.size()) {
			std::cerr << " >Warning! The containers' size are not equal. \n";
			std::cerr << " >Warning! x =" << x.size() << " y =" << y.size << " \n";
		}
		// inline data
		std::ostringstream sst;
		//
		sst << "plot \"-\" using 1:2 " << str;
		cmd(sst.str());
		sst.str("");
		typename X::const_iterator iterx = x.begin();
		typename Y::const_iterator itery = y.begin();
		if (x.size() >= y.size()) {
			for (; itery != y.end();) {
				sst << (*iterx) << " " << (*itery);
				sst << "\n";
				cmd(sst.str());
				sst.str("");
				iterx++;
				itery++;
			}
		} else {
			for (; iterx != x.end();) {
				sst << (*iterx) << " " << (*itery);
				sst << "\n";
				cmd(sst.str());
				sst.str("");
				iterx++;
				itery++;
			}
		}
		cmd("e\n");
		return *this;
	}
}
;

//------------------------------------------------------------------------------
//
// A string tokenizer taken from http://www.sunsite.ualberta.ca/Documentation/
// /Gnu/libstdc++-2.90.8/html/21_strings/stringtok_std_h.txt
//
template<typename Container>
void stringtok(Container &container, std::string const &in,
		const char * const delimiters = " \t\n") {
	const std::string::size_type len = in.length();
	std::string::size_type i = 0;

	while (i < len) {
		// eat leading whitespace
		i = in.find_first_not_of(delimiters, i);

		if (i == std::string::npos)
			return;   // nothing left but white space

		// find the end of the token
		std::string::size_type j = in.find_first_of(delimiters, i);

		// push token
		if (j == std::string::npos) {
			container.push_back(in.substr(i));
			return;
		} else
			container.push_back(in.substr(i, j - i));
		// set up for next loop
		i = j + 1;
	}
	return;
}

}

#endif