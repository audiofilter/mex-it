// Copyright (C) 2012 Massachusetts Institute of Technology, Lincoln Laboratory
// License: Boost Software License   See LICENSE.txt for the full license.
// Authors: Davis E. King (davis@dlib.net)

// Copyright (C) 2012 Massachusetts Institute of Technology, Lincoln Laboratory
// License: Boost Software License   See LICENSE.txt for the full license.
// Authors: Davis E. King (davis@dlib.net)

#include "mex_stuff.h"
#include "call_mex.h"

namespace mex_binding {

	// ----------------------------------------------------------------------------------------
	template <typename funct> void call_mex_function(const funct &f, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
		const long expected_nrhs = function_traits<funct>::get_number_of_inputs();
		const long expected_nlhs = function_traits<funct>::get_number_of_outputs();
		//const long expected_args = expected_nrhs + expected_nlhs;
		
		/* check for proper number of arguments */
		if (nrhs > expected_nrhs || nrhs < expected_nrhs) {
			std::ostringstream sout;
			sout << "Expected between " << expected_nrhs << " and " << expected_nrhs << " input arguments, got "
					 << nrhs << ".";
			mexErrMsgIdAndTxt("mex_function:nrhs", sout.str().c_str());
		}
		
		if (nlhs > expected_nlhs) {
			std::ostringstream sout;
			sout << "Expected at most " << expected_nlhs << " output arguments, got " << nlhs << ".";
			mexErrMsgIdAndTxt("mex_function:nlhs", sout.str().c_str());
		}
		
		try {
			call_mex_helper<funct> helper;
			helper.call_wrapper(nlhs, plhs, nrhs, prhs);
		} catch (invalid_args_exception &e) {
			mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", ("Input" + e.msg).c_str());
		} catch (...) {
			mexErrMsgIdAndTxt("mex_function:error","mex XXXX error");
		}
	}
	
	// ----------------------------------------------------------------------------------------
	
	class mex_streambuf : public std::streambuf {
	public:
		mex_streambuf() {
			buf.resize(1000);
			setp(&buf[0], &buf[0] + buf.size() - 2);
			
			// make cout send data to mex_streambuf
			std::cout.rdbuf(this);
		}
		
	protected:
		int sync() {
			int num = static_cast<int>(pptr() - pbase());
			if (num != 0) {
				buf[num] = 0;  // null terminate the string
				mexPrintf("%s", &buf[0]);
				mexEvalString("drawnow");  // flush print to screen
				pbump(-num);
			}
			return 0;
		}
		
		int_type overflow(int_type c) {
			if (c != EOF) {
				*pptr() = c;
				pbump(1);
			}
			sync();
			return c;
		}
		
	private:
		std::vector<char> buf;
	};
	
}

// ----------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------

/* The gateway function called by MATLAB*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
// Only remap cout if we aren't using octave since octave already does this.
#if !defined(OCTAVE_IMPORT) && !defined(OCTAVE_API)
    // make it so cout prints to mexPrintf()
    static mex_binding::mex_streambuf sb;
#endif

    mex_binding::call_mex_function(mex_function, nlhs, plhs, nrhs, prhs);
}

// ----------------------------------------------------------------------------------------
