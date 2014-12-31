#ifndef CALL_MEX_H_
#define CALL_MEX_H_

// Some original code from here
// Copyright (C) 2012 Massachusetts Institute of Technology, Lincoln Laboratory
// License: Boost Software License  
// Authors: Davis E. King (davis@dlib.net)

// C++11 additions and Variadic Template support
// Copyright (c) 2014 Tony Kirke


#if defined(_MSC_VER)
#define DLL_EXPORT_SYM __declspec(dllexport)
#endif

#include "mex.h"
#include <iostream>
#include <sstream>

namespace mex_binding {
    // -------------------------------------------------------

    void assign_function_handle(const long arg_idx, function_handle &dest, const mxArray *src) {
        const_cast<void *&>(dest.h) = (void *)src;
    }

    template <typename T> void assign_function_handle(const long arg_idx, T &, const mxArray *) {
        std::ostringstream sout;
        sout << "mex_function has some bug in it related to processing input argument " << arg_idx + 1;
        mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", sout.str().c_str());
    }

    // -------------------------------------------------------

    template <typename T>
    typename enable_if_cond<is_array_type<T>>::type assign_std_vector(const long arg_idx, T &dest, const mxArray *src) {
        const long nr = mxGetM(src);
        const long nc = mxGetN(src);

        typedef typename inner_type<T>::type type;

        if (!mxIsCell(src)) {
            std::ostringstream sout;
            sout << " argument " << arg_idx + 1 << " must be a cell array";
            throw invalid_args_exception(sout.str());
        }
        if (nr != 1 && nc != 1) {
            std::ostringstream sout;
            sout << " argument " << arg_idx + 1 << " must be a cell array with exactly "
                                                   "1 row or 1 column (i.e. a row or "
                                                   "column vector)";
            throw invalid_args_exception(sout.str());
        }

        const long size = nr * nc;
        dest.resize(size);

        for (unsigned long i = 0; i < dest.size(); ++i) {
            try {
                validate_and_populate_arg(i, mxGetCell(src, i), dest[i]);
            } catch (invalid_args_exception &e) {
                std::ostringstream sout;
                sout << "Error in argument " << arg_idx + 1 << ": element " << i + 1 << " of cell array not the expected type.\n";
                sout << "\t" << e.msg;
                throw invalid_args_exception(sout.str());
            }
        }
    }

    template <typename T> typename disable_if_cond<is_array_type<T>>::type assign_std_vector(const long arg_idx, T &, const mxArray *) {
        std::ostringstream sout;
        sout << "mex_function has some bug in it related to processing input argument " << arg_idx + 1;
        mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", sout.str().c_str());
    }

    // ----------------------------------------------------------------------------------------

    template <typename T> typename enable_if_cond<is_matrix<T>>::type assign_to_matlab(mxArray *&plhs, const T &item) {
        typedef typename T::type type;

        type *mat = 0;
        if (std::is_same<double, type>::value) {
            plhs = mxCreateDoubleMatrix(item.nr(), item.nc(), mxREAL);
            mat = (type *)mxGetPr(plhs);
        } else if (std::is_same<float, type>::value) {
            plhs = mxCreateNumericMatrix(item.nr(), item.nc(), mxSINGLE_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<bool, type>::value) {
            plhs = mxCreateLogicalMatrix(item.nr(), item.nc());
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<uint8_t, type>::value) {
            plhs = mxCreateNumericMatrix(item.nr(), item.nc(), mxUINT8_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<int8_t, type>::value) {
            plhs = mxCreateNumericMatrix(item.nr(), item.nc(), mxINT8_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<int16_t, type>::value ||
                   (std::is_same<short, type>::value && sizeof(short) == sizeof(int16_t))) {
            plhs = mxCreateNumericMatrix(item.nr(), item.nc(), mxINT16_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<uint16_t, type>::value ||
                   (std::is_same<unsigned short, type>::value && sizeof(unsigned short) == sizeof(uint16_t))) {
            plhs = mxCreateNumericMatrix(item.nr(), item.nc(), mxUINT16_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<int32_t, type>::value ||
                   (std::is_same<long, type>::value && sizeof(long) == sizeof(int32_t))) {
            plhs = mxCreateNumericMatrix(item.nr(), item.nc(), mxINT32_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<uint32_t, type>::value ||
                   (std::is_same<unsigned long, type>::value && sizeof(unsigned long) == sizeof(uint32_t))) {
            plhs = mxCreateNumericMatrix(item.nr(), item.nc(), mxUINT32_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<uint64_t, type>::value ||
                   (std::is_same<unsigned long, type>::value && sizeof(unsigned long) == sizeof(uint64_t))) {
            plhs = mxCreateNumericMatrix(item.nr(), item.nc(), mxUINT64_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<int64_t, type>::value ||
                   (std::is_same<long, type>::value && sizeof(long) == sizeof(int64_t))) {
            plhs = mxCreateNumericMatrix(item.nr(), item.nc(), mxINT64_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else {
            mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", "mex_function uses unsupported output argument type");
        }

        for (long c = 0; c < item.nc(); ++c) {
            for (long r = 0; r < item.nr(); ++r) {
                *mat++ = item(r, c);
            }
        }
    }

	// ----------------------------------------------------------------------------------------

	void assign_to_matlab(mxArray *&plhs, const std::string &item) { plhs = mxCreateString(item.c_str()); }

	/*
    template <typename T> void assign_to_matlab(mxArray *&plhs, const array2d<T> &item) {
		assign_to_matlab(plhs, array_to_matrix(item));
    }
	*/
    
	template <typename T> typename enable_if_cond<is_array_type<T>>::type assign_to_matlab(mxArray *&plhs, const T &item) {
		mwSize dims[1] = {item.size()};
		plhs = mxCreateCellArray(1, dims);
		for (unsigned long i = 0; i < item.size(); ++i) {
			mxArray *next = 0;
			assign_to_matlab(next, item[i]);
			mxSetCell(plhs, i, next);
		}
	}

	template <typename T>
	typename disable_if<is_matrix<T>::value || is_array_type<T>::value || std::is_same<T, function_handle>::value>::type
	assign_to_matlab(mxArray *&plhs, const T &item) {
		plhs = mxCreateDoubleScalar(item);
	}


	void assign_to_matlab(mxArray *&plhs, const char *str) { assign_to_matlab(plhs, std::string(str)); }
	void assign_to_matlab(mxArray *&plhs, const function_handle &h) {}

	// ----------------------------------------------------------------------------------------

	template <typename T, typename U>
	typename std::enable_if<std::is_arithmetic<T>::value || std::is_same<T, bool>::value>::type assign_scalar(const long arg_idx, T &dest, const U &src) {
		if (std::is_signed<U>::value && src < 0 && std::is_unsigned<T>::value) {
			std::ostringstream sout;
			sout << "Error, input argument " << arg_idx + 1 << " must be a non-negative number.";
			mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", sout.str().c_str());
		} else {
			dest = src;
		}
	}
	
	template <typename T, typename U>
	typename disable_if<std::is_arithmetic<T>::value || std::is_same<T, bool>::value>::type assign_scalar(
																																																				const long arg_idx, T &, const U &) {
		std::ostringstream sout;
		sout << "mex_function has some bug in it related to processing input argument " << arg_idx + 1;
		mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", sout.str().c_str());
	}
	
	// ----------------------------------------------------------------------------------------

	template <typename T> void validate_and_populate_arg(long arg_idx, const mxArray *prhs, T &arg) {
		if (std::is_arithmetic<T>::value || std::is_same<T, bool>::value) {
			if (!(mxIsDouble(prhs) || mxIsSingle(prhs) || mxIsLogical(prhs)) || mxIsComplex(prhs) ||
					mxGetNumberOfElements(prhs) != 1) {
				std::ostringstream sout;
				sout << " argument " << arg_idx + 1 << " must be a scalar";
				throw invalid_args_exception(sout.str());
			}
			assign_scalar(arg_idx, arg, mxGetScalar(prhs));
		} else if (is_array_type<T>::value) {
			assign_std_vector(arg_idx, arg, prhs);
		} else if (std::is_same<T, function_handle>::value) {
			if (!mxIsClass(prhs, "function_handle")) {
				std::ostringstream sout;
				sout << " argument " << arg_idx + 1 << " must be a function handle.";
				throw invalid_args_exception(sout.str());
			}
			assign_function_handle(arg_idx, arg, prhs);
		} else {
			mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", "mex_function uses unsupported input argument type");
		}
	}

	void validate_and_populate_arg(long arg_idx, const mxArray *prhs, std::string &arg) {
		if (!mxIsChar(prhs)) {
			std::ostringstream sout;
			sout << " argument " << arg_idx + 1 << " must be a char string";
			throw invalid_args_exception(sout.str());
		}
		const long nr = mxGetM(prhs);
		const long nc = mxGetN(prhs);
		const long size = nr * nc;
		arg.resize(size + 1);
		if (mxGetString(prhs, &arg[0], arg.size())) {
			std::ostringstream sout;
			sout << " argument " << arg_idx + 1 << " encountered an error while calling mxGetString()";
			throw invalid_args_exception(sout.str());
		}
		arg.resize(size);
	}


	// Go through std::tuple arguments and populate with values from matlab
	template<typename funct, std::size_t I=0, std::size_t N, typename Ts> inline typename std::enable_if< (I==N), void>::type
	validate_args(const mxArray *array[], int& arg_idx, Ts& arg) { };
	
	template<typename funct, std::size_t I=0, std::size_t N, typename Ts> inline typename std::enable_if< (I<N), void>::type
	validate_args(const mxArray *array[], int& arg_idx, Ts& arg) {
		if (function_traits<funct>::template is_input<I>::value) {
			validate_and_populate_arg(arg_idx, array[arg_idx], std::get<I>(arg));
			arg_idx++;
		}
		validate_args<funct,I+1,N>(array,arg_idx,arg);
	};

	// Go through std::tuple arguments and populate matlab with values from arguments
	template<typename funct, std::size_t I=0, std::size_t N, typename Ts> inline typename std::enable_if< (I==N), void>::type
	assign_args(mxArray *array[], int& arg_idx, const Ts& args) { };

	template<typename funct, std::size_t I=0, std::size_t N, typename Ts> inline typename std::enable_if< (I<N), void>::type
	assign_args(mxArray *array[], int& arg_idx, const Ts& args) {
		if (function_traits<funct>::template is_output<I>::value) {
			assign_to_matlab(array[arg_idx], std::get<I>(args));
			arg_idx++;
		}
		assign_args<funct,I+1,N>(array,arg_idx,args);
	};


	template <typename T> struct call_mex_helper;

	// This is where all of the variadic template magic happens
	template <typename R, typename... Args> struct call_mex_helper<R(Args...)> {
		void call_wrapper(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) const {

			// Strip out const, reference, etc from types so they can be used to create local variables
			// that can be passed to the various functions to either get data from matlab or send to matlab
			typedef typename function_args_tuple<sizeof...(Args),Args...>::root_type decayed_types;
			// This is a tuple of the 'raw' arguments
			decayed_types Local_Args;
			
			// Iterate through each argument and populate with values from matlab if it is an input type
			int i = 0;
			validate_args<R(Args...),0,sizeof...(Args)>(prhs,i,Local_Args);

			// This gets the parameter list from 'mex_function' and creates a type that can take the
			// parameter list as a std::tuple
			auto mex_using_tuple = auto_unpack(mex_function);
			// Call mex_function, using tuple as inputs/outputs
			mex_using_tuple(Local_Args);

			i = 0;
			// Iterate through each argument and populate matlab with values if it is an output type
			assign_args<R(Args...),0,sizeof...(Args)>(plhs,i,Local_Args);
		}
	};

// ----------------------------------------------------------------------------------------
} // end namespace
#endif
