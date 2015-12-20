// Copyright (C) 2012 Massachusetts Institute of Technology, Lincoln Laboratory
// License: Boost Software License   See LICENSE.txt for the full license.
// Authors: Davis E. King (davis@dlib.net)

// Copyright (C) 2012 Massachusetts Institute of Technology, Lincoln Laboratory
// License: Boost Software License   See LICENSE.txt for the full license.
// Authors: Davis E. King (davis@dlib.net)

#include <typeinfo> // std::bad_cast
#include <vector>
#include <string>  
#include <algorithm>
#include <new>      
#include <cstdlib>
#include <type_traits>
#include <tuple>
#include <functional>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstring> // for memcpy linux
#include <cstdint> // for MSVC int8_t, etc

#if defined(_MSC_VER)
#define DLL_EXPORT_SYM __declspec(dllexport)
#endif

#include "mex.h"

namespace mex_binding {

	template <bool B, class T = void> struct disable_if { typedef T type; };
	template <class T> struct disable_if<true, T> {};

	template <class Cond, class T = void> struct enable_if_cond  : public std::enable_if<Cond::value, T> {};
	template <class Cond, class T = void> struct disable_if_cond : public disable_if<Cond::value, T> {};
	
	// ----------------------------------------------------------------------------------------
	// This is slightly different than std::is_const, since it is true for const T& also
	// this is what we want to grepping input types
	template <typename T> struct is_const_type { static const bool value = false; };
	template <typename T> struct is_const_type<const T> { static const bool value = true; };
	template <typename T> struct is_const_type<const T&> { static const bool value = true; };

	// ----------------------------------------------------------------------------------------

	// For creating a std::tuple by appending types together
	template<typename, typename> struct append_to_type_seq { };
	template<typename T, typename... Ts> struct append_to_type_seq<T, std::tuple<Ts...>>
	{
		using type = std::tuple<T,Ts...>;
	};

	// ----------------------------------------------------------------------------------------

	// Apply std::decay on elements of a tuple to create a tuple of the basic decayed types
	template <int N,typename... Args> struct function_args_tuple;
	
	template <int N, typename T> struct function_args_tuple<N,T> {
		typedef typename std::tuple<typename std::decay<T>::type> root_type;
	};
	
	template <int N, typename T, typename ...InTypes> struct function_args_tuple<N, T, InTypes...> {
		typedef typename append_to_type_seq<typename std::decay<T>::type,typename function_args_tuple<N-1, InTypes...>::root_type>::type root_type;
	};
	// ----------------------------------------------------------------------------------------

	// Index sequence stuff -------------------------------------------------------------------
  template<int...> struct seq {};
  //  for generating a sequence of indices.
  template<int N, int... S> struct gen_seq : gen_seq<N-1, N-1, S...> {};
  template<int... S> struct gen_seq<0, S...> { typedef seq<S...> type; };

  // Constructs a sequence of indices for a variadic pack of template parameters.
  template<typename... T> inline auto make_seq() -> typename gen_seq<sizeof...(T)>::type {
    return typename gen_seq<sizeof...(T)>::type();
  }
	// ----------------------------------------------------------------------------------------

  // Constructs a sequence of indices for a tuple.
  template<typename... T>
  inline auto make_seq(const std::tuple<T...>&) -> typename gen_seq<sizeof...(T)>::type {
    return typename gen_seq<sizeof...(T)>::type();
  }

	//---------------------------------------------------------------------------------------------------

	/// Automatic unpacking of function parameters.
	// Next 3 functions go together. There are also other ways to do this but 
	// this is from github/funtup library
	//
	// Whever a single tuple is passed in as parameter, it is
	// automatically unpacked and the content is forwarded to the
	// wrapped function as parameters. Any other configuration of
	// parameters is forwarded to the wrapped function as is.
	//
	// ---------------------------------------------------------------------- //
	//
	// A function that unpacks its second argument and calls its
	// first argument with the unpacked parameter list. The third
	// parameter is needed for unpacking the second.
	//
	template <typename Func, typename Args, int... I>
	inline auto unpack_and_apply(  /// The function to call
															 Func &&func,
															 /// A packed representation of the parameters to
															 /// call the function with
															 Args &&args,
															 /// An index sequence needed to unpack the
															 /// parameters
															 seq<I...> args_s) -> decltype(func(std::get<I>(args)...)) {
		return func(std::get<I>(args)...);
	}
	
	//
	// A wrapper for a function to provide the automatic unpacking of
	// a single tuple into a parameter list.
	//
	template <typename Func> class apply_unpack_t {
	public:
		inline apply_unpack_t(Func &&func) : func_m(std::forward<Func>(func)) {}
		template <typename... Args>
		inline auto operator()(Args &&... args) const -> decltype(std::declval<Func>()(std::forward<Args>(args)...)) {
			return func_m(std::forward<Args>(args)...);
		}
		template <typename... Args>
		inline auto operator()(std::tuple<Args...> &args) const
      -> decltype(unpack_and_apply(std::declval<Func>(), args, make_seq(args))) {
			return unpack_and_apply(func_m, args, make_seq(args));
		}
		template <typename... Args>
		inline auto operator()(const std::tuple<Args...> &args) const
      -> decltype(unpack_and_apply(std::declval<Func>(), args, make_seq(args))) {
			return unpack_and_apply(func_m, args, make_seq(args));
		}
		template <typename... Args>
		inline auto operator()(std::tuple<Args...> &&args) const
      -> decltype(unpack_and_apply(std::declval<Func>(), args, make_seq(args))) {
			return unpack_and_apply(func_m, args, make_seq(args));
		}
		
	private:
		Func func_m;
	};  // apply_unpack_t
	
	// Transforms a function so that it automatically unpacks a tuple into a list of
	// parameters.
	// Useful when a function that returns a tuple needs to be piped to a function
	// taking multiple arguments.
	template <typename Func> inline apply_unpack_t<Func> auto_unpack(Func &&func) {
		return apply_unpack_t<Func>(std::forward<Func>(func));
	}

// The MIT License (MIT) Copyright (c) Tony Kirke 2014
// Nothing Matlab/Mex specific here but these are all helper classes and functions
// for the mex process using Variadic Templates, etc


	struct default_is_kind_value { static const bool value = false; };
	// ----------------------------------------------------------------------------------------

	template <typename T> struct is_std_vector : public default_is_kind_value      {    };
	template <typename T, typename alloc> struct is_std_vector<std::vector<T,alloc> >         { const static bool value = true; };
	template <typename T> struct is_std_vector<T&>      { const static bool value = is_std_vector<T>::value; };
	template <typename T> struct is_std_vector<const T&>{ const static bool value = is_std_vector<T>::value; };
	template <typename T> struct is_std_vector<const T> { const static bool value = is_std_vector<T>::value; };

	template <typename T> struct is_eigen_vector : public default_is_kind_value      {    };
#ifdef EIGEN_MAJOR_VERSION
	template <typename T> struct is_eigen_vector<Eigen::Matrix<T,Eigen::Dynamic,1> >  {
    typedef T type;
		const static bool value = true;
  };
	template <typename T> struct is_eigen_vector<T&>      { const static bool value = is_eigen_vector<T>::value; };
	template <typename T> struct is_eigen_vector<const T&>{ const static bool value = is_eigen_vector<T>::value; };
	template <typename T> struct is_eigen_vector<const T> { const static bool value = is_eigen_vector<T>::value; };
#endif

	// ----------------------------------------------------------------------------------------
	template <typename T, typename helper = void> struct is_matrix : public default_is_kind_value      {
		static_assert(std::is_same<helper,void>::value, "Need same types!");
	};

	template <typename T, typename helper = void> struct is_eigen_matrix : public default_is_kind_value      {
	};


	// ----------------------------------------------------------------------------------------
	// Array2d or Array are obsolete / not supported
	template <typename T>    struct is_array2d : public default_is_kind_value      {    };
	template <typename T>    struct is_array : public default_is_kind_value      {    };
	// true if T is std::vector or array
	template <typename T>    struct is_array_type {	const static bool value = is_std_vector<T>::value || is_array<T>::value
#ifdef EIGEN_MAJOR_VERSION
			|| is_eigen_vector<T>::value
#endif
			; };
	template <typename T>    struct is_pair : public default_is_kind_value      {    };

	// ----------------------------------------------------------------------------------------



	// For checking if mex_function's arguments are either inputs or outputs (based on type)
	// ----------------------------------------------------------------------------------------

	template <typename T> struct is_input_type {
		const static unsigned long value =
			(!std::is_same<void, T>::value && (!std::is_reference<T>::value || is_const_type<T>::value)) ? 1 : 0;
	};

	template <typename T> struct is_output_type {
		const static unsigned long value =
			(!std::is_same<void, T>::value && std::is_reference<T>::value && !is_const_type<T>::value) ? 1 : 0;
	};

	// ----------------------------------------------------------------------------------------

	// Variadic template to calculate number of input types used ------------------------------------

	template <typename... Args> struct get_num_inputs;
 
	template <typename T> struct get_num_inputs<T> {
		enum {value = is_input_type<T>::value};
	};
 
	template <typename T, typename... Args> struct get_num_inputs<T, Args...> {
		enum {value = is_input_type<T>::value + get_num_inputs<Args...>::value};
	};

	// ----------------------------------------------------------------------------------------

	// Variadic template to calculate number of output types used -----------------------------

	template <typename... Args> struct get_num_outputs;
 
	template <typename T> struct get_num_outputs<T> {
		enum {value = is_output_type<T>::value};
	};
 
	template <typename T, typename... Args> struct get_num_outputs<T, Args...> {
		enum {value = is_output_type<T>::value + get_num_outputs<Args...>::value};
	};

	// ----------------------------------------------------------------------------------------

	// Variadic template to work on a function pointer to garner information about the function
	// not all functions are used here
	
	template <typename T> struct function_traits;

	template <typename R, typename... Args> struct function_traits<R(Args...)> {
		static const size_t nargs = sizeof...(Args);

		typedef R result_type;
		typedef std::tuple<Args...> ttype;

		template <size_t i> struct arg { typedef typename std::tuple_element<i, std::tuple<Args...>>::type type; };

		template <size_t i> struct is_input {
			const static unsigned long value = is_input_type<typename std::tuple_element<i, std::tuple<Args...>>::type>::value;
		};

		template <size_t i> struct is_output {
			const static unsigned long value = is_output_type<typename std::tuple_element<i, std::tuple<Args...>>::type>::value;
		};

		template <size_t i> struct stripped {
			typedef typename std::decay<typename std::tuple_element<i,std::tuple<Args...>>::type>::type type;
		};

		const static unsigned long num_args = sizeof...(Args);
		const static int get_number_of_inputs() { return get_num_inputs<Args...>::value; }
		const static int get_number_of_outputs() { return get_num_outputs<Args...>::value; }
	};

	// ----------------------------------------------------------------------------------------

	// For getting the data type within a Matrix/Array2d/Array/Std::Vector

	template <typename T, typename enabled = void> struct inner_type { typedef T type; };
	
	template <typename T>
	struct inner_type<T,
										typename std::enable_if<is_matrix<T>::value || is_array2d<T>::value || is_array<T>::value>::type> {
		typedef typename T::type type;
	};

	template <typename T> struct inner_type<T, typename enable_if_cond<is_std_vector<T>>::type> {
		typedef typename T::value_type type;
	};
	
	// ----------------------------------------------------------------------------------------

	class bad_any_cast : public std::bad_cast {
	public:
		virtual const char* what() const throw() { return "bad_any_cast"; }
	};

	struct invalid_args_exception {
		invalid_args_exception(const std::string &msg_) : msg(msg_) {}
		std::string msg;
	};

	// ----------------------------------------------------------------------------------------

	// function handle for matlab.  
	struct function_handle {
    function_handle():h(0){}
    void* const h;
	};


// Some original code from here
// Copyright (C) 2012 Massachusetts Institute of Technology, Lincoln Laboratory
// License: Boost Software License  
// Authors: Davis E. King (davis@dlib.net)

// C++11 additions and Variadic Template support
// Copyright (c) 2014 Tony Kirke

	template <typename T>	void populate_to_eigen_mat(const long arg_idx, T& m, const mxArray* src, long nc, long nr) {
		// generic - should never happen since below specialization should be used for Eigen::Matrix
		std::ostringstream sout;
		sout << "mex_function has some bug in it related to processing input argument " << arg_idx + 1 << " on line " << __LINE__ << "\n";
		mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", sout.str().c_str());
	}


#ifdef EIGEN_MAJOR_VERSION
	template <typename T> struct inner_type<T, typename enable_if_cond<is_eigen_vector<T>>::type> {
		typedef typename is_eigen_vector<T>::type type;
	};

	template <typename T> struct is_eigen_matrix<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > {
		typedef T type;
		static const bool value = true;
	};
  
	template <typename T>	struct inner_type<T,typename std::enable_if<is_eigen_matrix<T>::value>::type> {
		typedef typename is_eigen_matrix<T>::type type;
	};

	template <>	void populate_to_eigen_mat(const long arg_idx, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& m, 
																				 const mxArray* src, long nc, long nr) {
		assert(nr > 0 && nc > 0);
		m.resize(nr,nc);
		memcpy(m.data(),mxGetPr(src),nr*nc*sizeof(double));
	}

	template <typename T>
	void populate_to_eigen_mat(const long arg_idx, Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>& m, 
													 const mxArray* src, long nc, long nr) {
		assert(nr > 0 && nc > 0);
		m.resize(nr,nc);
		memcpy(m.data(),(const T*)mxGetData(src),nr*nc*sizeof(double));
	}
#endif

	template <typename T>	void populate_to_vector(const long arg_idx, T& m, const mxArray* src, long nc) {
		// generic - should never happen since below specialization should be used for std::vector
		std::ostringstream sout;
		sout << "mex_function has some bug in it related to processing input argument " << arg_idx + 1 << " on line " << __LINE__ << "\n";
		mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", sout.str().c_str());
	}

	template <>	void populate_to_vector(const long arg_idx, std::vector<double>& m, 
																					const mxArray* src, long nc) {
		assert(nc > 0);
		m.resize(nc);
		memcpy(m.data(),mxGetPr(src),nc*sizeof(double));
	}

#ifdef EIGEN_MAJOR_VERSION
	template <typename T>
	void populate_to_vector(const long arg_idx, Eigen::Matrix<T,Eigen::Dynamic,1>& m, 
													const mxArray* src, long nc) {
		assert(nc > 0);
		m.resize(nc);
		memcpy(m.data(),(const T*)mxGetPr(src),nc*sizeof(double));
	}
#endif
	
	template <typename T>
	void populate_to_vector(const long arg_idx, std::vector<T>& m, 
															const mxArray* src, long nc) {
		assert(nc > 0);
		m.resize(nc);
		memcpy(m.data(),(const T*)mxGetData(src),nc*sizeof(T));
	}
	
	// -------------------------------------------------------

    void assign_function_handle(const long arg_idx, function_handle &dest, const mxArray *src) {
        const_cast<void *&>(dest.h) = (void *)src;
    }

    template <typename T> void assign_function_handle(const long arg_idx, T &, const mxArray *) {
        std::ostringstream sout;
				sout << "mex_function has some bug in it related to processing input argument " << arg_idx + 1 << " on line " << __LINE__ << "\n";
        mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", sout.str().c_str());
    }

    // ----------------------------------------------------------------------------------------

    template <typename T> typename enable_if_cond<is_matrix<T>>::type assign_to_matlab(mxArray *&plhs, const T &item) {
			typedef typename is_matrix<T>::type type;
			//typedef typename T::type type;

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

   template <typename T> typename enable_if_cond<is_eigen_matrix<T>>::type assign_to_matlab(mxArray *&plhs, const T &item) {
		 typedef typename is_eigen_matrix<T>::type type;
        type *mat = 0;
        if (std::is_same<double, type>::value) {
            plhs = mxCreateDoubleMatrix(item.rows(), item.cols(), mxREAL);
            mat = (type *)mxGetPr(plhs);
        } else if (std::is_same<float, type>::value) {
            plhs = mxCreateNumericMatrix(item.rows(), item.cols(), mxSINGLE_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<bool, type>::value) {
            plhs = mxCreateLogicalMatrix(item.rows(), item.cols());
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<uint8_t, type>::value) {
            plhs = mxCreateNumericMatrix(item.rows(), item.cols(), mxUINT8_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<int8_t, type>::value) {
            plhs = mxCreateNumericMatrix(item.rows(), item.cols(), mxINT8_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<int16_t, type>::value ||
                   (std::is_same<short, type>::value && sizeof(short) == sizeof(int16_t))) {
            plhs = mxCreateNumericMatrix(item.rows(), item.cols(), mxINT16_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<uint16_t, type>::value ||
                   (std::is_same<unsigned short, type>::value && sizeof(unsigned short) == sizeof(uint16_t))) {
            plhs = mxCreateNumericMatrix(item.rows(), item.cols(), mxUINT16_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<int32_t, type>::value ||
                   (std::is_same<long, type>::value && sizeof(long) == sizeof(int32_t))) {
            plhs = mxCreateNumericMatrix(item.rows(), item.cols(), mxINT32_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<uint32_t, type>::value ||
                   (std::is_same<unsigned long, type>::value && sizeof(unsigned long) == sizeof(uint32_t))) {
            plhs = mxCreateNumericMatrix(item.rows(), item.cols(), mxUINT32_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<uint64_t, type>::value ||
                   (std::is_same<unsigned long, type>::value && sizeof(unsigned long) == sizeof(uint64_t))) {
            plhs = mxCreateNumericMatrix(item.rows(), item.cols(), mxUINT64_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else if (std::is_same<int64_t, type>::value ||
                   (std::is_same<long, type>::value && sizeof(long) == sizeof(int64_t))) {
            plhs = mxCreateNumericMatrix(item.rows(), item.cols(), mxINT64_CLASS, mxREAL);
            mat = (type *)mxGetData(plhs);
        } else {
            mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", "mex_function uses unsupported output argument type");
        }

        for (long c = 0; c < item.cols(); ++c) {
            for (long r = 0; r < item.rows(); ++r) {
                *mat++ = item(r, c);
            }
        }
    }

	// ----------------------------------------------------------------------------------------

	void assign_to_matlab(mxArray *&plhs, const std::string &item) { plhs = mxCreateString(item.c_str()); }

	template <typename T>
	typename disable_if<is_eigen_matrix<T>::value || is_matrix<T>::value || is_array_type<T>::value || std::is_same<T, function_handle>::value>::type
	assign_to_matlab(mxArray *&plhs, const T &item) {
		plhs = mxCreateDoubleScalar(item);
	}

	template <typename T> typename enable_if_cond<is_array_type<T>>::type assign_to_matlab(mxArray *&plhs, const T &item) {
		//mwSize dims[1] = {static_cast<mwSize>(item.size())};
		//std::cout << "In " << __FILE__ << " at line " << __LINE__ << " size = " << item.size() << "\n";
		plhs = mxCreateDoubleMatrix(item.size(),1, mxREAL);
		typedef double type;
		type *mat = (type *)mxGetPr(plhs);
		for (unsigned long i = 0; i < item.size(); ++i) {
			*mat++ = item[i];
		}
	}



	void assign_to_matlab(mxArray *&plhs, const char *str) { assign_to_matlab(plhs, std::string(str)); }
	void assign_to_matlab(mxArray *&plhs, const function_handle &h) {}

	// ----------------------------------------------------------------------------------------

	template <typename T, typename U>
	typename std::enable_if<std::is_arithmetic<T>::value || std::is_same<T, bool>::value>::type
  assign_scalar(const long arg_idx, T &dest, const U &src) {
		if (std::is_signed<U>::value && src < 0 && std::is_unsigned<T>::value) {
			std::ostringstream sout;
			sout << "Error, input argument " << arg_idx + 1 << " must be a non-negative number.";
			mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", sout.str().c_str());
		} else {
			dest = (T)src;
		}
	}
	
	template <typename T, typename U>
	typename disable_if<std::is_arithmetic<T>::value || std::is_same<T, bool>::value>::type
  assign_scalar(const long arg_idx, T &, const U &) {
		std::ostringstream sout;
		sout << "mex_function has some bug in it related to processing input argument " << arg_idx + 1 << " on line " << __LINE__ << "\n";
		mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", sout.str().c_str());
	}

	// ----------------------------------------------------------------------------------------

	template <typename T> void validate_and_populate_arg(long arg_idx, const mxArray *prhs, T &arg) {
		if (std::is_arithmetic<T>::value || std::is_same<T, bool>::value) {
      std::ostringstream sout;
			if (!(mxIsDouble(prhs) ||  mxIsSingle(prhs) ||
            (mxGetClassID(prhs) == mxINT8_CLASS) ||
            (mxGetClassID(prhs) == mxINT16_CLASS) ||
            (mxGetClassID(prhs) == mxINT32_CLASS) ||
            (mxGetClassID(prhs) == mxINT64_CLASS) ||
            (mxGetClassID(prhs) == mxUINT8_CLASS) ||
            (mxGetClassID(prhs) == mxUINT16_CLASS) ||
            (mxGetClassID(prhs) == mxUINT32_CLASS) ||
            (mxGetClassID(prhs) == mxUINT64_CLASS) ||
            mxIsLogical(prhs)) || 
          mxIsComplex(prhs) || mxGetNumberOfElements(prhs) != 1) {
				sout << " argument " << arg_idx + 1 << " must be a scalar, type is " << mxGetClassName(prhs) << "\n";
				throw invalid_args_exception(sout.str());
			}
      if (std::is_same<T, uint8_t>::value && (mxGetClassID(prhs) != mxUINT8_CLASS)) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected uint8_t\n";
				throw invalid_args_exception(sout.str());
      } else if (std::is_same<T, uint16_t>::value && (mxGetClassID(prhs) != mxUINT16_CLASS)) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected uint16_t\n";
				throw invalid_args_exception(sout.str());
      } else if (std::is_same<T, uint32_t>::value && (mxGetClassID(prhs) != mxUINT32_CLASS)) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected uint32_t\n";
				throw invalid_args_exception(sout.str());
      } else if (std::is_same<T, uint64_t>::value && (mxGetClassID(prhs) != mxUINT64_CLASS)) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected uint64_t\n";
				throw invalid_args_exception(sout.str());
      } else if (std::is_same<T, int8_t>::value && (mxGetClassID(prhs) != mxINT8_CLASS)) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected uint8_t\n";
				throw invalid_args_exception(sout.str());
      } else if (std::is_same<T, int16_t>::value && (mxGetClassID(prhs) != mxINT16_CLASS)) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected uint16_t\n";
				throw invalid_args_exception(sout.str());
      } else if (std::is_same<T, int32_t>::value && (mxGetClassID(prhs) != mxINT32_CLASS)) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected uint32_t\n";
				throw invalid_args_exception(sout.str());
      } else if (std::is_same<T, int64_t>::value && (mxGetClassID(prhs) != mxINT64_CLASS)) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected uint64_t\n";
				throw invalid_args_exception(sout.str());
      } else if (std::is_same<T, double>::value && (!mxIsDouble(prhs))) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected double\n";
				throw invalid_args_exception(sout.str());
      } else if (std::is_same<T, float>::value && (!mxIsSingle(prhs))) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected float\n";
				throw invalid_args_exception(sout.str());
      } else if (std::is_same<T, bool>::value && (!mxIsLogical(prhs))) {
				sout << " argument " << arg_idx + 1 << " types don't match!, input is " << mxGetClassName(prhs) << " expected bool\n";
				throw invalid_args_exception(sout.str());
      }
      // will cast to (T) type
			assign_scalar(arg_idx, arg, mxGetScalar(prhs));
		} else if (is_array_type<T>::value) {
			auto nr = mxGetM(prhs);
			auto nc = mxGetN(prhs);
			typedef typename inner_type<T>::type type;
			if (nr != 1 && nc != 1) {
				std::ostringstream sout;
				sout << " argument " << arg_idx + 1 << " must be a 1-D matrix (got a " << nr << "*" << nc  << " matrix)";
				throw invalid_args_exception(sout.str());
			}
			const long len = (long)std::max(nr,nc);
			if (std::is_same<type, double>::value) {
				if (!mxIsDouble(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of doubles";
					throw invalid_args_exception(sout.str());
				}
			} else if (std::is_same<type, float>::value) {
				if (!mxIsSingle(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of single/float";
					throw invalid_args_exception(sout.str());
				}
			} else if (std::is_same<type, bool>::value) {
				if (!mxIsLogical(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of logical elements.";
					throw invalid_args_exception(sout.str());
				}
			} else if (std::is_same<type, uint8_t>::value) {
				if (!mxIsUint8(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of uint8";
					throw invalid_args_exception(sout.str());
				}
			} else if (std::is_same<type, int8_t>::value) {
				if (!mxIsInt8(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of int8";
					throw invalid_args_exception(sout.str());
				}
			} else if (std::is_same<type, int16_t>::value ||
								 (std::is_same<type, short>::value && sizeof(short) == sizeof(int16_t))) {
				if (!mxIsInt16(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of int16";
					throw invalid_args_exception(sout.str());
				}
			} else if (std::is_same<type, uint16_t>::value ||
								 (std::is_same<type, unsigned short>::value && sizeof(unsigned short) == sizeof(uint16_t))) {
				if (!mxIsUint16(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of uint16";
					throw invalid_args_exception(sout.str());
				}
			} else if (std::is_same<type, int32_t>::value ||
								 (std::is_same<type, int>::value && sizeof(int) == sizeof(int32_t)) ||
								 (std::is_same<type, long>::value && sizeof(long) == sizeof(int32_t))) {
				if (!mxIsInt32(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of int32";
					throw invalid_args_exception(sout.str());
				}
			} else if (std::is_same<type, uint32_t>::value ||
								 (std::is_same<type, unsigned int>::value && sizeof(unsigned int) == sizeof(uint32_t)) ||
								 (std::is_same<type, unsigned long>::value && sizeof(unsigned long) == sizeof(uint32_t))) {
				if (!mxIsUint32(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of uint32";
					throw invalid_args_exception(sout.str());
				}
			} else if (std::is_same<type, uint64_t>::value ||
								 (std::is_same<type, unsigned int>::value && sizeof(unsigned int) == sizeof(uint64_t)) ||
								 (std::is_same<type, unsigned long>::value && sizeof(unsigned long) == sizeof(uint64_t))) {
				if (!mxIsUint64(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of uint64";
					throw invalid_args_exception(sout.str());
				}
			} else if (std::is_same<type, int64_t>::value ||
								 (std::is_same<type, int>::value && sizeof(int) == sizeof(int64_t)) ||
								 (std::is_same<type, long>::value && sizeof(long) == sizeof(int64_t))) {
				if (!mxIsInt64(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a vector of int64";
					throw invalid_args_exception(sout.str());
				}
      } else {
        std::ostringstream sout;
        sout << " argument " << arg_idx + 1 << " must be a vector of a pod type";
        throw invalid_args_exception(sout.str());
      }
			populate_to_vector(arg_idx, arg, prhs, len);
		} else if (is_eigen_matrix<T>::value) {
			typedef typename inner_type<T>::type type;
			const int num_dims = mxGetNumberOfDimensions(prhs);
			const long nr = (long)mxGetM(prhs);
			const long nc = (long)mxGetN(prhs);

			if (num_dims != 2) {
				std::ostringstream sout;
				sout << " argument " << arg_idx + 1 << " must be a 2-D matrix (got a " << num_dims << "-D matrix)";
				throw invalid_args_exception(sout.str());
			}
			
			if (std::is_same<type, double>::value) {
				if (!mxIsDouble(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of doubles";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else if (std::is_same<type, float>::value) {
				if (!mxIsSingle(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of single/float";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else if (std::is_same<type, bool>::value) {
				if (!mxIsLogical(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of logical elements.";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else if (std::is_same<type, uint8_t>::value) {
				if (!mxIsUint8(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of uint8";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else if (std::is_same<type, int8_t>::value) {
				if (!mxIsInt8(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of int8";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else if (std::is_same<type, int16_t>::value ||
								 (std::is_same<type, short>::value && sizeof(short) == sizeof(int16_t))) {
				if (!mxIsInt16(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of int16";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else if (std::is_same<type, uint16_t>::value ||
								 (std::is_same<type, unsigned short>::value && sizeof(unsigned short) == sizeof(uint16_t))) {
				if (!mxIsUint16(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of uint16";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else if (std::is_same<type, int32_t>::value ||
								 (std::is_same<type, int>::value && sizeof(int) == sizeof(int32_t)) ||
								 (std::is_same<type, long>::value && sizeof(long) == sizeof(int32_t))) {
				if (!mxIsInt32(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of int32";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else if (std::is_same<type, uint32_t>::value ||
								 (std::is_same<type, unsigned int>::value && sizeof(unsigned int) == sizeof(uint32_t)) ||
								 (std::is_same<type, unsigned long>::value && sizeof(unsigned long) == sizeof(uint32_t))) {
				if (!mxIsUint32(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of uint32";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else if (std::is_same<type, uint64_t>::value ||
								 (std::is_same<type, unsigned int>::value && sizeof(unsigned int) == sizeof(uint64_t)) ||
								 (std::is_same<type, unsigned long>::value && sizeof(unsigned long) == sizeof(uint64_t))) {
				if (!mxIsUint64(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of uint64";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else if (std::is_same<type, int64_t>::value ||
								 (std::is_same<type, int>::value && sizeof(int) == sizeof(int64_t)) ||
								 (std::is_same<type, long>::value && sizeof(long) == sizeof(int64_t))) {
				if (!mxIsInt64(prhs) || mxIsComplex(prhs)) {
					std::ostringstream sout;
					sout << " argument " << arg_idx + 1 << " must be a matrix of int64";
					throw invalid_args_exception(sout.str());
				}
				populate_to_eigen_mat(arg_idx, arg, prhs, nc, nr);
			} else {
				mexErrMsgIdAndTxt("mex_function:validate_and_populate_arg", "mex_function uses unsupported matrix type");
			}
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
		auto nr = mxGetM(prhs);
		auto nc = mxGetN(prhs);
		auto size = nr * nc;
		arg.resize(size + 1);
		if (mxGetString(prhs, &arg[0], static_cast<mwSize>(arg.size()))) {
			std::ostringstream sout;
			sout << " argument " << arg_idx + 1 << " encountered an error while calling mxGetString()";
			throw invalid_args_exception(sout.str());
		}
		arg.resize(size);
	}

	// Use enable_if SFINAE idiom, see here for explanation http://en.wikibooks.org/wiki/More_C++_Idioms/enable-if
	// or http://eli.thegreenplace.net/2014/sfinae-and-enable_if/
	// Go through std::tuple arguments and populate with values from matlab
	template<typename funct, std::size_t I=0, std::size_t N, typename Ts> inline typename std::enable_if< (I==N), void>::type
	validate_args(const mxArray *array[], int& arg_idx, Ts& arg) { };
	
	template<typename funct, std::size_t I=0, std::size_t N, typename Ts> inline typename std::enable_if< (I<N), void>::type
	validate_args(const mxArray *array[], int& arg_idx, Ts& arg) {
		if (function_traits<funct>::template is_input<I>::value) {
			//mexPrintf("calling v & p in loop arg_idx = %d\n",arg_idx);
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

	// ----------------------------------------------------------------------------------------

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
			// N3915 proposes `apply' for C++17
			auto mex_using_tuple = auto_unpack(mex_function);
			// Call mex_function, using tuple as inputs/outputs
			mex_using_tuple(Local_Args);

			i = 0;
			// Iterate through each argument and populate matlab with values if it is an output type
			assign_args<R(Args...),0,sizeof...(Args)>(plhs,i,Local_Args);
		}
	};

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
	//static mex_binding::mex_streambuf sb;
#endif

    mex_binding::call_mex_function(mex_function, nlhs, plhs, nrhs, prhs);
}

// ----------------------------------------------------------------------------------------
