#ifndef MEX_STUFF_H_
#define MEX_STUFF_H_

#include <typeinfo> // std::bad_cast
#include "mex_traits.h"
#include <vector>

// The MIT License (MIT) Copyright (c) Tony Kirke 2014
// Nothing Matlab/Mex specific here but these are all helper classes and functions
// for the mex process using Variadic Templates, etc

namespace mex_binding 
{

	struct default_is_kind_value { static const bool value = false; };
	// ----------------------------------------------------------------------------------------

	template <typename T> struct is_std_vector : public default_is_kind_value      {    };
	template <typename T, typename alloc> struct is_std_vector<std::vector<T,alloc> >         { const static bool value = true; };
	template <typename T> struct is_std_vector<T&>      { const static bool value = is_std_vector<T>::value; };
	template <typename T> struct is_std_vector<const T&>{ const static bool value = is_std_vector<T>::value; };
	template <typename T> struct is_std_vector<const T> { const static bool value = is_std_vector<T>::value; };


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
	template <typename T>    struct is_array_type {	const static bool value = is_std_vector<T>::value || is_array<T>::value;	};
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


} // namespace mex_binding
// -------------------------------------------------------
#endif
