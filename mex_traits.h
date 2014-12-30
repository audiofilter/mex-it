#ifndef MEX_TRAITS_
#define MEX_TRAITS_

// The MIT License (MIT) Copyright (c) Tony Kirke 2014

#include <string>  
#include <algorithm>
#include <new>      
#include <cstdlib>
#include <type_traits>
#include <tuple>
#include <functional>

// Mostly generic type traits, classes, etc for use in mex-it
namespace mex_binding {

	template <bool B, class T = void> struct disable_if { typedef T type; };
	template <class T> struct disable_if<true, T> {};

	template <class Cond, class T = void> struct enable_if_cond  : public std::enable_if<Cond::value, T> {};
	template <class Cond, class T = void> struct disable_if_cond : public disable_if<Cond::value, T> {};
	
	// ----------------------------------------------------------------------------------------
	template <typename T> struct is_const_type { static const bool value = false; };
	template <typename T> struct is_const_type<const T> { static const bool value = true; };
	template <typename T> struct is_const_type<const T&> { static const bool value = true; };

	// ----------------------------------------------------------------------------------------

	// For creating a std::tuple by appending types together
	template<typename, typename> struct append_to_type_seq { };
	template<typename T, typename... Ts> struct append_to_type_seq<T, std::tuple<Ts...>>
	{
		using type = std::tuple<Ts..., T>;
	};
	// ----------------------------------------------------------------------------------------


	// Apply std::decay on elements of a tuple to create a tuple of the basic decayed types
	template <int N,typename... Args> struct function_args_tuple;
	
	template <int N, typename T> struct function_args_tuple<N,T> {
		typedef typename std::tuple<typename std::decay<T>::type> root_type;
		root_type val;
	};
	
	template <int N, typename T, typename ...InTypes> struct function_args_tuple<N, T, InTypes...> {
		typedef typename append_to_type_seq<typename std::decay<T>::type,typename function_args_tuple<N-1, InTypes...>::root_type>::type root_type;
		root_type val;
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
	///
	/// Whever a single tuple is passed in as parameter, it is
	/// automatically unpacked and the content is forwarded to the
	/// wrapped function as parameters. Any other configuration of
	/// parameters is forwarded to the wrapped function as is.
	///
	// ---------------------------------------------------------------------- //
	///
	/// A function that unpacks its second argument and calls its
	/// first argument with the unpacked parameter list. The third
	/// parameter is needed for unpacking the second.
	///
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
	
	///
	/// A wrapper for a function to provide the automatic unpacking of
	/// a single tuple into a parameter list.
	///
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

}

#endif
