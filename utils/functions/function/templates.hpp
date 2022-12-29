#ifndef TEMPLATES_H
#define TEMPLATES_H

#include "../vectoroperators.hpp"
#include "../../debug/debugdef.hpp"
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <map>
#include <initializer_list>
#include <type_traits>

template <template <typename...>  typename InputClass,typename...Types>
struct BindTemplateLeft{
    template <typename...Args>
    using ResultType = InputClass<Types...,Args...>;
};

template <template <typename...>  typename InputClass,typename...Types>
struct BindTemplateRight{
    template <typename...Args>
    using ResultType = InputClass<Args...,Types...>;
};

template <template <typename...> typename BaseType,typename T>
struct is_template_base_of_helper{
    template <typename...Args>
    static std::true_type convertable(const volatile BaseType<Args...> *);

    static std::false_type convertable(const volatile void *);

    using value_type = typename decltype(convertable(static_cast<T*>(nullptr)))::type;
};

template <template <typename...> typename BaseType,typename T>
struct is_template_base_of: is_template_base_of_helper<BaseType,T>::value_type{};

template <typename T>
struct VectorType;

template <typename T>
struct VectorType<std::vector<T>>{
    using ResultType = T;
};


#endif
