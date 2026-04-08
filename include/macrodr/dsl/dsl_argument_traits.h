#pragma once

#include <functional>
#include <memory>
#include <type_traits>
#include <utility>

namespace macrodr::dsl::detail {

template <class T>
struct is_unique_ptr : std::false_type {};

template <class T, class Deleter>
struct is_unique_ptr<std::unique_ptr<T, Deleter>> : std::true_type {};

template <class Arg, class Enable = void>
struct function_argument_storage {
    using type = std::remove_cvref_t<Arg>;
};

template <class T>
struct function_argument_storage<const T&> {
    using type = std::reference_wrapper<const T>;
};

template <class T>
struct function_argument_storage<T&> {
    using type = std::reference_wrapper<T>;
};

template <class Arg>
using function_argument_storage_t = typename function_argument_storage<Arg>::type;

template <class Param, class Storage>
decltype(auto) adapt_argument_like(Storage& storage) {
    using Base = std::remove_reference_t<Param>;
    using StorageType = std::remove_reference_t<Storage>;

    if constexpr (std::is_same_v<StorageType, std::reference_wrapper<Base>> ||
                  std::is_same_v<StorageType,
                                 std::reference_wrapper<std::remove_const_t<Base>>>) {
        return static_cast<Param>(storage.get());
    } else if constexpr (std::is_lvalue_reference_v<Param>) {
        if constexpr (is_unique_ptr<StorageType>::value) {
            return static_cast<Param>(*storage);
        } else if constexpr (std::is_pointer_v<StorageType>) {
            return static_cast<Param>(*storage);
        } else {
            if constexpr (std::is_const_v<Base>) {
                return static_cast<const Base&>(storage);
            } else {
                return static_cast<Base&>(storage);
            }
        }
    } else if constexpr (std::is_rvalue_reference_v<Param>) {
        if constexpr (is_unique_ptr<StorageType>::value) {
            return static_cast<Base&&>(*storage);
        } else if constexpr (std::is_pointer_v<StorageType>) {
            return static_cast<Base&&>(*storage);
        } else {
            return static_cast<Base&&>(storage);
        }
    } else {
        if constexpr (is_unique_ptr<StorageType>::value) {
            return std::move(storage);
        } else if constexpr (std::is_pointer_v<StorageType>) {
            return storage;
        } else {
            return static_cast<Param>(std::move(storage));
        }
    }
}

} // namespace macrodr::dsl::detail