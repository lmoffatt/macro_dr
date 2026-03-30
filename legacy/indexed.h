#ifndef INDEXED_H
#define INDEXED_H

#include <parameters_derivative.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

#include "maybe_error.h"

namespace var {





    template <class Abstract>
std::vector<std::unique_ptr<Abstract>> clone_unique_vector(
    const std::vector<std::unique_ptr<Abstract>>& xs) {
    std::vector<std::unique_ptr<Abstract>> out;
    out.reserve(xs.size());
    for (const auto& x : xs) {
        out.push_back(x ? x->clone_unique() : nullptr);
    }
    return out;
}

struct AxisId {
    std::size_t value = 0;

    friend bool operator==(const AxisId& lhs, const AxisId& rhs) = default;
    friend bool operator<(const AxisId& lhs, const AxisId& rhs) {
        return lhs.value < rhs.value;
    }
};

struct AxisIndex {
    std::size_t value = 0;

    friend bool operator==(const AxisIndex& lhs, const AxisIndex& rhs) = default;
    friend bool operator<(const AxisIndex& lhs, const AxisIndex& rhs) {
        return lhs.value < rhs.value;
    }
};

struct AxisSize {
    std::size_t value = 0;

    friend bool operator==(const AxisSize& lhs, const AxisSize& rhs) = default;
    friend bool operator<(const AxisSize& lhs, const AxisSize& rhs) {
        return lhs.value < rhs.value;
    }
};  

struct Axis{
    AxisId m_id;
    AxisSize m_size;
    Axis(AxisId id, AxisSize n): m_id(id), m_size(n) {}
 };




struct IndexSpace;

struct Coordinate {
    std::vector<Axis>  axes;
    std::vector<AxisIndex>  indexes;
   public:
    auto& axis() {
        return axes;
    }

    auto& axis() const {
        return axes;
    }   
    auto& index() {
        return indexes;
    }   
    auto& index() const {
        return indexes;
    }   
    Coordinate(std::vector<Axis> axes, std::vector<AxisIndex> indexes) : axes(axes), indexes(indexes) {}

    auto validate() const -> Maybe_error<bool> {
        if (axes.size() != indexes.size()) {
            return error_message("coordinate axes/indexes size mismatch: ", axes.size(), " vs ",
                                 indexes.size());
        }
        for (std::size_t i = 0; i < axes.size(); ++i) {
            if (indexes[i].value >= axes[i].m_size.value) {
                return error_message("coordinate index ", indexes[i].value,
                                     " out of range for axis ", axes[i].m_id.value,
                                     " with size ", axes[i].m_size.value);
            }
            for (std::size_t j = i + 1; j < axes.size(); ++j) {
                if (axes[i].m_id == axes[j].m_id) {
                    return error_message("duplicate axis ", axes[i].m_id.value,
                                         " in coordinate");
                }
            }
        }
        return true;
    }

    Coordinate& next()
    {
        assert(axes.size() == indexes.size());
        for (std::size_t i = 0; i < indexes.size(); ++i) {
            if (indexes[i].value + 1 < axes[i].m_size.value) {
                ++indexes[i].value;
                return *this;
            } else {
                indexes[i].value = 0;
            }
        }
        // If we got here, we have wrapped around to the beginning
        return *this;
    }
    bool last() const {
        assert(axes.size() == indexes.size());
        for (std::size_t i = 0; i < indexes.size(); ++i) {
            if (indexes[i].value + 1 < axes[i].m_size.value) {
                return false;
            }
        }
        return true;
    }
    
    
    


    std::size_t flat_index() const {
        assert(axes.size() == indexes.size());
        std::size_t flat = 0;
        std::size_t multiplier = 1;
        for (std::size_t i = 0; i < indexes.size(); ++i) {
            assert (indexes[i].value < axes[i].m_size.value);
            flat += indexes[i].value * multiplier;
            multiplier *= axes[i].m_size.value;
        }
        return flat;
    }
  };





struct IndexSpace {
    std::vector<Axis> m_axes;

   bool includes(AxisId a)const {
        for (const auto& axis : m_axes) {
            if (axis.m_id == a) {
                return true;
            }
        }
        return false;
    }

    std::size_t size() const {
        std::size_t total = 1;
        for (const auto& axis : m_axes) {
            total *= axis.m_size.value;
        }
        return total;
    }

    auto validate() const -> Maybe_error<bool> {
        for (std::size_t i = 0; i < m_axes.size(); ++i) {
            for (std::size_t j = i + 1; j < m_axes.size(); ++j) {
                if (m_axes[i].m_id == m_axes[j].m_id) {
                    return error_message("duplicate axis ", m_axes[i].m_id.value,
                                         " in index space");
                }
            }
        }
        return true;
    }

    Coordinate begin() const{ 
        return Coordinate (m_axes, std::vector<AxisIndex>(m_axes.size(), AxisIndex{0}));
    }

    
    
    auto  local_coordinate(const Coordinate& coord)const->Maybe_error<Coordinate>  {
          auto valid_space = validate();
          if (!valid_space) {
            return valid_space.error();
          }
          auto valid_coord = coord.validate();
          if (!valid_coord) {
            return valid_coord.error();
          }
          std::vector<AxisIndex> my_indexes(m_axes.size());
          for (std::size_t i=0; i<m_axes.size(); ++i) {
            const auto& axis = m_axes[i];   
            bool found = false;
            for (std::size_t j=0; j<coord.axes.size(); ++j) {
                if (coord.axes[j].m_id == axis.m_id) {
                    if (coord.axes[j].m_size != axis.m_size) {
                        return error_message("axis ", axis.m_id.value,
                                             " size mismatch: expected ", axis.m_size.value,
                                             ", got ", coord.axes[j].m_size.value);
                    }
                    my_indexes[i] = coord.indexes[j];
                    found = true;
                    break;
                }
            }
            if (!found) {
                return error_message("coordinate axis " + std::to_string(axis.m_id.value) + " not found in index space");
            }
        }
        return Coordinate(m_axes, my_indexes);
    }


    auto flat_index(const Coordinate& coord) const -> Maybe_error<std::size_t> {
        auto local_coord = local_coordinate(coord);
        if (!local_coord) {
            return local_coord.error();
        }
        return local_coord.value().flat_index();
    }

};




template <class T>
class Indexed {
    IndexSpace m_index_space;
    std::vector<T> m_values;

   public:
    using value_type = T;

    Indexed() = default;
    explicit Indexed(IndexSpace index_space) : m_index_space(std::move(index_space)) {}
    Indexed(IndexSpace index_space, std::vector<T> values)
        : m_index_space(std::move(index_space)), m_values(std::move(values)) {}
    const IndexSpace& index_space() const {
        return m_index_space;
    }

    IndexSpace& index_space() {
        return m_index_space;
    }

    const std::vector<T>& values() const {
        return m_values;
    }

    std::vector<T>& values() {
        return m_values;
    }

    std::size_t size() const {
        return m_values.size();
    }

    auto validate() const -> Maybe_error<bool> {
        auto valid_space = m_index_space.validate();
        if (!valid_space) {
            return valid_space.error();
        }
        const auto expected = m_index_space.size();
        if (m_values.size() != expected) {
            return error_message("indexed values size mismatch: expected ", expected,
                                 ", got ", m_values.size());
        }
        return true;
    }

    bool empty() const {
        return m_values.empty();
    }

    void reserve(std::size_t n) {
        m_values.reserve(n);
    }


    const T& operator[](std::size_t i) const {
        return m_values[i];
    }

    T& operator[](std::size_t i) {
        return m_values[i];
    }

    Maybe_error<Coordinate> begin() const {
        auto valid = validate();
        if (!valid) {
            return valid.error();
        }
        return m_index_space.begin();
    }


    
    
    Maybe_error<std::reference_wrapper<const T>> at(const Coordinate& where) const {
        auto valid = validate();
        if (!valid) {
            return valid.error();
        }
        auto flat = m_index_space.flat_index(where);
        if (!flat) {
            return flat.error();
        }
        if (flat.value() >= m_values.size()) {
            return error_message("indexed value not found at requested coordinate");
        }
        return std::cref(m_values[flat.value()]);
    }

    Maybe_error<std::reference_wrapper<T>> at(const Coordinate& where) {
        auto valid = validate();
        if (!valid) {
            return valid.error();
        }
        auto flat = m_index_space.flat_index(where);
        if (!flat) {
            return flat.error();
        }
        if (flat.value() >= m_values.size()) {
            return error_message("indexed value not found at requested coordinate");
        }
        return std::ref(m_values[flat.value()]);
    }

    
};

template<class T>
auto get_IndexSpace(T&&/*unused*/) {
    return Nothing{};
}

template<class T>
auto get_IndexSpace(Indexed<T>& indexed) -> Maybe_error<IndexSpace> {
    auto valid = indexed.validate();
    if (!valid) {
        return valid.error();
    }
    return indexed.index_space();
}

template<class T>
auto get_IndexSpace(const Indexed<T>& indexed) -> Maybe_error<IndexSpace> {
    auto valid = indexed.validate();
    if (!valid) {
        return valid.error();
    }
    return indexed.index_space();
}

template<class T>
auto get_IndexSpace(Indexed<T>&& indexed) -> Maybe_error<IndexSpace> {
    auto valid = indexed.validate();
    if (!valid) {
        return valid.error();
    }
    return indexed.index_space();
}


template<class T> 
auto at_coordinate(const Indexed<T>& indexed, const Coordinate& where) {
    return indexed.at(where);
}

template<class T>
auto at_coordinate(Indexed<T>& indexed, const Coordinate& where) {
    return indexed.at(where);
}

template<class T>
auto at_coordinate(Indexed<T>&& indexed, const Coordinate& where) {
    return indexed.at(where);
}

template<class T> 
decltype(auto) at_coordinate(T&& unindexed, const Coordinate& /*where*/) {
    return std::forward<T>(unindexed);
}


template<class ...Args>
auto merge_IndexSpaces(const IndexSpace& a, const IndexSpace& b) -> Maybe_error<IndexSpace> {
    auto valid_a = a.validate();
    if (!valid_a) {
        return valid_a.error();
    }
    auto valid_b = b.validate();
    if (!valid_b) {
        return valid_b.error();
    }
    std::vector<Axis> out_axes = a.m_axes;
    for (const auto& axis_b : b.m_axes) {
        auto it = std::find_if(out_axes.begin(), out_axes.end(),
                               [&](const Axis& axis_a) { return axis_a.m_id == axis_b.m_id; });
        if (it == out_axes.end()) {
            out_axes.push_back(axis_b);
        } else if (it->m_size != axis_b.m_size) {
            return error_message("axis ", axis_b.m_id.value,
                                 " has conflicting sizes: ", it->m_size.value, " and ",
                                 axis_b.m_size.value);
        }
    }
    return IndexSpace{std::move(out_axes)};
}

inline auto merge_IndexSpaces(Nothing, Nothing) {
    return Nothing{};
}

inline auto merge_IndexSpaces(const IndexSpace& a, Nothing) {
    return a;
}

inline auto merge_IndexSpaces(Nothing, const IndexSpace& b) {
    return b;
}

inline auto merge_IndexSpaces(Maybe_error<IndexSpace> a, Nothing) {
    return a;
}

inline auto merge_IndexSpaces(Nothing, Maybe_error<IndexSpace> b) {
    return b;
}

inline auto merge_IndexSpaces(const IndexSpace& a, Maybe_error<IndexSpace> b)
    -> Maybe_error<IndexSpace> {
    if (!b) {
        return b.error();
    }
    return merge_IndexSpaces(a, b.value());
}

inline auto merge_IndexSpaces(Maybe_error<IndexSpace> a, const IndexSpace& b)
    -> Maybe_error<IndexSpace> {
    if (!a) {
        return a.error();
    }
    return merge_IndexSpaces(a.value(), b);
}

inline auto merge_IndexSpaces(Maybe_error<IndexSpace> a, Maybe_error<IndexSpace> b)
    -> Maybe_error<IndexSpace> {
    if (!a) {
        return a.error();
    }
    if (!b) {
        return b.error();
    }
    return merge_IndexSpaces(a.value(), b.value());
}

inline auto get_IndexSpaces() {
    return Nothing{};
}

template<class Arg>
auto get_IndexSpaces(Arg&& arg) {
    return get_IndexSpace(std::forward<Arg>(arg));
}

template<class Arg0, class Arg1, class... Args>
auto get_IndexSpaces(Arg0&& arg0, Arg1&& arg1, Args&&... args) {
    return merge_IndexSpaces(get_IndexSpace(std::forward<Arg0>(arg0)),
                             get_IndexSpaces(std::forward<Arg1>(arg1),
                                             std::forward<Args>(args)...));
}

inline auto as_Maybe_IndexSpace(IndexSpace x) {
    auto valid = x.validate();
    if (!valid) {
        return Maybe_error<IndexSpace>(valid.error());
    }
    return Maybe_error<IndexSpace>(std::move(x));
}

inline auto as_Maybe_IndexSpace(Maybe_error<IndexSpace> x) {
    return x;
}



template<class F,class ...Args>
auto apply_Index(F&& f, Args&&... args) {
    auto index_space = get_IndexSpaces(args...);
    if constexpr(std::is_same_v<decltype(index_space), Nothing>) {
        return std::invoke( std::forward<F>(f), std::forward<Args>(args)...);
    } else {
    using R= std::decay_t<std::invoke_result_t<F, decltype(get_value(at_coordinate(std::forward<Args>(args), std::declval<Coordinate>())))...>>; 
    auto maybe_index_space = as_Maybe_IndexSpace(std::move(index_space));
    if (!maybe_index_space) {
        return Maybe_error<Indexed<R>>(maybe_index_space.error());
    }
    auto actual_index_space = std::move(maybe_index_space.value());
    if (actual_index_space.size() == 0) {
        return Maybe_error<Indexed<R>>(Indexed<R>(std::move(actual_index_space), {}));
    }
    std::vector<Maybe_error<R> >  out;
    out.reserve(static_cast<std::size_t>(actual_index_space.size()));
    auto coord = actual_index_space.begin();
    while (true) {
        out.push_back(apply_on_Maybe_error(std::forward<F>(f),at_coordinate(std::forward<Args>(args), coord)...));
        if (coord.last()) {
            break;
        }
        coord.next();
    }
    auto out_valid = promote_Maybe_error(std::move(out));
    if (!out_valid) { return Maybe_error<Indexed<R>>(out_valid.error()); }  
    return Maybe_error<Indexed<R>>(
        Indexed<R>(std::move(actual_index_space), std::move(out_valid.value())));
}   
};  







}  // namespace var

#endif  // INDEXED_H
