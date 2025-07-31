#pragma once
#include <functional>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include "maybe_error.h"

namespace macrodr {

namespace core {
namespace df {

template <typename... Ts>
class DataFrame {
   public:
    using Column = std::variant<std::vector<Ts>...>;

   private:
    std::vector<std::string> column_names_;
    std::vector<Column> columns_;
    std::unordered_map<std::string, size_t> name_to_index_;

   public:
    bool has_column(const std::string& name) const {
        return name_to_index_.find(name) != name_to_index_.end();
    }

    auto num_rows() const {
        if (columns_.size() > 0)
            return columns_.front().size();
        else
            return 0ul;
    }

    size_t num_columns() const {
        return columns_.size();
    }

    Maybe_error<std::reference_wrapper<Column const>> column(const std::string& name) const {
        auto it = name_to_index_.find(name);
        if (it == name_to_index_.end())
            return error_message("No such column: " + name);
        return std::cref(columns_[it->second]);
    }

    const std::vector<std::string>& column_names() const {
        return column_names_;
    }
};
}  // namespace df
}  // namespace core
}  // namespace macrodr
