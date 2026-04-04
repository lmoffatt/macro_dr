#include <macrodr/cmd/indexed_construction.h>

namespace macrodr::cmd {

Maybe_error<var::Axis> axis(std::string name, std::vector<std::string> labels) {
    var::Axis out{var::AxisId{std::move(name)}, std::move(labels)};
    auto valid = out.validate();
    if (!valid) {
        return valid.error();
    }
    return out;
}

}  // namespace macrodr::cmd
