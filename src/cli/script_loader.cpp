#include <macrodr/cli/script_loader.h>

#include <fstream>
#include <string>

namespace macrodr::cli {

Maybe_error<std::string> append_file(std::string&& acc, const std::string& path,
                                     const macrodr::io::PathResolver& resolver) {
    auto may_abs = resolver.resolve_existing(path);
    if (!may_abs) {
        return error_message(may_abs.error()());
    }
    const auto& abs = may_abs.value();
    std::ifstream f(abs);
    if (!f) {
        return error_message(abs + " does not exist or cannot be opened");
    }
    std::string line;
    while (std::getline(f, line)) {
        acc += line;
        acc.push_back('\n');
    }
    return std::move(acc);
}

Maybe_error<std::string> assemble_script(const std::vector<CliEvent>& sequence,
                                         const macrodr::io::PathResolver& resolver) {
    std::string s;
    for (const auto& ev : sequence) {
        if (ev.kind == CliEventKind::File) {
            auto appended = append_file(std::move(s), ev.payload, resolver);
            if (!appended) return appended.error();
            s = std::move(appended.value());
        } else {
            s += ev.payload;
            s.push_back('\n');
        }
    }
    return s;
}

}  // namespace macrodr::cli

