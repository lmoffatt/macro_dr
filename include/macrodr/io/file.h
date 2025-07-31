#ifndef FILE_H
#define FILE_H

#include <filesystem>
#include <string>
namespace macrodr {
namespace core {
namespace io {

struct IO {};

bool file_exists(std::string path, std::string name) {
    auto full = std::filesystem::path(path) / name;
    return std::filesystem::exists(full);
}

bool file_exists(std::string full) {
    return std::filesystem::exists(full);
}

bool is_directory(const std::string& p) {
    return std::filesystem::is_directory(std::filesystem::path(p));
}

class Path {
    std::string path;

   public:
    Path(IO, std::string const& path) : path{path} {
    }
    auto operator()() const {
        return path;
    }
};

class File {
    std::string path;
    std::string filename;

   public:
    File(IO, std::string const& path, std::string const& filename)
        : path{path}, filename{filename} {
    }
    auto operator()() const {
        return (std::filesystem::path(path) / filename).generic_string();
    }
};

}  // namespace io

}  // namespace core
}  // namespace macrodr

#endif  // FILE_H
