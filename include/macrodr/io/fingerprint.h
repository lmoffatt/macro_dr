// Stable content hashes for CSV outputs and a small JSON sidecar that
// records the fingerprint of the inputs that produced the CSV. Used by
// the read/cache layer to detect stale on-disk artifacts without having
// to re-inspect the CSV body itself.
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>

#include "maybe_error.h"
#include <macrodr/io/json/minijson.h>

namespace macrodr::io::fingerprint {

class Fingerprint {
   public:
    constexpr Fingerprint() = default;
    explicit constexpr Fingerprint(std::uint64_t v) : value_(v) {
    }
    constexpr std::uint64_t raw() const {
        return value_;
    }
    friend bool operator==(const Fingerprint& a, const Fingerprint& b) {
        return a.value_ == b.value_;
    }
    friend bool operator!=(const Fingerprint& a, const Fingerprint& b) {
        return !(a == b);
    }

    std::string to_hex() const {
        char buf[17];
        std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(value_));
        return std::string(buf);
    }

    static Maybe_error<Fingerprint> from_hex(std::string_view s) {
        if (s.size() != 16) {
            return error_message("fingerprint hex must be 16 chars, got ",
                                 std::to_string(s.size()));
        }
        std::uint64_t v = 0;
        for (char c : s) {
            int d = 0;
            if (c >= '0' && c <= '9')
                d = c - '0';
            else if (c >= 'a' && c <= 'f')
                d = c - 'a' + 10;
            else if (c >= 'A' && c <= 'F')
                d = c - 'A' + 10;
            else
                return error_message("invalid hex digit in fingerprint: ", std::string(1, c));
            v = (v << 4) | static_cast<std::uint64_t>(d);
        }
        return Fingerprint(v);
    }

   private:
    std::uint64_t value_ = 0;
};

// Streaming accumulator. Uses FNV-1a 64, with a separator byte between
// fields so that "abc"+"" cannot collide with "ab"+"c". Not cryptographic;
// sufficient for cache invalidation.
class Hasher {
   public:
    void update_bytes(const void* data, std::size_t len) {
        const auto* p = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < len; ++i) {
            value_ ^= p[i];
            value_ *= prime_;
        }
    }

    void update(std::string_view s) {
        update_bytes(s.data(), s.size());
        write_separator();
    }
    void update(const std::string& s) {
        update(std::string_view(s));
    }
    void update(const char* s) {
        update(std::string_view(s));
    }

    void update(bool v) {
        unsigned char b = v ? 1 : 0;
        update_bytes(&b, 1);
        write_separator();
    }

    void update(double d) {
        // Normalize so +0.0 / -0.0 hash the same, and all NaNs canonicalize.
        if (std::isnan(d))
            d = std::numeric_limits<double>::quiet_NaN();
        else if (d == 0.0)
            d = 0.0;
        std::uint64_t bits = 0;
        std::memcpy(&bits, &d, sizeof(bits));
        update_bytes(&bits, sizeof(bits));
        write_separator();
    }

    template <class I>
        requires(std::is_integral_v<I> && std::is_signed_v<I> && !std::is_same_v<I, bool>)
    void update(I v) {
        std::int64_t ext = static_cast<std::int64_t>(v);
        update_bytes(&ext, sizeof(ext));
        write_separator();
    }

    template <class U>
        requires(std::is_integral_v<U> && std::is_unsigned_v<U> && !std::is_same_v<U, bool>)
    void update(U v) {
        std::uint64_t ext = static_cast<std::uint64_t>(v);
        update_bytes(&ext, sizeof(ext));
        write_separator();
    }

    void update(const Fingerprint& f) {
        auto raw = f.raw();
        update_bytes(&raw, sizeof(raw));
        write_separator();
    }

    // Named field helper: hashes the name and then the value, so two
    // structs with differently-named but identically-valued fields hash
    // differently. Keeps fingerprint scheme tolerant to field reordering
    // only when caller consistently uses this helper.
    template <class T>
    void update_field(std::string_view name, const T& value) {
        update(name);
        update(value);
    }

    Fingerprint finalize() const {
        return Fingerprint(value_);
    }

   private:
    void write_separator() {
        unsigned char sep = 0xff;
        value_ ^= sep;
        value_ *= prime_;
    }
    static constexpr std::uint64_t offset_basis_ = 0xcbf29ce484222325ULL;
    static constexpr std::uint64_t prime_ = 0x100000001b3ULL;
    std::uint64_t value_ = offset_basis_;
};

template <class T>
Fingerprint hash_of(const T& value) {
    Hasher h;
    h.update(value);
    return h.finalize();
}

// Sidecar record persisted next to a CSV. `scheme` names the producing
// command / version ("simulate:v1", "calc_dlikelihood_predictions:v1") so
// readers can refuse a CSV whose fingerprint protocol they don't
// recognise. `metadata` is free-form JSON captured at write time for
// human inspection — never consulted for cache equality.
struct Sidecar {
    Fingerprint fingerprint{};
    std::string scheme;
    io::json::Json metadata = io::json::Json::object();
};

inline std::filesystem::path sidecar_path_for(const std::filesystem::path& csv_path) {
    std::filesystem::path p = csv_path;
    p += ".meta.json";
    return p;
}

inline Maybe_error<std::string> write_sidecar(const std::filesystem::path& csv_path,
                                              const Sidecar& sc) {
    auto path = sidecar_path_for(csv_path);
    io::json::Json root = io::json::Json::object();
    root["fingerprint"] = io::json::Json::string(sc.fingerprint.to_hex());
    root["scheme"] = io::json::Json::string(sc.scheme);
    root["metadata"] = sc.metadata;
    std::ofstream f(path);
    if (!f.is_open()) {
        return error_message("cannot open sidecar for write: ", path.string());
    }
    f << io::json::dump(root, 2);
    if (!f.good()) {
        return error_message("failed writing sidecar: ", path.string());
    }
    return path.string();
}

inline Maybe_error<Sidecar> read_sidecar(const std::filesystem::path& csv_path) {
    auto path = sidecar_path_for(csv_path);
    std::ifstream f(path);
    if (!f.is_open()) {
        return error_message("sidecar not found: ", path.string());
    }
    std::stringstream ss;
    ss << f.rdbuf();
    auto parsed = io::json::parse(ss.str());
    if (!parsed) {
        return parsed.error();
    }
    const auto& root = parsed.value();
    if (root.type != io::json::Json::Type::Object) {
        return error_message("sidecar root must be JSON object: ", path.string());
    }
    const auto* fp_j = root.find("fingerprint");
    if (fp_j == nullptr || fp_j->type != io::json::Json::Type::String) {
        return error_message("sidecar missing string field 'fingerprint': ", path.string());
    }
    auto fp = Fingerprint::from_hex(fp_j->str);
    if (!fp) {
        return fp.error();
    }
    Sidecar sc;
    sc.fingerprint = fp.value();
    if (const auto* sch = root.find("scheme");
        sch != nullptr && sch->type == io::json::Json::Type::String) {
        sc.scheme = sch->str;
    }
    if (const auto* md = root.find("metadata"); md != nullptr) {
        sc.metadata = *md;
    }
    return sc;
}

inline bool sidecar_fingerprint_matches(const std::filesystem::path& csv_path,
                                        const Fingerprint& expected) {
    auto sc = read_sidecar(csv_path);
    if (!sc) return false;
    return sc.value().fingerprint == expected;
}

}  // namespace macrodr::io::fingerprint
