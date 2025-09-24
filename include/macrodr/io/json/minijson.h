// Minimal JSON value, dump, and parser for simple persistence tasks
#pragma once

#include <cctype>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "maybe_error.h"

namespace macrodr::io::json {

struct Json {
    enum class Type { Null, Bool, Number, String, Array, Object };

    Type type = Type::Null;
    bool b = false;
    double num = 0.0;
    std::string str{};
    std::vector<Json> arr{};
    std::map<std::string, Json> obj{};

    static Json null() { return Json{}; }
    static Json boolean(bool v) {
        Json j; j.type = Type::Bool; j.b = v; return j;
    }
    static Json number(double v) {
        Json j; j.type = Type::Number; j.num = v; return j;
    }
    static Json string(std::string v) {
        Json j; j.type = Type::String; j.str = std::move(v); return j;
    }
    static Json array(std::vector<Json> v = {}) {
        Json j; j.type = Type::Array; j.arr = std::move(v); return j;
    }
    static Json object(std::map<std::string, Json> v = {}) {
        Json j; j.type = Type::Object; j.obj = std::move(v); return j;
    }

    Json& operator[](const std::string& key) { return obj[key]; }
    const Json* find(const std::string& key) const {
        auto it = obj.find(key);
        return it == obj.end() ? nullptr : &it->second;
    }
};

inline void dump_impl(const Json& j, std::string& out, int indent, int level) {
    auto ind = [&](int l){ for(int i=0;i<l;i++) out.push_back(' '); };
    switch (j.type) {
        case Json::Type::Null: out += "null"; break;
        case Json::Type::Bool: out += (j.b ? "true" : "false"); break;
        case Json::Type::Number: {
            // Minimal number formatting
            out += std::to_string(j.num);
            break;
        }
        case Json::Type::String: {
            out.push_back('"');
            for (char c : j.str) {
                if (c == '"' || c == '\\') { out.push_back('\\'); out.push_back(c); }
                else if (c == '\n') { out += "\\n"; }
                else { out.push_back(c); }
            }
            out.push_back('"');
            break;
        }
        case Json::Type::Array: {
            out.push_back('[');
            if (!j.arr.empty()) {
                if (indent > 0) out.push_back('\n');
                for (std::size_t i=0;i<j.arr.size();++i) {
                    if (indent > 0) ind((level+1)*indent);
                    dump_impl(j.arr[i], out, indent, level+1);
                    if (i + 1 < j.arr.size()) out += ",";
                    if (indent > 0) out.push_back('\n');
                }
                if (indent > 0) ind(level*indent);
            }
            out.push_back(']');
            break;
        }
        case Json::Type::Object: {
            out.push_back('{');
            if (!j.obj.empty()) {
                if (indent > 0) out.push_back('\n');
                std::size_t i = 0; const std::size_t n = j.obj.size();
                for (auto const& kv : j.obj) {
                    if (indent > 0) ind((level+1)*indent);
                    dump_impl(Json::string(kv.first), out, indent, level+1);
                    out += indent>0? ": ":":";
                    dump_impl(kv.second, out, indent, level+1);
                    if (++i < n) out += ",";
                    if (indent > 0) out.push_back('\n');
                }
                if (indent > 0) ind(level*indent);
            }
            out.push_back('}');
            break;
        }
    }
}

inline std::string dump(const Json& j, int indent = 2) {
    std::string out; out.reserve(256);
    dump_impl(j, out, indent, 0);
    return out;
}

struct Parser {
    std::string_view s;
    std::size_t i = 0;
    explicit Parser(std::string_view sv) : s(sv) {}

    void ws() { while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i; }
    bool consume(char c) { ws(); if (i < s.size() && s[i]==c){ ++i; return true;} return false; }
    bool match(std::string_view t) {
        ws(); if (s.substr(i, t.size()) == t) { i += t.size(); return true; } return false;
    }

    Maybe_error<Json> parse_value() {
        ws(); if (i >= s.size()) return error_message("unexpected end of input");
        char c = s[i];
        if (c == 'n') { if (!match("null")) return error_message("invalid null"); return Json::null(); }
        if (c == 't') { if (!match("true")) return error_message("invalid true"); return Json::boolean(true); }
        if (c == 'f') { if (!match("false")) return error_message("invalid false"); return Json::boolean(false); }
        if (c == '"') return parse_string();
        if (c == '[') return parse_array();
        if (c == '{') return parse_object();
        return parse_number();
    }

    Maybe_error<Json> parse_string() {
        if (!consume('"')) return error_message("expected \" for string");
        std::string out;
        while (i < s.size()) {
            char c = s[i++];
            if (c == '"') break;
            if (c == '\\') {
                if (i >= s.size()) return error_message("bad escape");
                char e = s[i++];
                switch (e) {
                    case '"': out.push_back('"'); break;
                    case '\\': out.push_back('\\'); break;
                    case '/': out.push_back('/'); break;
                    case 'b': out.push_back('\b'); break;
                    case 'f': out.push_back('\f'); break;
                    case 'n': out.push_back('\n'); break;
                    case 'r': out.push_back('\r'); break;
                    case 't': out.push_back('\t'); break;
                    default: return error_message("unsupported escape");
                }
            } else {
                out.push_back(c);
            }
        }
        Json j = Json::string(std::move(out));
        return j;
    }

    Maybe_error<Json> parse_number() {
        ws();
        std::size_t start = i;
        if (i < s.size() && (s[i] == '-' || s[i] == '+')) ++i;
        while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) ++i;
        if (i < s.size() && s[i] == '.') { ++i; while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) ++i; }
        if (i < s.size() && (s[i] == 'e' || s[i] == 'E')) {
            ++i; if (i < s.size() && (s[i]=='+'||s[i]=='-')) ++i; while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i]))) ++i;
        }
        if (start == i) return error_message("invalid number");
        double val = std::strtod(std::string(s.substr(start, i-start)).c_str(), nullptr);
        return Json::number(val);
    }

    Maybe_error<Json> parse_array() {
        if (!consume('[')) return error_message("expected [");
        Json j = Json::array();
        ws();
        if (consume(']')) return j;
        while (true) {
            auto v = parse_value(); if (!v) return v.error(); j.arr.emplace_back(std::move(v.value()));
            ws(); if (consume(']')) break; if (!consume(',')) return error_message("expected , or ]");
        }
        return j;
    }

    Maybe_error<Json> parse_object() {
        if (!consume('{')) return error_message("expected {");
        Json j = Json::object();
        ws();
        if (consume('}')) return j;
        while (true) {
            auto k = parse_string(); if (!k) return k.error();
            if (!consume(':')) return error_message("expected :");
            auto v = parse_value(); if (!v) return v.error();
            j.obj.emplace(k.value().str, std::move(v.value()));
            ws(); if (consume('}')) break; if (!consume(',')) return error_message("expected , or }");
        }
        return j;
    }
};

inline Maybe_error<Json> parse(std::string_view s) {
    Parser p(s);
    auto v = p.parse_value();
    if (!v) return v.error();
    p.ws();
    if (p.i != s.size()) return error_message("trailing characters after JSON value");
    return v;
}

}  // namespace macrodr::io::json

