// Environment JSON persistence: save/load variables and registry metadata
#pragma once

#include <chrono>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <macrodr/dsl/grammar_Identifier.h>
#include <macrodr/dsl/grammar_typed.h>
#include <macrodr/dsl/lexer_typed.h>
#include <macrodr/dsl/type_name.h>

#include "minijson.h"

namespace macrodr::io::json::envio {

// Helper: try to run a typed expression to a concrete value and write into JSON entry
template <class Lexer, class Compiler>
inline bool encode_supported_expr(const macrodr::dsl::Environment<Lexer, Compiler>& env,
                                  const macrodr::dsl::base_typed_expression<Lexer, Compiler>* expr,
                                  Json& out_val, std::string& out_type) {
    using namespace macrodr::dsl;
    {
        auto ptr = dynamic_cast<const typed_expression<Lexer, Compiler, double>*>(expr);
        if (ptr) {
            auto m = ptr->run(env);
            if (!m) return false;
            out_type = "double";
            out_val = Json::number(m.value());
            return true;
        }
    }
    {
        auto ptr_int = dynamic_cast<const typed_expression<Lexer, Compiler, int64_t>*>(expr);
        if (ptr_int) {
            auto m = ptr_int->run(env);
            if (!m) return false;
            out_type = "int";
            out_val = Json::number(static_cast<double>(m.value()));
            return true;
        }
    }
    {
        auto ptr_uint = dynamic_cast<const typed_expression<Lexer, Compiler, unsigned long>*>(expr);
        if (ptr_uint) {
            auto m = ptr_uint->run(env);
            if (!m) return false;
            out_type = "unsigned long";
            out_val = Json::number(static_cast<double>(m.value()));
            return true;
        }
    }
    {
        auto ptr = dynamic_cast<const typed_expression<Lexer, Compiler, int64_t>*>(expr);
        if (ptr) {
            auto m = ptr->run(env);
            if (!m) return false;
            out_type = "int64";
            out_val = Json::number(static_cast<double>(m.value()));
            return true;
        }
    }
    {
        auto ptr = dynamic_cast<const typed_expression<Lexer, Compiler, bool>*>(expr);
        if (ptr) {
            auto m = ptr->run(env);
            if (!m) return false;
            out_type = "bool";
            out_val = Json::boolean(m.value());
            return true;
        }
    }
    {
        auto ptr = dynamic_cast<const typed_expression<Lexer, Compiler, std::string>*>(expr);
        if (ptr) {
            auto m = ptr->run(env);
            if (!m) return false;
            out_type = "string";
            out_val = Json::string(m.value());
            return true;
        }
    }
    return false;
}

template <class Lexer, class Compiler>
inline Maybe_error<std::string> save_environment_json(
    const std::filesystem::path& path,
    const macrodr::dsl::Environment<Lexer, Compiler>& env,
    const macrodr::dsl::Compiler& compiler,
    bool include_functions = true,
    bool include_identifiers = true,
    std::string mode = "end") {
    using namespace macrodr::dsl;
    Json root = Json::object();
    root["version"] = Json::number(1);

    // variables
    Json vars = Json::object();
    for (auto const& id : env.list_variables()) {
        auto maybe_ptr = env.get(id);
        if (!maybe_ptr) continue;  // skip missing
        const auto* expr = maybe_ptr.value();
        Json val;
        std::string ty;
        if (encode_supported_expr(env, expr, val, ty)) {
            Json entry = Json::object();
            entry["type"] = Json::string(ty);
            entry["value"] = std::move(val);
            vars[id()] = std::move(entry);
        } else {
            Json entry = Json::object();
            entry["type"] = Json::string(expr->type_name());
            entry["value"] = Json::null();  // unsupported for now
            vars[id()] = std::move(entry);
        }
    }
    root["variables"] = std::move(vars);

    if (include_functions) {
        Json fns = Json::array();
        for (auto const& fid : compiler.list_functions()) {
            fns.arr.push_back(Json::string(fid()));
        }
        root["functions"] = std::move(fns);
    }
    if (include_identifiers) {
        Json ids = Json::array();
        for (auto const& iid : env.list_identifiers()) ids.arr.push_back(Json::string(iid()));
        root["identifiers"] = std::move(ids);
    }

    Json meta = Json::object();
    auto now = std::chrono::system_clock::now();
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(now.time_since_epoch()).count();
    meta["saved_at"] = Json::number(static_cast<double>(secs));
    meta["mode"] = Json::string(mode);
    root["meta"] = std::move(meta);

    std::error_code ec;
    std::filesystem::create_directories(path.parent_path(), ec);
    (void)ec;

    std::ofstream of(path);
    if (!of.is_open()) return error_message("could not open file for write: " + path.string());
    of << dump(root, 2);
    return std::string("saved environment to ") + path.string();
}

inline bool json_as_bool(const Json& j, bool def=false) {
    return (j.type == Json::Type::Bool) ? j.b : def;
}
inline double json_as_number(const Json& j, double def=0.0) {
    return (j.type == Json::Type::Number) ? j.num : def;
}
inline std::string json_as_string(const Json& j, std::string def={}) {
    return (j.type == Json::Type::String) ? j.str : def;
}

template <class Lexer, class Compiler>
inline Maybe_error<std::string> load_environment_json(
    const std::filesystem::path& path,
    macrodr::dsl::Environment<Lexer, Compiler>& env,
    std::string load_mode = "append") {
    using namespace macrodr::dsl;
    std::ifstream in(path);
    if (!in.is_open()) return error_message("could not open file for read: " + path.string());
    std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    auto maybe = macrodr::io::json::parse(content);
    if (!maybe) return maybe.error();
    const auto& root = maybe.value();
    const Json* vars = root.find("variables");
    if (!vars || vars->type != Json::Type::Object) return error_message("missing variables object");

    if (load_mode == "replace") {
        env.clear_identifiers();
        env.clear_variables();
    }

    std::size_t skipped = 0;

    for (auto const& kv : vars->obj) {
        const std::string& k = kv.first;
        const Json& entry = kv.second;
        if (entry.type != Json::Type::Object) continue;
        const Json* tyj = entry.find("type");
        const Json* valj = entry.find("value");
        if (!tyj || tyj->type != Json::Type::String || !valj) continue;
        std::string ty = tyj->str;

        auto may_id = to_Identifier<Lexer>(k);
        if (!may_id) continue;  // skip invalid identifiers
        auto id = may_id.value();

        if (ty == "double" && valj->type == Json::Type::Number) {
            auto lit = std::make_unique<typed_literal<Lexer, Compiler, double>>(valj->num);
            env.insert(id, std::unique_ptr<base_typed_expression<Lexer, Compiler>>(lit->clone()));
            env.push_back(id, new Identifier_compiler<Lexer, Compiler, double>(lit.release()));
        } else if ((ty == "int" || ty == "int64") && valj->type == Json::Type::Number) {
            auto lit = std::make_unique<typed_literal<Lexer, Compiler, int64_t>>(static_cast<int64_t>(valj->num));
            env.insert(id, std::unique_ptr<base_typed_expression<Lexer, Compiler>>(lit->clone()));
            env.push_back(id, new Identifier_compiler<Lexer, Compiler, int64_t>(lit.release()));
        } else if (ty == "bool" && (valj->type == Json::Type::Bool || valj->type == Json::Type::Number)) {
            bool b = (valj->type == Json::Type::Bool) ? valj->b : (valj->num != 0.0);
            auto lit = std::make_unique<typed_literal<Lexer, Compiler, bool>>(b);
            env.insert(id, std::unique_ptr<base_typed_expression<Lexer, Compiler>>(lit->clone()));
            env.push_back(id, new Identifier_compiler<Lexer, Compiler, bool>(lit.release()));
        } else if (ty == "string" && valj->type == Json::Type::String) {
            auto lit = std::make_unique<typed_literal<Lexer, Compiler, std::string>>(valj->str);
            env.insert(id, std::unique_ptr<base_typed_expression<Lexer, Compiler>>(lit->clone()));
            env.push_back(id, new Identifier_compiler<Lexer, Compiler, std::string>(lit.release()));
        } else {
            skipped += 1;
            continue;
        }
    }

    auto message = std::string("loaded environment from ") + path.string();
    if (skipped > 0) {
        message += " (skipped " + std::to_string(skipped) + " unsupported variable";
        message += (skipped == 1 ? ")" : "s)");
    }
    return message;
}

}  // namespace macrodr::io::json::envio
