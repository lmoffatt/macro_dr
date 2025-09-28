// Environment JSON persistence: save/load variables and registry metadata
#pragma once

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <macrodr/dsl/grammar_Identifier.h>
#include <macrodr/dsl/grammar_typed.h>
#include <macrodr/dsl/lexer_typed.h>
#include <macrodr/dsl/type_name.h>
#include <macrodr/cmd/load_model.h>

#include "minijson.h"
#include "convert.h"
#include "qmodel.h"

namespace macrodr::io::json::envio {

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
    Json schemas = Json::object();
    for (auto const& id : env.list_variables()) {
        auto maybe_ptr = env.get(id);
        if (!maybe_ptr) continue;  // skip missing
        const auto* expr = maybe_ptr.value();
        Json entry = Json::object();
        try {
            auto serialized = expr->serialize_json(env, macrodr::io::json::conv::TagPolicy::Minimal);
            if (serialized) {
                const auto& data = serialized.value();
                // Normalize model pointers to a stable 'model' type with id+binary
                const std::string model_ptr_type =
                    macrodr::dsl::type_name<std::unique_ptr<macrodr::interface::IModel<var::Parameters_values>>>();
                if (data.type == model_ptr_type) {
                    entry["type"] = Json::string("model");
                    entry["value"] = data.value;  // contains {model_id,binary}
                } else {
                    entry["type"] = Json::string(data.type);
                    entry["value"] = data.value;
                }
                // Capture schemas referenced as variables
                if (data.type == macrodr::dsl::type_name<var::Parameters_Transformations>()) {
                    const auto* jid = data.value.find("id");
                    if (jid && jid->type == Json::Type::String) {
                        schemas.obj.emplace(jid->str, data.value);
                    }
                }
            } else {
                entry["type"] = Json::string(expr->type_name());
                entry["value"] = Json::null();
                entry["error"] = Json::string(serialized.error()());
            }
        } catch (const std::exception& e) {
            entry["type"] = Json::string(expr->type_name());
            entry["value"] = Json::null();
            entry["error"] = Json::string(std::string{"exception during serialization: "} + e.what());
        } catch (...) {
            entry["type"] = Json::string(expr->type_name());
            entry["value"] = Json::null();
            entry["error"] = Json::string("unknown exception during serialization");
        }
        vars[id()] = std::move(entry);
    }
    root["variables"] = std::move(vars);
    if (!schemas.obj.empty()) root["schemas"] = std::move(schemas);

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
    const auto& compiler = env.compiler();
    std::ifstream in(path);
    if (!in.is_open()) return error_message("could not open file for read: " + path.string());
    std::string content((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    auto maybe = macrodr::io::json::parse(content);
    if (!maybe) return maybe.error();
    const auto& root = maybe.value();
    // Phase 1: load schemas if present
    const Json* sch = root.find("schemas");
    if (sch && sch->type == Json::Type::Object) {
        for (const auto& [sid, sobj] : sch->obj) {
            var::Parameters_Transformations pt;
            auto st = macrodr::io::json::conv::from_json(sobj, pt, std::string("$.schemas.") + sid,
                                                         macrodr::io::json::conv::TagPolicy::None);
            if (st) {
                env.register_parameter_schema(sid, std::make_shared<var::Parameters_Transformations>(pt));
            }
        }
    }

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
        const std::string value_path = std::string("$.variables.") + k + ".value";
        Maybe_error<void> load_status;
        if (ty == "model") {
            // Load model by id using factory; store as ModelPtr
            const Json* mid = valj->find("model_id");
            if (!mid || mid->type != Json::Type::String) {
                skipped += 1;
                continue;
            }
            auto handle = macrodr::cmd::load_model(mid->str);
            if (!handle) {
                skipped += 1;
                continue;
            }
            using ModelPtr = std::unique_ptr<macrodr::interface::IModel<var::Parameters_values>>;
            ModelPtr model = std::move(handle.value());
            auto literal = std::make_unique<typed_literal<Lexer, Compiler, ModelPtr>>(std::move(model));
            env.insert(id, std::unique_ptr<base_typed_expression<Lexer, Compiler>>(literal->clone()));
            env.push_back(id, new Identifier_compiler<Lexer, Compiler, ModelPtr>(literal.release()));
            load_status = {};
        } else {
            load_status = compiler.load_variable_from_json(
                ty, *valj, value_path, macrodr::io::json::conv::TagPolicy::Minimal, id, env);
        }
        if (!load_status) {
            bool recovered = false;
            if (ty == "int" || ty == "int64") {
                load_status = compiler.load_variable_from_json(
                    macrodr::dsl::type_name<int64_t>(), *valj, value_path,
                    macrodr::io::json::conv::TagPolicy::Minimal, id, env);
                recovered = static_cast<bool>(load_status);
            } else if (ty == "string") {
                load_status = compiler.load_variable_from_json(
                    macrodr::dsl::type_name<std::string>(), *valj, value_path,
                    macrodr::io::json::conv::TagPolicy::Minimal, id, env);
                recovered = static_cast<bool>(load_status);
            }
            if (!recovered) {
                skipped += 1;
                continue;
            }
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
