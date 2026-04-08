
#pragma once
#include <macrodr/dsl/dsl_forward.h>


namespace macrodr::dsl {

template<class Lexer>
class Compiler {
    using Json = macrodr::io::json::Json;
    using TagPolicy = macrodr::io::json::conv::TagPolicy;

    struct TypeEntry {
        using LoadFn = Maybe_error<void> (*)(const Json&, const std::string&, TagPolicy,
                                             const Identifier<Lexer>&,
                                             Environment<Lexer, Compiler>&);
        LoadFn load = nullptr;
    };

    // Support overloaded functions: multiple compilers per identifier
    std::map<Identifier<Lexer>,
             std::vector<std::unique_ptr<base_function_compiler<Lexer, Compiler>>>>
        m_func;
    std::map<std::string, TypeEntry> m_type_registry;

   public:
    Compiler() { register_builtin_types(); }
    // Compatibility ctor: accept a single-implementation map and wrap entries as one-overload vectors
    Compiler(std::map<Identifier<Lexer>, std::unique_ptr<base_function_compiler<Lexer, Compiler>>>&&
                 func) {
        for (auto& [name, fn] : func) {
            auto& vec = m_func[name];
            std::unique_ptr<base_function_compiler<Lexer, Compiler>> moved = std::move(fn);
            if (moved) {
                moved->register_types(*this);
            }
            vec.emplace_back(std::move(moved));
        }
        register_builtin_types();
    }

    Compiler(const Compiler& cm) : m_type_registry{cm.m_type_registry} {
        // Deep-clone the overload vectors per identifier
        for (const auto& [name, vec] : cm.m_func) {
            auto& dst_vec = m_func[name];
            dst_vec.reserve(vec.size());
            for (const auto& fn : vec) {
                dst_vec.emplace_back(fn ? fn->clone() : nullptr);
            }
        }
    }

    ~Compiler() = default;

    Compiler& operator=(const Compiler& cm) {
        if (this != &cm) {
            m_type_registry = cm.m_type_registry;
            m_func.clear();
            for (const auto& [name, vec] : cm.m_func) {
                auto& dst_vec = m_func[name];
                dst_vec.reserve(vec.size());
                for (const auto& fn : vec) {
                    dst_vec.emplace_back(fn ? fn->clone() : nullptr);
                }
            }
        }
        return *this;
    }

    Compiler(Compiler&& other) noexcept = default;

    Compiler& operator=(Compiler&& other) noexcept = default;

    [[nodiscard]] Maybe_error<base_function_compiler<Lexer, Compiler> const*> get_function(
        const Identifier<Lexer>& id) const {
        auto it = m_func.find(id);
        if (it == m_func.end()) {
            return error_message(id() + " function is not defined");
        }
        if (it->second.empty() || it->second.front().get() == nullptr) {
            return error_message(id() + " function is null");
        }
        return it->second.front().get();
    }

    [[nodiscard]] Maybe_error<std::vector<base_function_compiler<Lexer, Compiler> const*>>
        get_functions(const Identifier<Lexer>& id) const {
        auto it = m_func.find(id);
        if (it == m_func.end()) {
            return error_message(id() + " function is not defined");
        }
        std::vector<base_function_compiler<Lexer, Compiler> const*> out;
        out.reserve(it->second.size());
        for (const auto& fn : it->second) {
            out.push_back(fn.get());
        }
        return out;
    }

    bool has_registered_type(const std::string& type_name) const {
        return m_type_registry.find(type_name) != m_type_registry.end();
    }

    Maybe_error<bool> push_function(std::string id_candidate,
                                    base_function_compiler<Lexer, Compiler>* fun) {
        auto may_id = to_Identifier<Lexer>(std::move(id_candidate));
        if (!may_id) {
            return may_id.error();
        }
        fun->register_types(*this);
        m_func[may_id.value()].emplace_back(fun);
        return true;
    }
    Maybe_error<bool> push_function(
        std::string id_candidate, std::unique_ptr<base_function_compiler<Lexer, Compiler>>&& fun) {
        auto may_id = to_Identifier<Lexer>(std::move(id_candidate));
        if (!may_id) {
            return may_id.error();
        }
        fun->register_types(*this);
        m_func[may_id.value()].emplace_back(std::move(fun));
        return true;
    }
    void merge(const Compiler& other) {
        // Merge functions: clone and append overloads
        for (const auto& [name, vec] : other.m_func) {
            auto& dst_vec = m_func[name];
            for (const auto& fn : vec) {
                if (fn) {
                    auto cloned =
                        std::unique_ptr<base_function_compiler<Lexer, Compiler>>(fn->clone());
                    cloned->register_types(*this);
                    dst_vec.emplace_back(std::move(cloned));
                }
            }
        }
        for (const auto& [key, entry] : other.m_type_registry) {
            if (m_type_registry.find(key) == m_type_registry.end()) {
                m_type_registry.emplace(key, entry);
            }
        }
    }
    void merge(Compiler&& other) {
        // Move functions: append overload vectors
        for (auto& [name, vec] : other.m_func) {
            auto& dst_vec = m_func[name];
            for (auto& fn : vec) {
                if (fn) {
                    fn->register_types(*this);
                }
                dst_vec.emplace_back(std::move(fn));
            }
        }
        for (auto& [key, entry] : other.m_type_registry) {
            if (m_type_registry.find(key) == m_type_registry.end()) {
                m_type_registry.emplace(key, entry);
            }
        }
    }

    // Introspection helper for environment persistence
    [[nodiscard]] std::vector<Identifier<Lexer>> list_functions() const {
        std::vector<Identifier<Lexer>> out;
        out.reserve(m_func.size());
        for (const auto& kv : m_func) out.push_back(kv.first);
        return out;
    }

   private:
    template <class Value>
    static TypeEntry make_entry() {
        return TypeEntry{
            &load_literal_from_json_helper<Lexer, Compiler, std::remove_cvref_t<Value>>};
    }

    void register_type_alias(const std::string& alias, const TypeEntry& entry) {
        if (!alias.empty() && m_type_registry.find(alias) == m_type_registry.end()) {
            m_type_registry.emplace(alias, entry);
        }
    }

    void register_builtin_types() {
        ensure_type_registered<double>();
        ensure_type_registered<int64_t>();
        ensure_type_registered<unsigned long>();
        ensure_type_registered<bool>();
        ensure_type_registered<std::string>();
        ensure_type_registered<Matrix<double>>();
        ensure_type_registered<DiagonalMatrix<double>>();

        const auto string_key = type_name<std::string>();
        auto string_it = m_type_registry.find(string_key);
        if (string_it != m_type_registry.end()) {
            register_type_alias("string", string_it->second);
        }

        const auto int_key = type_name<int64_t>();
        auto int_it = m_type_registry.find(int_key);
        if (int_it != m_type_registry.end()) {
            register_type_alias("int", int_it->second);
            register_type_alias("int64", int_it->second);
        }
    }

   public:
    template <class T>
    void ensure_type_registered() {
        using Value = std::remove_cvref_t<T>;
        if constexpr (std::is_void_v<Value>) {
            return;
        } else if constexpr (!macrodr::io::json::conv::has_json_codec_v<Value>) {
            return;
        } else {
            const std::string key = type_name<Value>();
            if (m_type_registry.find(key) == m_type_registry.end()) {
                m_type_registry.emplace(key, make_entry<Value>());
            }
        }
    }

    Maybe_error<void> load_variable_from_json(const std::string& type_name, const Json& value,
                                              const std::string& path, TagPolicy policy,
                                              const Identifier<Lexer>& id,
                                              Environment<Lexer, Compiler>& env) const {
        auto it = m_type_registry.find(type_name);
        if (it == m_type_registry.end() || it->second.load == nullptr) {
            return error_message(path + ": unknown type '" + type_name + "'");
        }
        return it->second.load(value, path, policy, id, env);
    }
};
}