#ifndef GRAMMAR_UNTYPED_H
#define GRAMMAR_UNTYPED_H

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "fold.h"
#include "grammar_Identifier.h"
#include "grammar_typed.h"
#include "lexer_typed.h"
#include "maybe_error.h"
namespace macrodr::dsl {

template <class Abstract>
auto clone(const std::unique_ptr<Abstract>& x) {
    return std::unique_ptr<Abstract>(x->clone());
}
template <class Abstract>
std::vector<std::unique_ptr<Abstract>> clone(const std::vector<std::unique_ptr<Abstract>>& x) {
    return Map(x, [](auto& e) { return std::unique_ptr<Abstract>(e->clone()); });
}

template <class... Abstract>
std::tuple<std::unique_ptr<Abstract>...> clone_tuple(
    const std::tuple<std::unique_ptr<Abstract>...>& x) {
    return std::apply(
        [](auto const&... e) -> std::tuple<std::unique_ptr<Abstract>...> {
            return std::tuple<std::unique_ptr<Abstract>...>(
                std::unique_ptr<Abstract>(e->clone())...);
        },
        x);
}

template <class Lexer, class Compiler>
class base_typed_statement;
template <class Lexer, class Compiler>
class base_typed_expression;

template <class Lexer, class Compiler>
class untyped_statement {
   public:
    virtual ~untyped_statement() = default;
    [[nodiscard]] virtual std::string str() const = 0;
    [[nodiscard]] virtual untyped_statement* clone() const = 0;
    [[nodiscard]] virtual std::string type_of_statement() const = 0;
    virtual Maybe_unique<base_typed_statement<Lexer, Compiler>> compile_statement(
        Environment<Lexer, Compiler>& cm) const = 0;
};

template <class Lexer, class Compiler>
class untyped_expression : public untyped_statement<Lexer, Compiler> {
   public:
    ~untyped_expression() override = default;
    [[nodiscard]] untyped_expression* clone() const override = 0;

    virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_expression(
        Environment<Lexer, Compiler>& cm) const = 0;

    Maybe_unique<base_typed_statement<Lexer, Compiler>> compile_statement(
        Environment<Lexer, Compiler>& cm) const override {
        auto Maybe_exp = compile_expression(cm);
        if (!Maybe_exp) {
            return Maybe_exp.error();
        }
        return Maybe_exp.value().release();
    }
};

template <class Lexer, class Compiler>
class untyped_program {
    std::string m_unformatted;
    std::vector<std::unique_ptr<untyped_statement<Lexer, Compiler>>> m_statements;

   public:
    untyped_program() = default;
    untyped_program(const untyped_program& other)
        : m_unformatted{other.m_unformatted}, m_statements{clone(other.m_statements)} {}

    untyped_program(untyped_program&& other) noexcept
        : m_unformatted{std::move(other.m_unformatted)},
          m_statements{std::move(other.m_statements)} {}

    auto& push_back(untyped_statement<Lexer, Compiler>* t_expr) {
        m_statements.emplace_back(t_expr);
        return *this;
    }

    auto& push_back(const std::string& s) {
        m_unformatted += s + std::string(Lexer::statement_sep);
        return *this;
    }

    [[nodiscard]] std::string str() const {
        std::string out;
        for (auto& elem : m_statements) {
            out += elem->str() + std::string(Lexer::statement_sep);
        }
        return out;
    };

    [[nodiscard]] auto& unformatted() const { return m_unformatted; }
    friend std::ostream& operator<<(std::ostream& os, untyped_program const& p) {
        os << p.str();
        os << "\n [unfomatted]\n" << p.unformatted() << "\n";
        return os;
    }
    [[nodiscard]] auto& statements() const { return m_statements; }
};

template <class Lexer, class Compiler>
class untyped_literal : public untyped_expression<Lexer, Compiler> {
    std::string m_expression;

   public:
    ~untyped_literal() override = default;
    [[nodiscard]] untyped_literal* clone() const override { return new untyped_literal(*this); };

    [[nodiscard]] std::string type_of_statement() const override { return statement_type(); }
    static std::string statement_type() { return "untyped_simple_expression"; }

    untyped_literal(std::string t_expression) : m_expression{std::move(t_expression)} {}

    [[nodiscard]] std::string str() const override { return m_expression; };
    auto& operator()() const { return m_expression; };

    Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> compile_expression(
        Environment<Lexer, Compiler>& /*unused*/) const override {
        return error_message("this should not happen");
    }
};

template <class Lexer, class Compiler>
class untyped_numeric_literal : public untyped_literal<Lexer, Compiler> {
   public:
    using untyped_literal<Lexer, Compiler>::untyped_literal;
    [[nodiscard]] untyped_numeric_literal* clone() const override {
        return new untyped_numeric_literal(*this);
    };
    ~untyped_numeric_literal() override = default;

    [[nodiscard]] std::string type_of_statement() const override { return statement_type(); }
    static std::string statement_type() { return "untyped_numeric_value"; }

    Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> compile_expression(
        Environment<Lexer, Compiler>& /*unused*/) const override {
        double const out = std::stod(this->str());
        return new typed_literal<Lexer, Compiler, double>(out);
    }
};

template <class Lexer, class Compiler>
class untyped_string_literal : public untyped_literal<Lexer, Compiler> {
   public:
    using untyped_literal<Lexer, Compiler>::untyped_literal;

    [[nodiscard]] untyped_string_literal* clone() const override {
        return new untyped_string_literal(*this);
    };
    ~untyped_string_literal() override = default;

    [[nodiscard]] virtual Maybe_error<base_typed_statement<Lexer, Compiler>*> compile(
        const Environment<Lexer, Compiler>& /*cm*/) const {
        return error_message("");
    };
    [[nodiscard]] std::string type_of_statement() const override { return statement_type(); }
    static std::string statement_type() { return "untyped_string_value"; }

    Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> compile_expression(
        Environment<Lexer, Compiler>& /*unused*/) const override {
        return new typed_literal<Lexer, Compiler, std::string>(this->str());
    }
};

template <class Lexer, class Compiler>
class untyped_literal_comment : public untyped_literal<Lexer, Compiler> {
   public:
    using untyped_literal<Lexer, Compiler>::untyped_literal;

    [[nodiscard]] untyped_literal_comment* clone() const override {
        return new untyped_literal_comment(*this);
    };

    [[nodiscard]] virtual Maybe_error<base_typed_statement<Lexer, Compiler>*> compile(
        const Environment<Lexer, Compiler>& /*cm*/) const {
        return error_message("");
    };

    [[nodiscard]] std::string str() const override {
        return std::string(Lexer::comment_start) + untyped_literal<Lexer, Compiler>::str();
    }
    ~untyped_literal_comment() override = default;
    [[nodiscard]] std::string type_of_statement() const override { return statement_type(); }
    static std::string statement_type() { return "untyped_comment"; }
};

template <class Lexer, class Compiler>
class untyped_identifier : public untyped_expression<Lexer, Compiler> {
    Identifier<Lexer> m_id;

   public:
    ~untyped_identifier() override = default;
    untyped_identifier(Identifier<Lexer> t_id) : m_id{std::move(std::move(t_id))} {}
    [[nodiscard]] untyped_identifier* clone() const override {
        return new untyped_identifier(*this);
    };

    [[nodiscard]] std::string str() const final { return m_id(); }

    auto& operator()() const { return m_id; }

    // untyped_statement interface

    [[nodiscard]] std::string type_of_statement() const override { return statement_type(); }
    static std::string statement_type() { return "untyped_identifier"; }

    Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> compile_expression(
        Environment<Lexer, Compiler>& cm) const override {
        return cm.get_Identifier(m_id);
    }
};

template <class Lexer, class Compiler>
class untyped_assignment : public untyped_statement<Lexer, Compiler> {
    std::unique_ptr<untyped_identifier<Lexer, Compiler>> m_id;
    std::unique_ptr<untyped_expression<Lexer, Compiler>> m_expr;

   public:
    ~untyped_assignment() override = default;
    untyped_assignment(untyped_identifier<Lexer, Compiler>* t_id,
                       untyped_expression<Lexer, Compiler>* t_expr)
        : m_id{t_id}, m_expr{t_expr} {}

    untyped_assignment() = default;
    untyped_assignment(untyped_assignment&& other) noexcept
        : m_id{std::move(other.m_id)}, m_expr{std::move(other.m_expr)} {}
    untyped_assignment(untyped_assignment const& other)
        : m_id{other.m_id->clone()}, m_expr{other.m_expr->clone()} {}

    [[nodiscard]] auto& id() const { return *m_id; }
    [[nodiscard]] auto& expression() const { return *m_expr; }

    [[nodiscard]] std::string str() const final {
        return m_id->str() + std::string(Lexer::assignment_operator) +
               std::string(Lexer::assignment_sep) + m_expr->str();
    }
    [[nodiscard]] untyped_assignment* clone() const override {
        return new untyped_assignment(*this);
    };

    [[nodiscard]] std::string type_of_statement() const override { return statement_type(); }
    static std::string statement_type() { return "untyped_assignment"; }

    Maybe_unique<dsl::base_typed_statement<Lexer, Compiler>> compile_statement(
        Environment<Lexer, Compiler>& cm) const override {
        auto Maybe_typed_expr = expression().compile_expression(cm);
        if (!Maybe_typed_expr) {
            return Maybe_typed_expr.error();
        }
        auto& typed_expr = Maybe_typed_expr.value();
        cm.push_back(id()(), typed_expr->compile_identifier_unique());
        return typed_expr->compile_assigment_unique(id()());
    }
};

template <class Lexer, class Compiler>
class untyped_argument_list : public untyped_expression<Lexer, Compiler> {
    std::vector<std::unique_ptr<untyped_statement<Lexer, Compiler>>> m_list;

   public:
    untyped_argument_list(
        const std::vector<std::unique_ptr<untyped_statement<Lexer, Compiler>>>& t_list)
        : m_list{dsl::clone(t_list)} {}

    untyped_argument_list(std::vector<std::unique_ptr<untyped_statement<Lexer, Compiler>>>&& t_list)
        : m_list{std::move(t_list)} {}
    untyped_argument_list(const untyped_argument_list& other)
        : m_list{dsl::clone_vector(other.m_list)} {}

    untyped_argument_list(untyped_argument_list&& other) noexcept
        : m_list{std::move(other.m_list)} {}

    untyped_argument_list() = default;

    [[nodiscard]] auto& arg() const { return m_list; }

    untyped_argument_list& push_back(untyped_statement<Lexer, Compiler>* t_expr) {
        m_list.emplace_back(t_expr);
        return *this;
    }
    [[nodiscard]] untyped_argument_list* clone() const override {
        return new untyped_argument_list(*this);
    }

    [[nodiscard]] std::string type_of_statement() const override { return statement_type(); }
    static std::string statement_type() { return "untyped_argument_list"; }
    [[nodiscard]] std::string str() const final {
        std::string out = std::string(Lexer::argument_list_start);
        if (m_list.size() > 0) {
            out += m_list[0]->str();
            for (auto i = 1UL; i < m_list.size(); ++i) {
                out += std::string(Lexer::argument_sep) + m_list[i]->str();
            }
        }
        out += Lexer::argument_list_end;
        return out;
    }

    Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> compile_expression(
        Environment<Lexer, Compiler>& cm) const override {
        auto out = compile_argument_list(cm);
        if (!out) {
            return out.error();
        }
        auto list = std::move(out.value());
        return std::unique_ptr<base_typed_expression<Lexer, Compiler>>(std::move(list));
    }

    virtual Maybe_unique<typed_argument_list<Lexer, Compiler>> compile_argument_list(
        Environment<Lexer, Compiler>& cm) const {
        auto out = std::make_unique<typed_argument_list<Lexer, Compiler>>();
        auto* current = out.get();
        for (auto& e : arg()) {
            auto v_c = e->compile_statement(cm);
            if (!v_c) {
                return v_c.error();
            }
            auto May_out = v_c.value()->compile_argument_list(current);
            if (!May_out) {
                return May_out.error();
            }
            current = May_out.value();
        }
        return out;
    }
};

template <class Lexer, class Compiler>
class untyped_function_evaluation : public untyped_expression<Lexer, Compiler> {
    std::unique_ptr<untyped_identifier<Lexer, Compiler>> m_fid;
    std::unique_ptr<untyped_argument_list<Lexer, Compiler>> m_arg;

   public:
    untyped_function_evaluation(untyped_identifier<Lexer, Compiler>* t_fid,
                                untyped_argument_list<Lexer, Compiler>* t_arg)
        : m_fid{t_fid}, m_arg{t_arg} {}
    ~untyped_function_evaluation() override = default;
    untyped_function_evaluation(untyped_function_evaluation&& other) noexcept
        : m_fid{std::move(other.m_fid)}, m_arg{std::move(other.m_arg)} {}

    untyped_function_evaluation(untyped_function_evaluation const& other)
        : m_fid{other.m_fid->clone()}, m_arg{other.m_arg->clone()} {}

    [[nodiscard]] auto& fid() const { return *m_fid; }
    [[nodiscard]] auto& args() const { return *m_arg; }

    [[nodiscard]] std::string type_of_statement() const override { return statement_type(); }
    static std::string statement_type() { return "untyped_function_evaluation"; }

    [[nodiscard]] untyped_function_evaluation* clone() const override {
        return new untyped_function_evaluation(*this);
    }
    [[nodiscard]] std::string str() const final { return m_fid->str() + m_arg->str(); }

    // untyped_expression interface

    Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_expression(
        Environment<Lexer, Compiler>& cm) const override {
        // Try all overloads for this function identifier in registration order
        auto Maybe_fns = cm.compiler().get_functions(fid()());
        if (!Maybe_fns) {
            return Maybe_fns.error();
        }
        const auto& overloads = Maybe_fns.value();
        std::string errors;
        std::size_t i_overload = 0;
        for (auto* fn : overloads) {
            if (fn == nullptr) {
                continue;
            }
            auto attempt = fn->compile_function_evaluation(cm, args());
            if (attempt) {
                return attempt;  // success on this overload
            }
            // Accumulate diagnostics (best-effort, do not explode)

            errors += std::to_string(i_overload) + " overload: " + attempt.error()();
            ++i_overload;
            errors += "\n\n";
        }
        return error_message(std::string{"\nno matching overload for function '"} + fid().str() +
                             "' with provided arguments:"+args().str()+"\n" + errors);
    }
};  

template <class Lexer, class Compiler>
class untyped_vector_construction : public untyped_expression<Lexer, Compiler> {
    std::unique_ptr<untyped_argument_list<Lexer, Compiler>> m_arg;

   public:
    untyped_vector_construction(untyped_argument_list<Lexer, Compiler>* t_arg)
        :  m_arg{t_arg} {}
    ~untyped_vector_construction() override = default;
    untyped_vector_construction(untyped_vector_construction&& other) noexcept
        :  m_arg{std::move(other.m_arg)} {}

    untyped_vector_construction(untyped_vector_construction const& other)
        :  m_arg{other.m_arg->clone()} {}

    [[nodiscard]] auto& args() const { return *m_arg; }

    [[nodiscard]] std::string type_of_statement() const override { return statement_type(); }
    static std::string statement_type() { return "untyped_vector_construction"; }

    [[nodiscard]] untyped_vector_construction* clone() const override {
        return new untyped_vector_construction(*this);
    }
    [[nodiscard]] std::string str() const final { return std::string(Lexer::vector_list_start) + m_arg->str()+ std::string(Lexer::vector_list_end) ; }

    // untyped_expression interface

    Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_expression(
        Environment<Lexer, Compiler>& /*cm*/) const override {
        return error_message("type unknown at this point");
    }
};

template <class Lexer, class Compiler>
class untyped_tuple_construction : public untyped_expression<Lexer, Compiler> {
    std::unique_ptr<untyped_argument_list<Lexer, Compiler>> m_arg;

   public:
    untyped_tuple_construction(untyped_argument_list<Lexer, Compiler>* t_arg)
        :  m_arg{t_arg} {}
    ~untyped_tuple_construction() override = default;
    untyped_tuple_construction(untyped_tuple_construction&& other) noexcept
        :  m_arg{std::move(other.m_arg)} {}

    untyped_tuple_construction(untyped_tuple_construction const& other)
        :  m_arg{other.m_arg->clone()} {}

    [[nodiscard]] auto& args() const { return *m_arg; }

    [[nodiscard]] std::string type_of_statement() const override { return statement_type(); }
    static std::string statement_type() { return "untyped_tuple_construction"; }

    [[nodiscard]] untyped_tuple_construction* clone() const override {
        return new untyped_tuple_construction(*this);
    }
    [[nodiscard]] std::string str() const final { return std::string(Lexer::tuple_list_start) + m_arg->str()+ 
        std::string(Lexer::tuple_list_end) ; }

    // untyped_expression interface

    Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_expression(
        Environment<Lexer, Compiler>& /*cm*/) const override {
        return error_message("type unknown at this point");
    }
};


template <class Lexer, class Compiler>
inline Maybe_error<typed_program<Lexer, Compiler>> compile_program(
    Environment<Lexer, Compiler>& cm, const untyped_program<Lexer, Compiler>& s) {
    typed_program<Lexer, Compiler> out;
    for (auto& e : s.statements()) {
        auto Maybe_compiled_statement = e->compile_statement(cm);
        if (!Maybe_compiled_statement) {
            return error_message(std::string("in ") + e->str() +
                                 "\n error:" + Maybe_compiled_statement.error()());
        }

        out.push_back(std::move(Maybe_compiled_statement.value()));
    }
    return out;
}

template <class Lexer, class Compiler>
inline Maybe_error<typed_program<Lexer, Compiler>> compile_and_run_program_line_by_line(
    Environment<Lexer, Compiler>& cm, const untyped_program<Lexer, Compiler>& s) {
    typed_program<Lexer, Compiler> out;
    for (const auto& e : s.statements()) {
        auto Maybe_compiled_statement = e->compile_statement(cm);
        if (!Maybe_compiled_statement) {
            return Maybe_compiled_statement.error();
        }
        out.push_back(std::move(Maybe_compiled_statement.value()));
    }
    return out;
}

}  // namespace macrodr::dsl

#endif  // GRAMMAR_UNTYPED_H
