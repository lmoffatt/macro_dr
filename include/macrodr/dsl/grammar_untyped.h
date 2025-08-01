#ifndef GRAMMAR_UNTYPED_H
#define GRAMMAR_UNTYPED_H

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "fold.h"
#include "grammar_Identifier.h"
#include "grammar_typed.h"
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
    virtual ~untyped_statement() {};
    virtual std::string str() const = 0;
    virtual untyped_statement* clone() const = 0;
    virtual std::string type_of_statement() const = 0;
    virtual Maybe_unique<base_typed_statement<Lexer, Compiler>> compile_statement(
        Compiler& cm) const = 0;
};

template <class Lexer, class Compiler>
class untyped_expression : public untyped_statement<Lexer, Compiler> {
   public:
    virtual ~untyped_expression() {};
    virtual untyped_expression* clone() const = 0;

    virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_expression(
        Compiler& cm) const = 0;

    virtual Maybe_unique<base_typed_statement<Lexer, Compiler>> compile_statement(
        Compiler& cm) const override {
        auto Maybe_exp = compile_expression(cm);
        if (!Maybe_exp)
            return Maybe_exp.error();
        else
            return Maybe_exp.value().release();
    }
};

template <class Lexer, class Compiler>
class untyped_program {
    std::string m_unformatted;
    std::vector<std::unique_ptr<untyped_statement<Lexer, Compiler>>> m_statements;

   public:
    untyped_program() {
    }
    untyped_program(const untyped_program& other)
        : m_unformatted{other.m_unformatted}, m_statements{clone(other.m_statements)} {
    }

    untyped_program(untyped_program&& other)
        : m_unformatted{std::move(other.m_unformatted)},
          m_statements{std::move(other.m_statements)} {
    }

    auto& push_back(untyped_statement<Lexer, Compiler>* t_expr) {
        m_statements.emplace_back(t_expr);
        return *this;
    }

    auto& push_back(const std::string& s) {
        m_unformatted += s + std::string(Lexer::statement_sep);
        return *this;
    }

    std::string str() const {
        std::string out = "";
        for (auto& elem : m_statements) out += elem->str() + std::string(Lexer::statement_sep);
        return out;
    };

    auto& unformatted() const {
        return m_unformatted;
    }
    friend std::ostream& operator<<(std::ostream& os, untyped_program const& p) {
        os << p.str();
        os << "\nunfomatted\n" << p.unformatted() << "\n";
        return os;
    }
    auto& statements() const {
        return m_statements;
    }
};

template <class Lexer, class Compiler>
class untyped_literal : public untyped_expression<Lexer, Compiler> {
    std::string m_expression;

   public:
    virtual ~untyped_literal() {};
    virtual untyped_literal* clone() const override {
        return new untyped_literal(*this);
    };

    virtual std::string type_of_statement() const {
        return statement_type();
    }
    static std::string statement_type() {
        return "untyped_simple_expression";
    }

    untyped_literal(std::string t_expression) : m_expression{std::move(t_expression)} {
    }

    virtual std::string str() const override {
        return m_expression;
    };
    auto& operator()() const {
        return m_expression;
    };

   public:
    virtual Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> compile_expression(
        Compiler&) const override {
        return error_message("this should not happen");
    }
};

template <class Lexer, class Compiler>
class untyped_numeric_literal : public untyped_literal<Lexer, Compiler> {
   public:
    using untyped_literal<Lexer, Compiler>::untyped_literal;
    virtual untyped_numeric_literal* clone() const override {
        return new untyped_numeric_literal(*this);
    };
    virtual ~untyped_numeric_literal() {};

    virtual std::string type_of_statement() const {
        return statement_type();
    }
    static std::string statement_type() {
        return "untyped_numeric_value";
    }

   public:
    virtual Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> compile_expression(
        Compiler&) const override {
        double out = std::stod(this->str());
        return new typed_literal<Lexer, Compiler, double>(out);
    }
};

template <class Lexer, class Compiler>
class untyped_string_literal : public untyped_literal<Lexer, Compiler> {
   public:
    using untyped_literal<Lexer, Compiler>::untyped_literal;

    virtual untyped_string_literal* clone() const override {
        return new untyped_string_literal(*this);
    };
    virtual ~untyped_string_literal() {};

    virtual Maybe_error<base_typed_statement<Lexer, Compiler>*> compile(const Compiler& cm) const {
        return error_message("");
    };
    virtual std::string type_of_statement() const {
        return statement_type();
    }
    static std::string statement_type() {
        return "untyped_string_value";
    }

   public:
    virtual Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> compile_expression(
        Compiler&) const override {
        return new typed_literal<Lexer, Compiler, std::string>(this->str());
    }
};

template <class Lexer, class Compiler>
class untyped_literal_comment : public untyped_literal<Lexer, Compiler> {
   public:
    using untyped_literal<Lexer, Compiler>::untyped_literal;

    virtual untyped_literal_comment* clone() const override {
        return new untyped_literal_comment(*this);
    };

    virtual Maybe_error<base_typed_statement<Lexer, Compiler>*> compile(const Compiler& cm) const {
        return error_message("");
    };

    virtual std::string str() const {
        return std::string(Lexer::comment_start) + untyped_literal<Lexer, Compiler>::str();
    }
    virtual ~untyped_literal_comment() {};
    virtual std::string type_of_statement() const {
        return statement_type();
    }
    static std::string statement_type() {
        return "untyped_comment";
    }
};

template <class Lexer, class Compiler>
class untyped_identifier : public untyped_expression<Lexer, Compiler> {
    Identifier<Lexer> m_id;

   public:
    virtual ~untyped_identifier() {};
    untyped_identifier(Identifier<Lexer> t_id) : m_id{t_id} {
    }
    virtual untyped_identifier* clone() const override {
        return new untyped_identifier(*this);
    };

    virtual std::string str() const override final {
        return m_id();
    }

    auto& operator()() const {
        return m_id;
    }

    // untyped_statement interface
   public:
    virtual std::string type_of_statement() const {
        return statement_type();
    }
    static std::string statement_type() {
        return "untyped_identifier";
    }

   public:
    virtual Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> compile_expression(
        Compiler& cm) const override {
        return cm.get_Identifier(m_id);
    }
};

template <class Lexer, class Compiler>
class untyped_assignment : public untyped_statement<Lexer, Compiler> {
    std::unique_ptr<untyped_identifier<Lexer, Compiler>> m_id;
    std::unique_ptr<untyped_expression<Lexer, Compiler>> m_expr;

   public:
    virtual ~untyped_assignment() override {
    }
    untyped_assignment(untyped_identifier<Lexer, Compiler>* t_id,
                       untyped_expression<Lexer, Compiler>* t_expr)
        : m_id{t_id}, m_expr{t_expr} {
    }

    untyped_assignment() {
    }
    untyped_assignment(untyped_assignment&& other)
        : m_id{std::move(other.m_id)}, m_expr{std::move(other.m_expr)} {
    }
    untyped_assignment(untyped_assignment const& other)
        : m_id{other.m_id->clone()}, m_expr{other.m_expr->clone()} {
    }

    auto& id() const {
        return *m_id;
    }
    auto& expression() const {
        return *m_expr;
    }

    virtual std::string str() const override final {
        return m_id->str() + std::string(Lexer::assignment_operator) +
               std::string(Lexer::assignment_sep) + m_expr->str();
    }
    virtual untyped_assignment* clone() const override {
        return new untyped_assignment(*this);
    };

    virtual std::string type_of_statement() const {
        return statement_type();
    }
    static std::string statement_type() {
        return "untyped_assignment";
    }

   public:
    virtual Maybe_unique<dsl::base_typed_statement<Lexer, Compiler>> compile_statement(
        Compiler& cm) const override {
        auto Maybe_typed_expr = expression().compile_expression(cm);
        if (!Maybe_typed_expr)
            return Maybe_typed_expr.error();
        else {
            auto cid = Maybe_typed_expr.value()->compile_identifier();
            cm.push_back(id()(), cid);
            return std::unique_ptr<base_typed_statement<Lexer, Compiler>>(
                Maybe_typed_expr.value()->compile_assigment(id()()));
        }
    }
};

template <class Lexer, class Compiler>
class untyped_argument_list : public untyped_expression<Lexer, Compiler> {
    std::vector<std::unique_ptr<untyped_statement<Lexer, Compiler>>> m_list;

   public:
    untyped_argument_list(
        const std::vector<std::unique_ptr<untyped_statement<Lexer, Compiler>>>& t_list)
        : m_list{dsl::clone(t_list)} {
    }

    untyped_argument_list(std::vector<std::unique_ptr<untyped_statement<Lexer, Compiler>>>&& t_list)
        : m_list{t_list} {
    }
    untyped_argument_list(const untyped_argument_list& other)
        : m_list{dsl::clone_vector(other.m_list)} {
    }

    untyped_argument_list(untyped_argument_list&& other) : m_list{std::move(other.t_list)} {
    }

    untyped_argument_list() {
    }

    auto& arg() const {
        return m_list;
    }

    untyped_argument_list& push_back(untyped_statement<Lexer, Compiler>* t_expr) {
        m_list.emplace_back(t_expr);
        return *this;
    }
    virtual untyped_argument_list* clone() const override {
        return new untyped_argument_list(*this);
    }

    virtual std::string type_of_statement() const {
        return statement_type();
    }
    static std::string statement_type() {
        return "untyped_argument_list";
    }
    virtual std::string str() const override final {
        std::string out = std::string(Lexer::argument_list_start);
        if (m_list.size() > 0) {
            out += m_list[0]->str();
            for (auto i = 1ul; i < m_list.size(); ++i) {
                out += std::string(Lexer::argument_sep) + m_list[i]->str();
            }
        }
        out += Lexer::argument_list_end;
        return out;
    }

   public:
    virtual Maybe_unique<dsl::base_typed_expression<Lexer, Compiler>> compile_expression(
        Compiler& cm) const override {
        auto out = compile_argument_list(cm);
        if (!out)
            return out.error();
        else
            return out.value();
    }

    virtual Maybe_error<typed_argument_list<Lexer, Compiler>*> compile_argument_list(
        Compiler& cm) const {
        auto out = new typed_argument_list<Lexer, Compiler>{};
        for (auto& e : arg()) {
            auto v_c = e->compile_statement(cm);
            if (!v_c)
                return v_c.error();
            else {
                auto May_out = v_c.value()->compile_argument_list(out);
                if (!May_out)
                    return May_out.error();
                else
                    out = May_out.value();
            }
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
        : m_fid{t_fid}, m_arg{t_arg} {
    }
    virtual ~untyped_function_evaluation() {};
    untyped_function_evaluation(untyped_function_evaluation&& other)
        : m_fid{std::move(other.m_fid)}, m_arg{std::move(other.m_arg)} {
    }

    untyped_function_evaluation(untyped_function_evaluation const& other)
        : m_fid{other.m_fid->clone()}, m_arg{other.m_arg->clone()} {
    }

    auto& fid() const {
        return *m_fid;
    }
    auto& args() const {
        return *m_arg;
    }

    virtual std::string type_of_statement() const override {
        return statement_type();
    }
    static std::string statement_type() {
        return "untyped_function_evaluation";
    }

    virtual untyped_function_evaluation* clone() const override {
        return new untyped_function_evaluation(*this);
    }
    virtual std::string str() const override final {
        return m_fid->str() + m_arg->str();
    }

    // untyped_expression interface
   public:
    virtual Maybe_unique<base_typed_expression<Lexer, Compiler>> compile_expression(
        Compiler& cm) const override {
        auto Maybe_fn = cm.get_function(fid()());
        if (!Maybe_fn)
            return Maybe_fn.error();
        else {
            auto cmc = cm;
            //   auto maybe_args = args().com<pile_argument_list(cmc);
            // if (!maybe_args)
            //   return maybe_args.error();
            // else
            // {
            return Maybe_fn.value()->compile_function_evaluation(cmc, args());
            // }
        }
    }
};

template <class Lexer, class Compiler>
inline Maybe_error<typed_program<Lexer, Compiler>> compile_program(
    Compiler& cm, const untyped_program<Lexer, Compiler>& s) {
    typed_program<Lexer, Compiler> out;
    for (auto& e : s.statements()) {
        auto Maybe_compiled_statement = e->compile_statement(cm);
        if (!Maybe_compiled_statement)
            return Maybe_compiled_statement.error();
        else {
            out.push_back(Maybe_compiled_statement.value().release());
        }
    }
    return out;
}

inline Maybe_error<typed_program<Lexer, Compiler>> compile_and_run_program_line_by_line(
    Compiler& cm, const untyped_program<Lexer, Compiler>& s) {
    typed_program<Lexer, Compiler> out;
    for (auto& e : s.statements()) {
        auto Maybe_compiled_statement = e->compile_statement(cm);
        if (!Maybe_compiled_statement)
            return Maybe_compiled_statement.error();
        else
            out.push_back(Maybe_compiled_statement.value().release());
    }
    return out;
}


}  // namespace macrodr::dsl

#endif  // GRAMMAR_UNTYPED_H
