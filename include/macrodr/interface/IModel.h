#pragma once

#include <memory>

#include "IObject.h"
#include "parameters.h"
#include "qmodel.h"

namespace macrodr::interface {

template <typename ParamValues>
struct IParameterApplication {
    IParameterApplication(const IParameterApplication&) = delete;
    IParameterApplication(IParameterApplication&&) = default;
    IParameterApplication& operator=(const IParameterApplication&) = delete;
    IParameterApplication& operator=(IParameterApplication&&) = default;
    using Result = macrodr::Transfer_Op_to<ParamValues, macrodr::Patch_Model>;

    [[nodiscard]] virtual Maybe_error<Result> operator()(const ParamValues& parameters) const = 0;

    [[nodiscard]] virtual Maybe_error<Result> operator()(const Result& t_p) const = 0;
    virtual ~IParameterApplication() = default;
    virtual void report_model(save_Parameter<ParamValues>& paramters_file) = 0;
};

template <typename... ParamValues>
struct IModel : public IObject, public IParameterApplication<ParamValues>... {
    IModel(const IModel&) = delete;
    IModel(IModel&&) = default;
    IModel& operator=(const IModel&) = delete;
    IModel& operator=(IModel&&) = default;
    // --- Model metadata ---
    [[nodiscard]] virtual std::string model_name() const = 0;
    [[nodiscard]] virtual std::vector<std::string> const& names() const = 0;

    [[nodiscard]] virtual var::Parameters_Transformations const& parameters_transformations()
        const = 0;

    [[nodiscard]] virtual const macrodr::Q0_formula& get_Q0_formula() const = 0;
    [[nodiscard]] virtual const macrodr::Qa_formula& get_Qa_formula() const = 0;
    [[nodiscard]] virtual const macrodr::g_formula& get_g_formula() const = 0;
    [[nodiscard]] virtual std::size_t number_of_states() const = 0;

    // --- Optional introspection / postcondition hooks ---
    // virtual bool check_invariants(const ParamValues& p) const { return true; }
    virtual std::unique_ptr<IModel<ParamValues...>> clone() const = 0;

    virtual ~IModel() = default;
};

template <typename Scheme>
struct ConcreteBaseScheme {
   private:
    Scheme scheme_;  // private: no tidy warning

   protected:
    // read-only view for derived classes
    const Scheme& scheme() const { return scheme_; }
    // (optional) mutable access if you truly need it
    //  Scheme& scheme() { return scheme_; }

   public:
    explicit ConcreteBaseScheme(Scheme sch) noexcept(std::is_nothrow_move_constructible_v<Scheme>)
        : scheme_(std::move(sch)) {}
};

template <typename Scheme, typename ParamValues>
struct ConcreteParameterApplication : public ConcreteBaseScheme<Scheme>,
                                      public IParameterApplication<ParamValues> {
    using Result = macrodr::Transfer_Op_to<ParamValues, macrodr::Patch_Model>;

    using ConcreteBaseScheme<Scheme>::ConcreteBaseScheme;  // inherit ctor

    Maybe_error<Result> operator()(const ParamValues& parameter) const override {
        return this->scheme()(parameter);  // protected accessor
    }
};

template <typename Scheme, typename... ParamValues>
class ConcreteModel : public IModel<ParamValues...>,
                      public ConcreteParameterApplication<Scheme, ParamValues>... {
    using Base = ConcreteBaseScheme<Scheme>;

   public:
    using ConcreteParameterApplication<Scheme, ParamValues>::operator()...;

    explicit ConcreteModel(Scheme scheme)
        : Base(std::move(scheme)),
          ConcreteParameterApplication<Scheme, ParamValues>(Base::scheme)... {}

    // --- metadata passthrough ---
    [[nodiscard]] std::string model_name() const override { return Base::scheme().model_name(); }
    [[nodiscard]] const std::vector<std::string>& names() const override {
        return Base::scheme().parameter_names();
    }
    [[nodiscard]] const var::Parameters_Transformations& parameters_transformations()
        const override {
        return Base::scheme().get_transformations();
    }
    [[nodiscard]] const macrodr::Q0_formula& get_Q0_formula() const override {
        return Base::scheme().get_Q0_formula();
    }
    [[nodiscard]] const macrodr::Qa_formula& get_Qa_formula() const override {
        return Base::scheme().get_Qa_formula();
    }
    [[nodiscard]] const macrodr::g_formula& get_g_formula() const override {
        return Base::scheme().get_g_formula();
    }
    [[nodiscard]] std::size_t number_of_states() const override {
        return Base::scheme().number_of_states();
    }

    [[nodiscard]] std::unique_ptr<IModel<ParamValues...>> clone() const override {
        return std::make_unique<ConcreteModel>(*this);
    }
};

template <typename Scheme, typename... ParamValues>
auto make_model_interface(Scheme scheme) {
    return std::make_unique<ConcreteModel<Scheme, ParamValues...>>(std::move(scheme));
}
}  // namespace macrodr::interface
