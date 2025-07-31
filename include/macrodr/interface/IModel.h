#pragma once

#include "IObject.h"
#include "parameters.h"
#include "qmodel.h"

namespace macrodr {

namespace interface {

template <typename ParamValues>
struct IParameterApplication {
    using Result = macrodr::Transfer_Op_to<ParamValues, macrodr::Patch_Model>;  // or a better-named

    virtual Maybe_error<Result> operator()(const ParamValues& p) const = 0;
    virtual ~IParameterApplication() = default;
};

template <typename... ParamValues>
struct IModel : public IObject, public IParameterApplication<ParamValues>... {
    // --- Model metadata ---
    virtual std::string model_name() const = 0;
    virtual const std::vector<std::string>& parameter_names() const = 0;
    virtual const var::Parameters_Transformations& get_transformations() const = 0;

    // --- Internal structure ---
    virtual const macrodr::Q0_formula& get_Q0_formula() const = 0;
    virtual const macrodr::Qa_formula& get_Qa_formula() const = 0;
    virtual const macrodr::g_formula& get_g_formula() const = 0;
    virtual std::size_t number_of_states() const = 0;

    // --- Optional introspection / postcondition hooks ---
    //virtual bool check_invariants(const ParamValues& p) const { return true; }
    virtual std::unique_ptr<IModel<ParamValues...>> clone() const = 0;

    virtual ~IModel() = default;
};

template <typename Scheme>
struct ConcreteBaseScheme {
   protected:
    Scheme scheme;
    template <class S>
    explicit ConcreteBaseScheme(S&& s) : scheme(std::forward<S>(s)) {
    }

   public:
    // Forward metadata helpers so derived classes can just use them.
    const auto& _scheme() const {
        return scheme;
    }
};

template <typename Scheme, typename ParamValues>
struct ConcreteParameterApplication : public virtual ConcreteBaseScheme<Scheme>,
                                      public IParameterApplication<ParamValues> {
    using Result = macrodr::Transfer_Op_to<ParamValues, macrodr::Patch_Model>;
    using ConcreteBaseScheme<Scheme>::ConcreteBaseScheme;  // inherit ctor

    Maybe_error<Result> operator()(const ParamValues& p) const override {
        return ConcreteBaseScheme<Scheme>::scheme(p);
    }
};

template <typename Scheme, typename... ParamValues>
class ConcreteModel : public IModel<ParamValues...>,
                      public ConcreteParameterApplication<Scheme, ParamValues>... {
    using Base = ConcreteBaseScheme<Scheme>;

   public:
    using ConcreteParameterApplication<Scheme, ParamValues>::operator()...;

    template <class S>
    explicit ConcreteModel(S&& scheme)
        : Base(std::forward<S>(scheme)),
          ConcreteParameterApplication<Scheme, ParamValues>(Base::scheme)... {
    }

    // --- metadata passthrough ---
    std::string model_name() const override {
        return Base::scheme.model_name();
    }
    const std::vector<std::string>& parameter_names() const override {
        return Base::scheme.parameter_names();
    }
    const var::Parameters_Transformations& get_transformations() const override {
        return Base::scheme.get_transformations();
    }
    const macrodr::Q0_formula& get_Q0_formula() const override {
        return Base::scheme.get_Q0_formula();
    }
    const macrodr::Qa_formula& get_Qa_formula() const override {
        return Base::scheme.get_Qa_formula();
    }
    const macrodr::g_formula& get_g_formula() const override {
        return Base::scheme.get_g_formula();
    }
    std::size_t number_of_states() const override {
        return Base::scheme.number_of_states();
    }

    std::unique_ptr<IModel<ParamValues...>> clone() const override {
        return std::make_unique<ConcreteModel>(*this);
    }
};
}  // namespace interface
}  // namespace macrodr
