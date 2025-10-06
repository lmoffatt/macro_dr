#pragma once

#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "IObject.h"
#include "parameters.h"
#include "qmodel.h"

namespace macrodr::interface {

// ─────────────────────────────────────────────────────────────────────────────
// 1) Per-parameter interface (pure virtual)
//    Result is defined as requested.
// ─────────────────────────────────────────────────────────────────────────────
template <typename ParamValues, typename PatchModel>
struct IParameterApplication {
    using Result = PatchModel;

    IParameterApplication() = default;
    IParameterApplication(const IParameterApplication&) = default;
    IParameterApplication(IParameterApplication&&) = default;
    IParameterApplication& operator=(const IParameterApplication&) = default;
    IParameterApplication& operator=(IParameterApplication&&) = default;
    virtual ~IParameterApplication() = default;

    [[nodiscard]] virtual Maybe_error<PatchModel> operator()(
        const ParamValues& parameters) const = 0;
};

// ─────────────────────────────────────────────────────────────────────────────
// 2) Polymorphic model aggregating multiple IParameterApplication<T>.
// ─────────────────────────────────────────────────────────────────────────────
template <typename... ParamValues>
struct IModel : public IObject,
                public virtual IParameterApplication<
                    ParamValues, macrodr::Transfer_Op_to<ParamValues, macrodr::Patch_Model>>... {
    // inside template <typename... ParamValues> struct IModel ...
    IModel(const IModel&) = default;
    IModel& operator=(const IModel&) = default;
    virtual ~IModel() = default;  // if not already present
    IModel() = default;
    IModel(IModel&&) = default;
    IModel& operator=(IModel&&) = default;

    [[nodiscard]] virtual std::string model_name() const = 0;
    [[nodiscard]] virtual const std::vector<std::string>& names() const = 0;
    [[nodiscard]] virtual std::size_t number_of_states() const = 0;
    template <typename PV>
        requires((std::is_same_v<PV, ParamValues> || ...))
    [[nodiscard]] Maybe_error<macrodr::Transfer_Op_to<PV, macrodr::Patch_Model>> operator()(
        const PV& parameters) const {
        using App = IParameterApplication<PV, macrodr::Transfer_Op_to<PV, macrodr::Patch_Model>>;
        return static_cast<const App&>(*this).operator()(parameters);
    }

    [[nodiscard]] virtual std::unique_ptr<IModel<ParamValues...>> clone() const = 0;
};

// ─────────────────────────────────────────────────────────────────────────────
// 3) Single shared base that stores Scheme and exposes scheme().
//    Virtual base => avoids multiple subobjects in the diamond.
// ─────────────────────────────────────────────────────────────────────────────
template <typename Scheme>
class ConcreteBaseScheme {
   public:
    explicit ConcreteBaseScheme(Scheme s) : scheme_(std::move(s)) {}
    ConcreteBaseScheme(const ConcreteBaseScheme&) = default;
    ConcreteBaseScheme(ConcreteBaseScheme&&) = default;
    ConcreteBaseScheme& operator=(const ConcreteBaseScheme&) = default;
    ConcreteBaseScheme& operator=(ConcreteBaseScheme&&) = default;
    virtual ~ConcreteBaseScheme() = default;

   protected:
    [[nodiscard]] const Scheme& scheme() const { return scheme_; }
    [[nodiscard]] Scheme& scheme() { return scheme_; }

   private:
    Scheme scheme_;
};

// ─────────────────────────────────────────────────────────────────────────────
// 4) Implements IParameterApplication<T> by forwarding to Scheme.
//    Virtual-inherit ConcreteBaseScheme so all specializations share one.
// ─────────────────────────────────────────────────────────────────────────────
template <typename Scheme, typename ParamValues, typename PatchModel>
class ConcreteParameterApplication : public virtual ConcreteBaseScheme<Scheme>,
                                     public virtual IParameterApplication<ParamValues, PatchModel> {
    using SchemeBase = ConcreteBaseScheme<Scheme>;
    using Result = PatchModel;

   public:
    // after
    ConcreteParameterApplication() noexcept {}  // user-provided, empty

    ConcreteParameterApplication(const ConcreteParameterApplication&) = default;
    ConcreteParameterApplication(ConcreteParameterApplication&&) = default;
    ConcreteParameterApplication& operator=(const ConcreteParameterApplication&) = default;
    ConcreteParameterApplication& operator=(ConcreteParameterApplication&&) = default;
    explicit ConcreteParameterApplication(const Scheme& s) : SchemeBase{s} {}
    ~ConcreteParameterApplication() override = default;

    [[nodiscard]] Maybe_error<Result> operator()(const ParamValues& p) const override {
        return SchemeBase::scheme()(p);
    }

    // [[nodiscard]] Maybe_error<ParamValues> operator()(const Result& r) const override {
    //     return Base::scheme()(r);
    // }
};

// ─────────────────────────────────────────────────────────────────────────────
// 5) ConcreteModel: ties it all together
// ─────────────────────────────────────────────────────────────────────────────
template <typename Scheme, typename... ParamValues>
class ConcreteModel
    : public IModel<ParamValues...>,
      public virtual ConcreteBaseScheme<Scheme>,
      public ConcreteParameterApplication<
          Scheme, ParamValues, macrodr::Transfer_Op_to<ParamValues, macrodr::Patch_Model>>... {
    using SchemeBase = ConcreteBaseScheme<Scheme>;

   public:
    using ConcreteParameterApplication<
        Scheme, ParamValues,
        macrodr::Transfer_Op_to<ParamValues, macrodr::Patch_Model>>::operator()...;
    explicit ConcreteModel(Scheme s)
        : SchemeBase(std::move(s)),
          ConcreteParameterApplication<Scheme, ParamValues,
                                       macrodr::Transfer_Op_to<ParamValues, macrodr::Patch_Model>>(
              SchemeBase::scheme())... {}
    // In ConcreteModel<T...>

    ConcreteModel(const ConcreteModel&) = default;
    ConcreteModel& operator=(const ConcreteModel&) = default;
    // keep moves:
    ConcreteModel(ConcreteModel&&) = default;
    ConcreteModel& operator=(ConcreteModel&&) = default;

    // clone(): reconstruct from scheme (copying Scheme)
    [[nodiscard]] std::unique_ptr<IModel<ParamValues...>> clone() const override {
        return std::make_unique<ConcreteModel>(SchemeBase::scheme());
    }

    ~ConcreteModel() override = default;

    [[nodiscard]] std::string model_name() const override {
        return SchemeBase::scheme().model_name();
    }

    [[nodiscard]] const std::vector<std::string>& names() const override {
        return SchemeBase::scheme().names();
    }

    [[nodiscard]] std::size_t number_of_states() const override {
        return SchemeBase::scheme().number_of_states();
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// 6) Factory
// ─────────────────────────────────────────────────────────────────────────────
template <typename Scheme, typename... ParamValues>
auto make_model_interface(Scheme scheme) {
    return std::make_unique<ConcreteModel<Scheme, ParamValues...>>(std::move(scheme));
}

}  // namespace macrodr::interface
