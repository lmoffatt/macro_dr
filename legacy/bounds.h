#ifndef BOUNDS_H
#define BOUNDS_H

#include "moment_statistics.h"   // value_type_t<Id>; same utility used by count<Va>/Sum<Va>
#include "variables.h"           // var::Constant

// -----------------------------------------------------------------------------
// Generic upper/lower bound tag templates.
//
//   Max<Id>   : maximum observed (running upper bound during a search;
//               result-side, e.g. Max<logL> = max_logL achieved).
//   Min<Id>   : minimum observed (running lower bound).
//   Cap<Id>   : maximum allowed (configured upper limit;
//               e.g. Cap<lambda> = lambda_max in LM, Cap<N_Iterations> = max_iter).
//   Floor<Id> : minimum allowed (configured lower limit;
//               e.g. Floor<lambda> = lambda_min).
//
// Payload type is value_type_t<Id>, so the wrapped value automatically matches
// the Id's natural numeric type:
//   Cap<lambda>          → payload double
//   Cap<N_Iterations>    → payload std::size_t
//   Max<logL>            → payload double
//
// Same mechanic as count<Va> (line 217) and Sum<Va> (line 432) in
// moment_statistics.h: a tagged var::Constant whose payload is derived from Id.
//
// Not related to <limits> / std::numeric_limits — these are var::Constant tag
// types for use in Vector_Space-based config/result records, not standard-
// library numeric introspection.
// -----------------------------------------------------------------------------

template <class Id>
struct Max : public var::Constant<Max<Id>, value_type_t<Id>> {
    using value_type = value_type_t<Id>;
    using base_type  = var::Constant<Max<Id>, value_type>;
    using base_type::base_type;
    Max() : base_type{value_type{}} {}
    Max(value_type v) : base_type{std::move(v)} {}
    Max(Id const& x) : base_type{x()} {}
};

template <class Id>
struct Min : public var::Constant<Min<Id>, value_type_t<Id>> {
    using value_type = value_type_t<Id>;
    using base_type  = var::Constant<Min<Id>, value_type>;
    using base_type::base_type;
    Min() : base_type{value_type{}} {}
    Min(value_type v) : base_type{std::move(v)} {}
    Min(Id const& x) : base_type{x()} {}
};

template <class Id>
struct Cap : public var::Constant<Cap<Id>, value_type_t<Id>> {
    using value_type = value_type_t<Id>;
    using base_type  = var::Constant<Cap<Id>, value_type>;
    using base_type::base_type;
    Cap() : base_type{value_type{}} {}
    Cap(value_type v) : base_type{std::move(v)} {}
};

template <class Id>
struct Floor : public var::Constant<Floor<Id>, value_type_t<Id>> {
    using value_type = value_type_t<Id>;
    using base_type  = var::Constant<Floor<Id>, value_type>;
    using base_type::base_type;
    Floor() : base_type{value_type{}} {}
    Floor(value_type v) : base_type{std::move(v)} {}
};

#endif  // BOUNDS_H
