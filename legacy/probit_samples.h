#ifndef PROBIT_SAMPLES_H
#define PROBIT_SAMPLES_H

// Probit-percentile sampling of per-group MLE results.
//
// Instead of saving every per-group raw state (heavy: 2 MB × N_groups for
// the rich state), we save a small number of "representative" groups picked
// by ranking each group on a chosen scalar and snapshotting the groups
// closest to a few probit-percentile heights (typically {2.5%, 16%, 50%,
// 84%, 97.5%} — extremes, ±1σ, median).
//
// Multiple ranking variables can be applied to the same N_groups in one
// pass — each produces its own set of probit samples. Defaults:
//   * Maximum_logL              (proxy for "well-fit" vs "poorly-fit" group)
//   * F_Lambda_Min               (spectral safety; <0 = BIASED)
//   * Wald_T2 vs F-metric        (distance from θ_0 weighted by F)
//
// Memory per saved sample = sizeof(State) + a few hundred bytes of metadata.
// For State = dMacro_State_Hessian_minimal_param: ~1 KB per sample, ~5 KB
// per ranking variable, ~15 KB for 3 variables. Negligible.
//
// For State = dMacro_State_Ev_gradient_all_param (with per-step Evolution_of
// of length T × group_size): ~2 MB × group_size per sample. At group_size=1
// and 5 percentiles × 3 ranking variables: 30 MB. Still feasible.

#include <cstddef>
#include <string>
#include <vector>

#include "distributions.h"
#include "matrix.h"
#include "parameter_indexed.h"
#include "variables.h"

namespace macrodr {

// One representative group snapshot at a specific (ranking_variable, probit
// height) coordinate. State is the same template parameter the surrounding
// command was instantiated on (minimal or rich), so the snapshot resolution
// matches whatever the user paid for in the run.
template <class State>
struct Probit_Sample_Record {
    std::size_t              group_id;                       // index into the per-group result vector
    std::vector<std::size_t> recording_indices;                // groups[group_id]'s recording indices (refs into input)
    State                    state_at_theta_hat;               // dMacro_State_..._param for this group
    parameter_spd_payload    F;                                // numerical Fisher at θ̂_group
    Matrix<double>           F_eigenvalues_sorted_asc;         // p×1
    double                   F_lambda_min{0.0};                // = F_eigenvalues_sorted_asc(0)
    std::string              safety_categorization;            // "SAFE" / "MARGINAL" / "UNRELIABLE" / "BIASED"
    double                   ranking_value{0.0};               // the value of the ranking variable for this group
    double                   probit_height{0.0};               // the probit percentile this group represents (0.025, 0.5, …)
    std::string              ranking_variable_name;            // "Maximum_logL" / "F_Lambda_Min" / "Wald_T2" / …
};

// All probit samples collected at a given group_size. Multiple ranking
// variables can populate this vector — each contributes a small set of
// records (one per probit height). Distinguish by `ranking_variable_name`.
template <class State>
class Probit_Samples_at_Group_Size
    : public var::Constant<Probit_Samples_at_Group_Size<State>,
                           std::vector<Probit_Sample_Record<State>>> {
    using base_type = var::Constant<Probit_Samples_at_Group_Size<State>,
                                    std::vector<Probit_Sample_Record<State>>>;

   public:
    using base_type::base_type;
    Probit_Samples_at_Group_Size() = default;

    friend std::string className(Probit_Samples_at_Group_Size) {
        return "Probit_Samples_at_Group_Size";
    }
};

}  // namespace macrodr

#endif  // PROBIT_SAMPLES_H
