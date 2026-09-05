// Gates for legacy/acquisition_filter.h — phase 0 of the MacroIR-Bessel plan
// (theory/macroir/notes/macroir_bessel_plan.md, milestone M0).
//
// Gate map:
//   [coeffs]    exact reverse-Bessel coefficients (independent identities:
//               θ_n(0) = (2n−1)!!, a1 = a0).
//   [roots]     pole residuals, stability, conjugacy, root sums.
//   [cutoff]    −3 dB at f_c after scaling.
//   [dc]        H(0) = 1 from residues + ∫h by quadrature + Laplace spot checks.
//   [dd]        dd1/dd2 against quadrature of their defining integrals,
//               including the near-coincident stress; dd1_pair (real-pair,
//               derivative-mode seed) against complex dd1; ν=0 reproduces the
//               classic real divided differences (the qmodel Ee/E3 algebra).
//   [noise]     closed-form read variance / lag covariances against direct
//               quadrature, the Var+2ΣCov = S0/Δ conservation, and the
//               manuscript appendix bound table (61.5% deficit, 0.641 lag-1
//               correlation at f_c·Δ = 0.2; 2.9%/1.5% at f_c·Δ = 5).

#include <cmath>
#include <complex>

#include "acquisition_filter.h"
#include "catch_amalgamated.hpp"

using namespace macrodr::acqf;
using std::abs;

namespace {

// composite Simpson on [a,b]
template <class F>
double simpson(F&& f, double a, double b, int n) {
    if (n % 2)
        ++n;
    const double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; ++i) s += f(a + i * h) * (i % 2 ? 4.0 : 2.0);
    return s * h / 3.0;
}

cdouble simpson_c(auto&& f, double a, double b, int n) {
    if (n % 2)
        ++n;
    const double h = (b - a) / n;
    cdouble s = f(a) + f(b);
    for (int i = 1; i < n; ++i) s += f(a + i * h) * (i % 2 ? 4.0 : 2.0);
    return s * h / 3.0;
}

double dfact_odd(int n) {  // (2n−1)!!
    double r = 1.0;
    for (int k = 1; k <= 2 * n - 1; k += 2) r *= k;
    return r;
}

// classic real divided differences, transcribing the qmodel.h algebra
// (E1/Ee/E111 well-separated branches) as the ν=0 anchor at kernel level
double E1_ref(double x) {
    return std::abs(x) < 1e-14 ? 1.0 : std::expm1(x) / x;
}
double Ee_ref(double x, double y) {
    return std::abs(x - y) < 1e-9 ? std::exp(x) : (std::exp(x) - std::exp(y)) / (x - y);
}
double E111_ref(double x, double y, double z) {
    return std::exp(x) / ((x - y) * (x - z)) + std::exp(y) / ((y - x) * (y - z)) +
           std::exp(z) / ((z - x) * (z - y));
}

}  // namespace

TEST_CASE("reverse Bessel coefficients are the exact integers", "[acqf][coeffs]") {
    auto c4 = reverse_bessel_coeffs(4);
    REQUIRE(c4.size() == 5);
    CHECK(c4[0] == 105.0);
    CHECK(c4[1] == 105.0);
    CHECK(c4[2] == 45.0);
    CHECK(c4[3] == 10.0);
    CHECK(c4[4] == 1.0);
    for (int n : {2, 3, 4, 5, 6, 7, 8}) {
        auto c = reverse_bessel_coeffs(n);
        CHECK(c[0] == dfact_odd(n));  // θ_n(0) = (2n−1)!!
        CHECK(c[1] == c[0]);          // a1 = a0 (unit-delay normalization)
        CHECK(c[n] == 1.0);           // monic
    }
}

TEST_CASE("poles: residual, stability, conjugacy, sums", "[acqf][roots]") {
    for (int n : {4, 8}) {
        auto coeffs = reverse_bessel_coeffs(n);
        auto roots = poly_roots(coeffs);
        REQUIRE(roots.size() == std::size_t(n));
        double scale = coeffs[0];
        cdouble sum = 0.0;
        for (auto& p : roots) {
            CHECK(abs(poly_eval(coeffs, p)) < 1e-8 * scale);
            CHECK(p.real() < 0.0);
            sum += p;
            // conjugate partner present
            double best = 1e300;
            for (auto& q : roots) best = std::min(best, abs(q - std::conj(p)));
            CHECK(best < 1e-9);
        }
        CHECK(abs(sum - cdouble(-coeffs[n - 1], 0.0)) < 1e-9 * coeffs[n - 1]);
    }
}

TEST_CASE("-3 dB sits at f_c after scaling", "[acqf][cutoff]") {
    for (int n : {4, 8}) {
        const double fc = 10e3;
        auto f = make_bessel_filter(n, fc);
        REQUIRE(int(f.pairs.size()) == n / 2);
        // |H(i·2π·fc)|² = Π |p_k|² / Π |iω − p_k|²  with unit DC gain
        cdouble num = 1.0, den = 1.0;
        const cdouble iw(0.0, 2.0 * std::numbers::pi * fc);
        for (auto& p : f.poles) {
            num *= (-p);
            den *= (iw - p);
        }
        const double mag2 = std::norm(num / den);
        CHECK(abs(mag2 - 0.5) < 1e-9);
        auto err = validate(f);
        CHECK(!err.has_value());
    }
}

TEST_CASE("DC gain and impulse response", "[acqf][dc]") {
    auto f = make_bessel_filter(4, 10e3);
    CHECK(abs(impulse_area(f) - 1.0) < 1e-12);

    // ∫h dt = 1 by quadrature; h decays within ~1 ms at fc = 10 kHz
    const double T = 1.5e-3;
    const double area = simpson([&](double t) { return h_impulse(f, t); }, 0.0, T, 40000);
    CHECK(abs(area - 1.0) < 1e-6);

    // Laplace spot checks: ∫ h e^{−st} dt = θ(0)/θ_phys(s) at a few real s
    auto coeffs = reverse_bessel_coeffs(4);
    const double sc = 2.0 * std::numbers::pi * f.f_c_hz / f.omega3_norm;
    for (double s : {1e3, 2e4, 1e5}) {
        const double lap =
            simpson([&](double t) { return h_impulse(f, t) * std::exp(-s * t); }, 0.0, T, 40000);
        const double Href = (coeffs[0] / poly_eval(coeffs, cdouble(s / sc, 0.0))).real();
        CHECK(abs(lap - Href) < 1e-6);
    }
}

TEST_CASE("dd1 against its defining integral, incl. complex shift", "[acqf][dd]") {
    // ∫₀¹ e^{x s} e^{y(1−s)} ds = dd1[x, y]
    const double xs[] = {-5.0, -1.0, -0.1, 0.0};
    const cdouble ys[] = {{-4.0, 0.0}, {-0.5, 0.3}, {0.0, 3.0}, {-2.1, -2.65}};
    for (double x : xs)
        for (cdouble y : ys) {
            const cdouble quad =
                simpson_c([&](double s) { return std::exp(x * s) * std::exp(y * (1.0 - s)); }, 0.0,
                          1.0, 4000);
            CHECK(abs(dd1(cdouble(x, 0.0), y) - quad) < 1e-10);
        }
}

TEST_CASE("dd1 near-coincidence is smooth and exact", "[acqf][dd]") {
    const double x = -1.3;
    for (double d : {1e-14, 1e-10, 1e-7, 1e-4, 1e-2}) {
        // real coincidence
        const cdouble v = dd1(cdouble(x, 0.0), cdouble(x + d, 0.0));
        const double ref =
            std::exp(x + 0.5 * d) * (1.0 + d * d / 24.0 + d * d * d * d / 1920.0);  // e^{m}·sinhc(d/2)
        CHECK(abs(v.real() - ref) < 1e-13 * std::abs(ref) + 1e-16);
        // complex coincidence
        const cdouble vc = dd1(cdouble(x, 0.7), cdouble(x, 0.7 + d));
        CHECK(std::isfinite(vc.real()));
        CHECK(std::isfinite(vc.imag()));
        CHECK(abs(vc - std::exp(cdouble(x, 0.7 + 0.5 * d))) < 1e-3 * d + 1e-13);
    }
}

TEST_CASE("dd1_pair (real-pair algebra) equals complex dd1", "[acqf][dd]") {
    const double xs[] = {-6.0, -1.0, 0.0, -0.001};
    const double yrs[] = {-3.0, -0.5, 0.0};
    const double yis[] = {0.0, 1e-9, 0.4, 5.0, -2.7};
    for (double x : xs)
        for (double yr : yrs)
            for (double yi : yis) {
                double re, im;
                dd1_pair(x, yr, yi, re, im);
                const cdouble ref = dd1(cdouble(x, 0.0), cdouble(yr, yi));
                CHECK(abs(re - ref.real()) < 1e-14 + 1e-13 * abs(ref));
                CHECK(abs(im - ref.imag()) < 1e-14 + 1e-13 * abs(ref));
            }
}

TEST_CASE("dd2 against 2D quadrature (Hermite–Genocchi)", "[acqf][dd]") {
    // dd2[x,y,z] = ∫₀¹∫₀^{t} e^{x u + y (t−u) + z (1−t)} du dt
    auto quad2 = [](cdouble x, cdouble y, cdouble z) {
        const int N = 400;
        cdouble s = 0.0;
        const double h = 1.0 / N;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j <= i; ++j) {
                const double t = (i + 0.5) * h, u = (j + 0.5) * h;
                if (u <= t)
                    s += std::exp(x * u + y * (t - u) + z * (1.0 - t));
            }
        return s * h * h;
    };
    const cdouble args[] = {{-2.0, 0.0}, {-0.3, 1.7}, {0.0, 0.0}, {-1.1, -2.4}};
    for (cdouble x : args)
        for (cdouble y : args)
            for (cdouble z : args)
                CHECK(abs(dd2(x, y, z) - quad2(x, y, z)) < 5e-3);  // midpoint 2D is O(h²)

    // sharp identities
    CHECK(abs(dd2(cdouble(-1.0, 0), cdouble(-1.0, 0), cdouble(-1.0, 0)) -
              0.5 * std::exp(-1.0)) < 1e-14);
    // dd2[x,y,z] = (dd1[x,z] − dd1[y,z])/(x−y), well-separated
    const cdouble x(-2.0, 0.4), y(-0.5, -1.0), z(-3.3, 2.0);
    CHECK(abs(dd2(x, y, z) - (dd1(x, z) - dd1(y, z)) / (x - y)) < 1e-12);
}

TEST_CASE("dd2 near-coincidence stress", "[acqf][dd]") {
    const cdouble x(-1.7, 0.9);
    for (double d : {1e-13, 1e-9, 1e-5, 1e-2}) {
        // dd2 ≈ (1/2)·e^{centroid}·(1 + O(spread²)); centroid = x + d/6 here
        const cdouble a = dd2(x, x + d, x - 0.5 * d);
        CHECK(abs(a - 0.5 * std::exp(x + d / 6.0)) < d * d + 1e-12);
        // continuity against the fully coincident value
        const cdouble b = dd2(x, x, x);
        CHECK(abs(a - b) < 2.0 * d);
    }
}

TEST_CASE("zero shift reproduces the classic real divided differences", "[acqf][dd]") {
    // dd1 ≡ Ee, dd2 ≡ E111 on real, well-separated arguments: the ν = 0
    // anchor that chains the new kernels to the box member's algebra.
    const double xs[] = {-4.0, -1.0, -0.2, 0.0};
    for (double a : xs)
        for (double b : xs) {
            if (abs(a - b) < 1e-9)
                continue;
            CHECK(abs(dd1(cdouble(a, 0), cdouble(b, 0)).real() - Ee_ref(a, b)) <
                  1e-13 * (1.0 + abs(Ee_ref(a, b))));
        }
    CHECK(abs(dd2(cdouble(-4, 0), cdouble(-1, 0), cdouble(-0.2, 0)).real() -
              E111_ref(-4.0, -1.0, -0.2)) < 1e-13);
    CHECK(abs(dd1(cdouble(-3, 0), cdouble(0, 0)).real() - E1_ref(-3.0)) < 1e-14);
}

TEST_CASE("exp_rot equals the complex exponential", "[acqf][rot]") {
    for (double sdt : {0.0, 0.3, 2.0})
        for (double wdt : {0.0, 0.5, 4.0}) {
            auto M = exp_rot(sdt, wdt);
            const cdouble e = std::exp(cdouble(-sdt, wdt));
            // (u,v) = (2Re x, 2Im x): u' = Re(e·x)·2, matches M·(u,v)
            CHECK(abs(M[0] - e.real()) < 1e-15);
            CHECK(abs(M[1] + e.imag()) < 1e-15);
            CHECK(abs(M[2] - e.imag()) < 1e-15);
            CHECK(abs(M[3] - e.real()) < 1e-15);
        }
}

TEST_CASE("noise closed forms against quadrature", "[acqf][noise]") {
    auto f = make_bessel_filter(4, 10e3);
    const double S0 = 1.0;

    // R(0) = S0·∫ h² dt
    const double h2 =
        simpson([&](double t) { double h = h_impulse(f, t); return h * h; }, 0.0, 1.5e-3, 40000);
    CHECK(abs(noise_autocov(f, S0, 0.0) - h2) < 1e-6 * h2);

    // R(τ) vs ∫ h(u)h(u+τ)du at a few τ
    for (double tau : {1e-5, 5e-5, 2e-4}) {
        const double q = simpson([&](double u) { return h_impulse(f, u) * h_impulse(f, u + tau); },
                                 0.0, 1.5e-3, 40000);
        CHECK(abs(noise_autocov(f, S0, tau) - q) < 1e-6 * h2);
    }

    // Var(z) = (2/Δ²)∫₀^Δ (Δ−τ)R(τ)dτ
    const double dt = 0.2 / f.f_c_hz;  // f_c·Δ = 0.2
    const double varq =
        simpson([&](double tau) { return (dt - tau) * noise_autocov(f, S0, tau); }, 0.0, dt,
                20000) *
        2.0 / (dt * dt);
    CHECK(abs(noise_read_variance(f, S0, dt) - varq) < 1e-8 * S0 / dt);

    // Cov(z_j, z_{j+1}) = (1/Δ²)∫₀^Δ∫₀^Δ R(Δ+u−v) du dv
    const int N = 400;
    double covq = 0.0;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            covq += noise_autocov(f, S0, dt + (i + 0.5) * dt / N - (j + 0.5) * dt / N);
    covq /= double(N) * double(N);
    CHECK(abs(noise_read_cov(f, S0, dt, 1) - covq) < 1e-4 * S0 / dt);
}

TEST_CASE("manuscript bound table and conservation", "[acqf][noise][paper]") {
    auto f = make_bessel_filter(4, 10e3);
    const double S0 = 1.0;

    auto deficit = [&](double fcdt) {
        const double dt = fcdt / f.f_c_hz;
        return 1.0 - noise_read_variance(f, S0, dt) * dt / S0;
    };
    auto rho1 = [&](double fcdt) {
        const double dt = fcdt / f.f_c_hz;
        return noise_read_cov(f, S0, dt, 1) / noise_read_variance(f, S0, dt);
    };

    // appendix table (08_appendix_derivation.tex:1016-1017)
    CHECK(abs(deficit(0.2) - 0.615) < 0.002);
    CHECK(abs(rho1(0.2) - 0.641) < 0.002);
    CHECK(abs(deficit(5.0) - 0.029) < 0.002);
    CHECK(abs(rho1(5.0) - 0.015) < 0.002);
    CHECK(abs(deficit(20.0) - 0.0074) < 0.0005);
    CHECK(abs(rho1(20.0) - 0.0037) < 0.0005);
    // asymptotic ~1/(f_c·Δ) decay continues past the table
    CHECK(deficit(50.0) < 0.5 * deficit(20.0));
    CHECK(rho1(50.0) < 0.5 * rho1(20.0));

    // conservation: Var + 2·Σ_{m≥1} Cov_m = S0/Δ (grouping identity)
    for (double fcdt : {0.2, 1.0, 5.0}) {
        const double dt = fcdt / f.f_c_hz;
        double total = noise_read_variance(f, S0, dt);
        for (int m = 1; m <= 200; ++m) total += 2.0 * noise_read_cov(f, S0, dt, m);
        CHECK(abs(total - S0 / dt) < 1e-6 * S0 / dt);
    }
}
