#include <catch_amalgamated.hpp>

#include <macrodr/io/fingerprint.h>

#include <cstdio>
#include <filesystem>
#include <limits>

namespace fp = macrodr::io::fingerprint;

TEST_CASE("fingerprint hashing basics", "[fingerprint]") {
    SECTION("same inputs hash same") {
        REQUIRE(fp::hash_of(std::string("hello")) == fp::hash_of(std::string("hello")));
        REQUIRE(fp::hash_of(42) == fp::hash_of(42));
        REQUIRE(fp::hash_of(3.14) == fp::hash_of(3.14));
    }

    SECTION("different inputs hash different") {
        REQUIRE(fp::hash_of(std::string("hello")) != fp::hash_of(std::string("world")));
        REQUIRE(fp::hash_of(1) != fp::hash_of(2));
        REQUIRE(fp::hash_of(1.0) != fp::hash_of(2.0));
    }

    SECTION("separator prevents string-boundary collisions") {
        fp::Hasher a;
        a.update(std::string("abc"));
        a.update(std::string(""));
        fp::Hasher b;
        b.update(std::string("ab"));
        b.update(std::string("c"));
        REQUIRE(a.finalize() != b.finalize());
    }

    SECTION("-0.0 hashes same as +0.0") {
        REQUIRE(fp::hash_of(0.0) == fp::hash_of(-0.0));
    }

    SECTION("NaNs canonicalize") {
        double n1 = std::nan("1");
        double n2 = std::nan("2");
        REQUIRE(fp::hash_of(n1) == fp::hash_of(n2));
    }

    SECTION("field helper differentiates names") {
        fp::Hasher a;
        a.update_field("alpha", 1.0);
        fp::Hasher b;
        b.update_field("beta", 1.0);
        REQUIRE(a.finalize() != b.finalize());
    }

    SECTION("hash composition chains") {
        auto inner = fp::hash_of(std::string("payload"));
        fp::Hasher outer;
        outer.update_field("child", inner);
        outer.update_field("seed", 0);
        auto result = outer.finalize();
        REQUIRE(result.raw() != 0);
    }
}

TEST_CASE("fingerprint hex round-trip", "[fingerprint]") {
    auto f = fp::hash_of(std::string("payload"));
    auto hex = f.to_hex();
    REQUIRE(hex.size() == 16);
    auto parsed = fp::Fingerprint::from_hex(hex);
    REQUIRE(parsed);
    REQUIRE(parsed.value() == f);

    REQUIRE_FALSE(fp::Fingerprint::from_hex("short"));
    REQUIRE_FALSE(fp::Fingerprint::from_hex("zzzzzzzzzzzzzzzz"));
}

TEST_CASE("sidecar write and read round-trip", "[fingerprint]") {
    auto tmp = std::filesystem::temp_directory_path() / "macrodr_fingerprint_test.csv";
    std::filesystem::remove(fp::sidecar_path_for(tmp));

    fp::Sidecar w;
    w.fingerprint = fp::hash_of(std::string("payload"));
    w.scheme = "test_scheme:v1";
    w.metadata["note"] = macrodr::io::json::Json::string("hand-written test");
    w.metadata["count"] = macrodr::io::json::Json::number(7);

    auto wrote = fp::write_sidecar(tmp, w);
    REQUIRE(wrote);

    auto r = fp::read_sidecar(tmp);
    REQUIRE(r);
    REQUIRE(r.value().fingerprint == w.fingerprint);
    REQUIRE(r.value().scheme == "test_scheme:v1");

    REQUIRE(fp::sidecar_fingerprint_matches(tmp, w.fingerprint));
    auto other = fp::hash_of(std::string("different"));
    REQUIRE_FALSE(fp::sidecar_fingerprint_matches(tmp, other));

    std::filesystem::remove(fp::sidecar_path_for(tmp));
}

TEST_CASE("read_sidecar reports missing file", "[fingerprint]") {
    auto tmp = std::filesystem::temp_directory_path() / "macrodr_fingerprint_absent.csv";
    std::filesystem::remove(fp::sidecar_path_for(tmp));
    auto r = fp::read_sidecar(tmp);
    REQUIRE_FALSE(r);
    REQUIRE_FALSE(fp::sidecar_fingerprint_matches(tmp, fp::Fingerprint{}));
}
