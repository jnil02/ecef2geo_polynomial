// Copyright (c) 2024, John-Olof Nilsson
// Distributed under the terms of the BSD 2-Clause License.

#include <random>

// This tells Catch to provide a main() - only do this in one cpp file
#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "util_geo.hpp"

#include "../transformation_implementations/ecef2geo_polynomial.hpp"
#include "../coefficient_computations/util_mpgeo.hpp"

// Precision used for the multiprecision transformations.
constexpr mp_prec_t prec = 200;
// Mpfr consts with appropriate precision.
wgs84_mpfr_constants util_mpgeo::mpfr_consts = wgs84_mpfr_constants(prec);

using namespace util_geo;

/** ECEF to lla multiprecision reference transformation wrapper.
 *
 * @param ecef {x, y, z}
 * @return {latitude, longitude, altitude}
 */
static lla ecef2lla(xyz ecef) {
	lla_mpfr lla_mpfr;
	mpfr_init2(lla_mpfr.lat, prec);
	mpfr_init2(lla_mpfr.lon, prec);
	mpfr_init2(lla_mpfr.alt, prec);
	xyz_mpfr ecef_mpfr;
	mpfr_init2(ecef_mpfr.x, prec);
	mpfr_init2(ecef_mpfr.y, prec);
	mpfr_init2(ecef_mpfr.z, prec);

	mpfr_set_d(ecef_mpfr.x, ecef.x, MPFR_RNDN);
	mpfr_set_d(ecef_mpfr.y, ecef.y, MPFR_RNDN);
	mpfr_set_d(ecef_mpfr.z, ecef.z, MPFR_RNDN);
	util_mpgeo::ecef2lla(lla_mpfr, ecef_mpfr);
	lla _lla = {mpfr_get_d1(lla_mpfr.lat), mpfr_get_d1(lla_mpfr.lon),
				mpfr_get_d1(lla_mpfr.alt)};

	mpfr_clear(ecef_mpfr.x);
	mpfr_clear(ecef_mpfr.y);
	mpfr_clear(ecef_mpfr.z);
	mpfr_clear(lla_mpfr.lat);
	mpfr_clear(lla_mpfr.lon);
	mpfr_clear(lla_mpfr.alt);
	return _lla;
}

/** ECEF to nva multiprecision reference transformation wrapper.
 *
 * @param ecef {x, y, z}
 * @return {n-vector, altitude}
 */
static nva ecef2nva(xyz _ecef) {
	nva_mpfr nva_mpfr;
	mpfr_init2(nva_mpfr.n.x, prec);
	mpfr_init2(nva_mpfr.n.y, prec);
	mpfr_init2(nva_mpfr.n.z, prec);
	mpfr_init2(nva_mpfr.alt, prec);
	xyz_mpfr ecef;
	mpfr_init2(ecef.x, prec);
	mpfr_init2(ecef.y, prec);
	mpfr_init2(ecef.z, prec);

	mpfr_set_d(ecef.x, _ecef.x, MPFR_RNDN);
	mpfr_set_d(ecef.y, _ecef.y, MPFR_RNDN);
	mpfr_set_d(ecef.z, _ecef.z, MPFR_RNDN);
	util_mpgeo::ecef2nva(nva_mpfr, ecef);
	nva _lla = {{mpfr_get_d1(nva_mpfr.n.x), mpfr_get_d1(nva_mpfr.n.y),
				 mpfr_get_d1(nva_mpfr.n.z)}, mpfr_get_d1(nva_mpfr.alt)};

	mpfr_clear(ecef.x);
	mpfr_clear(ecef.y);
	mpfr_clear(ecef.z);
	mpfr_clear(nva_mpfr.n.x);
	mpfr_clear(nva_mpfr.n.y);
	mpfr_clear(nva_mpfr.n.z);
	mpfr_clear(nva_mpfr.alt);
	return _lla;
}

static void test_lla(lla geo, lla ref, double horz_margin, double alt_margin) {
	REQUIRE(geo.lat == Approx(ref.lat).margin(horz_margin));
	REQUIRE(geo.lon == Approx(ref.lon).margin(horz_margin));
	REQUIRE(geo.alt == Approx(ref.alt).margin(alt_margin));
}

static void test_nva(nva geo, nva ref, double margin_n, double margin_alt) {
	REQUIRE(geo.n.x == Approx(ref.n.x).margin(margin_n));
	REQUIRE(geo.n.y == Approx(ref.n.y).margin(margin_n));
	REQUIRE(geo.n.z == Approx(ref.n.z).margin(margin_n));
	REQUIRE(geo.alt == Approx(ref.alt).margin(margin_alt));
}

TEST_CASE("Transformation of ECEF to lla representation.") {
	using ecef2geo::OMEGA_H_MIN;
	using ecef2geo::OMEGA_H_MAX;

	auto test_trans = &ecef2geo::ecef2lla_asina<2,1,2,0,7,8>;
	auto ref_double = &vermeille::ecef2lla;
	auto ref_mpfr = &ecef2lla;
	double horz_margin = 0.41 / (wgs84::b + OMEGA_H_MIN);
	double alt_margin = 0.41;
	int nr_samples = 1e4;

	SECTION("Transformation of a single coordinate.") {
		xyz ecef = lla2ecef({27.986065 / 180. * M_PI, 86.922623 / 180. * M_PI,
							 8848.});
		lla geo = test_trans(ecef);  // Max error ~0.41m.
		lla geo_ref = ref_mpfr(ecef);  // Reference transformation.
		test_lla(geo, geo_ref, horz_margin, alt_margin);
	}

	SECTION("Transformation of coordinates within approximation limits.") {
//		std::random_device rd(0);
//		std::mt19937 gen(rd);
		std::mt19937 gen(0);
		std::uniform_real_distribution<> lat_dist(-M_PI_2, M_PI_2);
		std::uniform_real_distribution<> lon_dist(-M_PI, M_PI);
		std::uniform_real_distribution<> alt_dist(OMEGA_H_MIN, OMEGA_H_MAX);
		for (int i = 0; i < nr_samples; ++i) {
			xyz ecef = lla2ecef({lat_dist(gen), lon_dist(gen), alt_dist(gen)});
			lla geo = test_trans(ecef);  // Max error ~0.41m.
			lla geo_ref = ref_mpfr(ecef);  // Reference transformation.
			test_lla(geo, geo_ref, horz_margin, alt_margin);
		}
	}

	/*
	 * This is a very inefficient test of reference transformation. It is
	 * inefficient in the sense that samples will be much denser towards the
	 * center of Earth. However, this is fine for this test.
	 *
	 * Note, all values are tested except those i in the singular disc. The
	 * singular disc requires special tests. The polynomial approximations are
	 * not suitable for values in the singular disc so not testing them are
	 * fine.
	 */
	SECTION("Reference transformation test.") {
		double max_alt = 1e6;
//		std::random_device rd(0);
//		std::mt19937 gen(rd);
		std::mt19937 gen(0);
		std::uniform_real_distribution<> lat_dist(-M_PI_2, M_PI_2);
		std::uniform_real_distribution<> lon_dist(-M_PI, M_PI);
		std::uniform_real_distribution<> norm_alt_dist(0, 1);
		for (int i = 0; i < nr_samples; ++i) {
			double lat = lat_dist(gen);
			// All valid altitude including inside evolute.
			double alt = normAlt2Alt(norm_alt_dist(gen), max_alt, lat);
			lla geo = {lat, lon_dist(gen), alt};
			xyz ecef = lla2ecef(geo);
			lla geo_ref = ref_mpfr(ecef);  // Reference transformation.
			lla geo_ref2 = ref_double(ecef);  // Reference transformation.
			test_lla(geo, geo_ref, 0.0, 0.0);  // Perfect transformation.
			test_lla(geo, geo_ref2, 0.0, 0.0);  // Perfect transformation.
		}
	}

	SECTION("Transformation of special points.") {
		double err_margin_rad = 1e-15;
		double err_margin_alt = 1e-9;
		// North Pole.
		lla np_lla = {M_PI_2, 0.0, 0.0};
		xyz np_ecef = lla2ecef(np_lla);
		lla np_double = ref_double(np_ecef);
		lla np_mpfr = ref_mpfr(np_ecef);
		REQUIRE(np_mpfr.lat == Approx(np_lla.lat).margin(err_margin_rad));
		REQUIRE(np_mpfr.alt == Approx(np_lla.alt).margin(err_margin_alt));
		REQUIRE(np_double.lat == Approx(np_lla.lat).margin(err_margin_rad));
		REQUIRE(np_double.alt == Approx(np_lla.alt).margin(err_margin_alt));
		// South Pole.
		lla sp_lla = {-M_PI_2, 0.0, 0.0};
		xyz sp_ecef = lla2ecef(sp_lla);
		lla sp_double = ref_double(sp_ecef);
		lla sp_mpfr = ref_mpfr(sp_ecef);
		REQUIRE(sp_mpfr.lat == Approx(sp_lla.lat).margin(err_margin_rad));
		REQUIRE(sp_mpfr.alt == Approx(sp_lla.alt).margin(err_margin_alt));
		REQUIRE(sp_double.lat == Approx(sp_lla.lat).margin(err_margin_rad));
		REQUIRE(sp_double.alt == Approx(sp_lla.alt).margin(err_margin_alt));
		// Inside evolute.
		lla evo_lla = {M_PI_2 - 0.1, 0.0, -wgs84::b + 1000.};
		xyz evo_ecef = lla2ecef(evo_lla);
		lla evo_double = ref_double(evo_ecef);
		lla evo_mpfr = ref_mpfr(evo_ecef);
		test_lla(evo_lla, evo_double, err_margin_rad, err_margin_alt);
		test_lla(evo_lla, evo_mpfr, err_margin_rad, err_margin_alt);
		// Centre of Earth.
		lla ce_lla = {M_PI_2, 0.0, -wgs84::b};
		xyz ce_ecef = lla2ecef(ce_lla);
		lla ce_double = ref_double(ce_ecef);
		lla ce_mpfr = ref_mpfr(ce_ecef);
		REQUIRE(ce_mpfr.lat == Approx(ce_lla.lat).margin(err_margin_rad));
		REQUIRE(ce_mpfr.alt == Approx(ce_lla.alt).margin(err_margin_alt));
		REQUIRE(ce_double.lat == Approx(ce_lla.lat).margin(err_margin_rad));
		REQUIRE(ce_double.alt == Approx(ce_lla.alt).margin(err_margin_alt));
		// Point on the singular disk.
		// Assume lat = pi/4 and hence sin(lat) = sqrt(0.5)
		// This imply p = e4/(2 - e2) and h and x follows trivially.
		double p = wgs84::e2 * wgs84::e2 / (2. - wgs84::e2);
		double h = - wgs84::b / wgs84::e * std::sqrt(wgs84::e2 - p);
		double x = std::sqrt(p) * wgs84::a;
		lla sd_lla = {M_PI_4, 0.0, h};
		xyz sd_ecef = {x, 0.0, 0.0};
		lla sd_double = ref_double(sd_ecef);
		lla sd_mpfr = ref_mpfr(sd_ecef);
		test_lla(sd_lla, sd_double, err_margin_rad, err_margin_alt);
		test_lla(sd_lla, sd_mpfr, err_margin_rad, err_margin_alt);
	}
}

TEST_CASE("Transformation of ECEF to nva representation.") {
	using ecef2geo::OMEGA_H_MIN;
	using ecef2geo::OMEGA_H_MAX;

	int nr_samples = 1e4;
	auto test_trans = &ecef2geo::ecef2nva<2,2,1,2,0>;
	auto ref_double = &vermeille::ecef2nva;
	auto ref_mpfr = &ecef2nva;
	// Maximum error for test_trans.
	double nva_margin = 1.25 / wgs84::b;  // ~1.25m
	double alt_margin = 0.3;

	SECTION("Transformation of coordinates within approximation limits.") {
//	std::random_device rd(0);
//	std::mt19937 gen(rd);
		std::mt19937 gen(0);
		std::uniform_real_distribution<> lat_dist(-M_PI_2, M_PI_2);
		std::uniform_real_distribution<> lon_dist(-M_PI, M_PI);
		std::uniform_real_distribution<> alt_dist(OMEGA_H_MIN, OMEGA_H_MAX);
		for (int i = 0; i < nr_samples; ++i) {
			xyz ecef = lla2ecef({lat_dist(gen), lon_dist(gen), alt_dist(gen)});
			nva geo = test_trans(ecef);
			nva geo_ref = ref_double(ecef);
			REQUIRE(std::sqrt(
					(geo.n.x - geo_ref.n.x) * (geo.n.x - geo_ref.n.x)
					+ (geo.n.y - geo_ref.n.y) * (geo.n.y - geo_ref.n.y)
					+ (geo.n.z - geo_ref.n.z) * (geo.n.z - geo_ref.n.z)) <
					nva_margin);
			REQUIRE(geo.alt == Approx(geo_ref.alt).margin(alt_margin));
		}
	}

	SECTION("Reference transformation test.") {
		double max_alt = 1e6;
		double err_margin_rad = 1e-15;  // Roughly numerical limit for double.
		double err_margin_alt = 1e-9;  // Roughly numerical limit for double.
//		std::random_device rd(0);
//		std::mt19937 gen(rd);
		std::mt19937 gen(0);
		std::uniform_real_distribution<> lat_dist(-M_PI_2, M_PI_2);
		std::uniform_real_distribution<> lon_dist(-M_PI, M_PI);
		std::uniform_real_distribution<> norm_alt_dist(0, 1);
		for (int i = 0; i < nr_samples; ++i) {
			double lat = lat_dist(gen);
			// All valid altitude including inside evolute.
			double alt = normAlt2Alt(norm_alt_dist(gen), max_alt, lat);
			lla geo = {lat, lon_dist(gen), alt};
			nva geo_nva = lla2nva(geo);
			xyz ecef = lla2ecef(geo);
			nva geo_ref = ref_mpfr(ecef);
			nva geo_ref2 = ref_double(ecef);
			test_nva(geo_nva, geo_ref, err_margin_rad, err_margin_alt);
			test_nva(geo_nva, geo_ref2, err_margin_rad, err_margin_alt);
		}
	}


	SECTION("Transformation of special points.") {
		double err_margin_rad = 1e-15;  // Roughly numerical limit for double.
		double err_margin_alt = 1e-9;  // Roughly numerical limit for double.
		// North Pole.
		nva np_lla = {{0, 0, 1.0}, wgs84::b};
		xyz np_ecef = nva2ecef(np_lla);
		nva np_double = ref_double(np_ecef);
		nva np_mpfr = ref_mpfr(np_ecef);
		test_nva(np_mpfr, np_lla, err_margin_rad, err_margin_alt);
		test_nva(np_double, np_lla, err_margin_rad, err_margin_alt);
		// South Pole.
		nva sp_lla = {{0, 0, -1.0}, wgs84::b};
		xyz sp_ecef = nva2ecef(sp_lla);
		nva sp_double = ref_double(sp_ecef);
		nva sp_mpfr = ref_mpfr(sp_ecef);
		test_nva(sp_mpfr, sp_lla, err_margin_rad, err_margin_alt);
		test_nva(sp_double, sp_lla, err_margin_rad, err_margin_alt);
		// Inside evolute.
		lla evo_lla = {M_PI_2 - 0.1, 0.0, -wgs84::b + 1000.};
		nva evo_nva = lla2nva(evo_lla);
		xyz evo_ecef = nva2ecef(evo_nva);
		nva evo_double = ref_double(evo_ecef);
		nva evo_mpfr = ref_mpfr(evo_ecef);
		test_nva(evo_mpfr, evo_nva, err_margin_rad, err_margin_alt);
		test_nva(evo_double, evo_nva, err_margin_rad, err_margin_alt);
		// Centre of Earth.
		nva ce_lla = {{0, 0, 1}, -wgs84::b};
		xyz ce_ecef = nva2ecef(ce_lla);
		nva ce_double = ref_double(ce_ecef);
		nva ce_mpfr = ref_mpfr(ce_ecef);
		test_nva(ce_mpfr, ce_lla, err_margin_rad, err_margin_alt);
		test_nva(ce_double, ce_lla, err_margin_rad, err_margin_alt);
		// Point on the singular disk.
		// Assume lat = pi/4 and hence sin(lat) = sqrt(0.5)
		// This imply p = e4/(2 - e2) and h and x follows trivially.
		double p = wgs84::e2 * wgs84::e2 / (2. - wgs84::e2);
		double h = - wgs84::b / wgs84::e * std::sqrt(wgs84::e2 - p);
		double x = std::sqrt(p) * wgs84::a;
		lla sd_lla = {M_PI_4, 0.0, h};
		nva sd_nva = lla2nva(sd_lla);
		xyz sd_ecef = {x, 0.0, 0.0};
		nva sd_double = ref_double(sd_ecef);
		nva sd_mpfr = ref_mpfr(sd_ecef);
		test_nva(sd_mpfr, sd_nva, err_margin_rad, err_margin_alt);
		test_nva(sd_double, sd_nva, err_margin_rad, err_margin_alt);
	}
}
