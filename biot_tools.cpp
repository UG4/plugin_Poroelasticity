/*
 * Copyright (c) 2019-2020:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
 *
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "biot_tools.h"

namespace ug {
namespace Poroelasticity {

//! Bessel functions
double BesselJ0(double x)
{ return boost::math::cyl_bessel_j(0, x); }

double BesselJ1(double x)
{ return boost::math::cyl_bessel_j(1, x); }

#ifdef WITH_JSON
/**
{
	"subsets" = "subset1, subset2",
	"alpha" = 1.0,
	"kappa" = 1.0,
	"phi" = 0.0,
	"lambda" = 0.0,
	"mu" = 0.0,
	"beta" = 0.0,
}
*/
void to_json(nlohmann::json &j, const BiotSubsetParameters &p) {
	j = nlohmann::json{
		{"subsets", p.get_subsets()},
		{"alpha",  p.get_alpha()},
		{"kappa", p.get_kappa()},
		{"phi", p.get_phi()},
		{"lambda", p.get_lambda()},
		{"mu", p.get_mu()},
		{"beta", p.get_beta()}
	};
}

void from_json(const nlohmann::json &j, BiotSubsetParameters &p) {
	p.set_subsets(j.at("subsets").get<std::string>());
	p.set_alpha(j.at("alpha").get<number>());
	j.at("kappa").get_to(p.m_kappa);
	j.at("phi").get_to(p.m_phi);
	j.at("lambda").get_to(p.m_lambda);
	j.at("mu").get_to(p.m_mu);
	j.at("beta").get_to(p.m_beta_uzawa);
}
#endif

/** Value: $ \frac{l^2}{\kappa (\lambda + 2* \mu)}$ */
double DefaultCharTime(const BiotSubsetParameters& p, double length)
{
	nlohmann::json json;
	to_json(json, p);

	std::stringstream ss;
	ss << json;

	UG_LOG("JSON = "<< ss.str() << std::endl);
	double consolidation = p.get_kappa()* (p.get_lambda() + 2*p.get_mu());
	return (length*length)/consolidation;
}





}
}
