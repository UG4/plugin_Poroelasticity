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

#ifndef __BIOT_TOOLS_H__
#define __BIOT_TOOLS_H__


// std libs.
#include <string.h>
#include <vector>

// Boost
#include <boost/math/special_functions/bessel.hpp>

// UG4 dependencies.
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/linker/scale_add_linker.h"


// Plugin dependencies.
#include "../ConvectionDiffusion/fe/convection_diffusion_fe.h"
#include "../ConvectionDiffusion/fe/convection_diffusion_stab_fe.h"
#include "../ConvectionDiffusion/fv1/convection_diffusion_fv1.h"
#include "../SmallStrainMechanics/small_strain_mech.h"
#include "../SmallStrainMechanics/material_laws/hooke.h"

#define WITH_JSON
#ifdef WITH_JSON
#include <nlohmann/json.hpp>
#endif

namespace ug {
namespace Poroelasticity {

//! Bessel functions
double BesselJ0(double x);
double BesselJ1(double x);




struct BiotDiscConfig
{
	BiotDiscConfig(const char* uCmp, const char *pCmp)
	: m_uCmp(uCmp), m_pCmp(pCmp), m_uOrder(2), m_pOrder(1), m_dStab(0.0), m_bSteadyStateMechanics(true)
	{}

	BiotDiscConfig(const char* uCmp, int uorder, const char *pCmp, int porder)
	: m_uCmp(uCmp), m_pCmp(pCmp), m_uOrder(uorder), m_pOrder(porder), m_dStab(0.0), m_bSteadyStateMechanics(true) {}

	BiotDiscConfig(const char* uCmp, int uorder, const char *pCmp, int porder, bool bSteadyStateMechanics)
	: m_uCmp(uCmp), m_pCmp(pCmp), m_uOrder(uorder), m_pOrder(porder), m_dStab(0.0), m_bSteadyStateMechanics(bSteadyStateMechanics) {}


	/// Stabilization parameter (from [0,1]).
	void set_stabilization(double dStab)
	{ m_dStab = dStab; }

	int get_porder() const {return m_pOrder;}
	int get_uorder() const {return m_uOrder;}

	std::string m_uCmp, m_pCmp;
	int m_uOrder, m_pOrder;
	double m_dStab;

	bool m_bSteadyStateMechanics;
};

//! Class for Biot parameters (per subset)
/*
 *
 *
 {

  {"alpha":1.0,"beta":2.799999999999999e-06,"kappa":0.01,"lambda":142857.1428571429,"mu":35714.28571428572,"phi":0.0,"subsets":"INNER"}
 } */
class BiotSubsetParameters
{
public:

	//! Default constructor.
	BiotSubsetParameters() {}

	//! Create from table.
	BiotSubsetParameters(const char *s, number alpha, number kappa, number phi, number lambda, number mu, number beta)
	: m_subsets(s), m_alpha(alpha), m_kappa(kappa), m_phi(phi), m_lambda(lambda), m_mu(mu), m_beta_uzawa(beta)
	{}

	//! Subset
	std::string get_subsets() const { return m_subsets; }
	void set_subsets(std::string subset) { m_subsets = subset; }

	//! Biot coefficient
	number get_alpha() const { return m_alpha; }
	void set_alpha(number alpha)  { m_alpha = alpha; }

	//! Storativity
	number get_phi() const { return m_phi; }
	void set_phi(number phi) { m_phi = phi; }

	//! Diffusion
	number get_kappa() const { return m_kappa; }
	void set_kappa(number kappa) { m_kappa = kappa; }

	//! Elasticity lambda
	number get_lambda() const { return m_lambda; }
	void set_lambda(number lambda) { m_lambda = lambda; }

	//! Elasticity mu
	number get_mu() const { return m_mu; }
	void set_mu(number mu) { m_mu = mu; }

	//! OPTIONAL: Fixed-stress beta
	number get_beta() const { return m_beta_uzawa; }
	void set_beta(number beta) { m_beta_uzawa = beta; }

#ifdef WITH_JSON
	friend void from_json(const nlohmann::json& j, BiotSubsetParameters& p);
#endif

protected:
	std::string m_subsets;

	number m_alpha;
	number m_kappa;

	number m_phi;
	number m_lambda;
	number m_mu;

	number m_beta_uzawa;
};

#ifdef WITH_JSON
void to_json(nlohmann::json &j, const BiotSubsetParameters &p);
void from_json(const nlohmann::json &j, BiotSubsetParameters &p);
#endif



//! Compute characteristic time
double DefaultCharTime(const BiotSubsetParameters& p, double length=1.0);

// Container for data (TODO: Extend to real IElemDisc?)
template <typename TDomain>
class BiotElemDisc
{
public:
	static const int dim = TDomain::dim;
	typedef IElemDisc<TDomain> TElemDisc;
	typedef SmallStrainMechanics::SmallStrainMechanicsElemDisc<TDomain> TSmallStrainMechanics;
	typedef ConvectionDiffusionPlugin::ConvectionDiffusionFE<TDomain> TConvectionDiffusion;
	typedef ScaleAddLinker<number, dim, number> TScaleAddLinkerNumber;

	BiotElemDisc()
	: flowEqDisc(SPNULL), displacementEqDisc(SPNULL) {}

	BiotElemDisc(SmartPtr<TConvectionDiffusion> pDisc, SmartPtr<TSmallStrainMechanics> uDisc)
	: flowEqDisc(pDisc), displacementEqDisc(uDisc) {}

	// SmartPtr<TConvectionDiffusion> pressure_disc()
 	SmartPtr<TElemDisc> pressure_disc()
	{return flowEqDisc; }

	// SmartPtr<TSmallStrainMechanics> displacement_disc()
	SmartPtr<TElemDisc> displacement_disc()
	{return displacementEqDisc; }

	ConstSmartPtr<TScaleAddLinkerNumber> compression_linker()
	{return compressionLinker; }

	ConstSmartPtr<TScaleAddLinkerNumber> divergence()
	{return divLinker; }

protected:
	// Container data.
	SmartPtr<TConvectionDiffusion> flowEqDisc;
	SmartPtr<TSmallStrainMechanics> displacementEqDisc;


	SmartPtr<TScaleAddLinkerNumber> divLinker;
	SmartPtr<TScaleAddLinkerNumber> compressionLinker;

public:
	void CreateElemDiscs(const BiotSubsetParameters &param,
						 const BiotDiscConfig &config)
	{
		CreateElemDiscs(param, config.m_uCmp.c_str(), config.m_uOrder,
							   config.m_pCmp.c_str(), config.m_pOrder,
							   config.m_bSteadyStateMechanics);
	}
	// This is the main function for constructing from a parameter set.
	void CreateElemDiscs(const BiotSubsetParameters &param,
							const char *ucmps, int uorder,
							const char *pcmp, int porder,
							bool bSteadyStateMechanics=true)
		{

			// Create main objects.
			flowEqDisc = make_sp(new TConvectionDiffusion(pcmp, param.get_subsets().c_str()));
			displacementEqDisc = make_sp(new TSmallStrainMechanics(ucmps, param.get_subsets().c_str()));

			 // Do not scale with tau?
			displacementEqDisc->set_stationary(bSteadyStateMechanics);

			// A) Specify displacement eq (for u)
			// (Note: This corresponds to plane strain in 2D)
			auto matLaw = make_sp(new SmallStrainMechanics::HookeLaw<TDomain>());
			matLaw->set_hooke_elasticity_tensor(param.get_lambda(), param.get_mu());

			displacementEqDisc->set_material_law(matLaw);
			displacementEqDisc->set_mass_scale(0.0);

			// Add divergence.
			divLinker = make_sp(new TScaleAddLinkerNumber());
			divLinker->add(param.get_alpha(), flowEqDisc->value());
			displacementEqDisc->set_div_factor(divLinker);

			// B) Specify flow eq (for p).
			compressionLinker = make_sp(new TScaleAddLinkerNumber());
			compressionLinker->add(param.get_alpha(), displacementEqDisc->divergence());

			flowEqDisc->set_mass(compressionLinker);
			flowEqDisc->set_mass_scale(param.get_phi()); 	 //  Storativity 1.0/M = S ‰
			flowEqDisc->set_diffusion((number) param.get_kappa());

			// TODO: Adjust quadrature order
			if (dim==2) {

				if (porder==1) { flowEqDisc->set_quad_order(4); }
				if (uorder==2) { displacementEqDisc->set_quad_order(4); }

			}
			else if (dim == 3) {

				if (uorder == 1) {
				    // displacementEqDisc->set_quad_order(2);
				} else if (uorder == 2) {
				    displacementEqDisc->set_quad_order(5);
				    flowEqDisc->set_quad_order(3);
				}
			}

			// Print info.
			// UG_LOG(flowEqDisc->config_string());
			UG_LOG(displacementEqDisc->config_string());


		}



};

//! This class generates element discretizations
template <typename TDomain>
class BiotElemDiscFactory{
public:
	static const int dim = TDomain::dim;
	typedef BiotElemDisc<TDomain> TBiotDisc;

	BiotElemDiscFactory(const char *ucmps, int uorder, const char *pcmp, int porder,
			bool bSteadyStateMechanics=true)
	: m_config(ucmps, uorder, pcmp, porder, bSteadyStateMechanics){}

	BiotElemDiscFactory(const BiotDiscConfig& config) : m_config(config) {}

	// Create new disc container.
	SmartPtr<TBiotDisc> create_elem_discs(const BiotSubsetParameters &param) const
	{
		SmartPtr<TBiotDisc>  biot = make_sp(new TBiotDisc ());
		biot->CreateElemDiscs(param, m_config);
		return biot;
	}
protected:
	// BiotDiscConfig& config() {return m_config;}
	const BiotDiscConfig& config() {return m_config;} const
	BiotDiscConfig m_config;

};

/* JSON object representation
{
	"discretization" = {
		"ucmp" = "ux, uy, uz"
		"uorder" = 2,

		"pcmp" = "p"
		"porder" = 1,

		"stab" = 1.0
	}


	"gridname" = "mygrid"


	"params" : [
			{ "subset" = "subset1", "alpha" = 1.0, "mu" = 	},
			{ "subset" = "subset2", "alpha" = 1.0	},
	]
}
*/

//! A Biot problem consists of several element discs plus boundary conditions.
/*! For each problem, we can add multiple elemDiscs and boundary conditions. */
template <typename TDomain, typename TAlgebra>
class BiotProblem
{
public:
	static const int dim = TDomain::dim;
	typedef IElemDisc<TDomain> TElemDisc;
	typedef DomainDiscretization<TDomain, TAlgebra> TDomainDisc;
	typedef GridFunction<TDomain,TAlgebra> TGridFunction;
protected:
	typedef ConvectionDiffusionPlugin::ConvectionDiffusionStabFE<TDomain> TConvectionDiffusionStab;
	typedef ConvectionDiffusionPlugin::ConvectionDiffusionFV1<TDomain> TConvectionDiffusionFV1;
public:

	/// CTOR (default orders)
	BiotProblem(const char* uCmp, const char *pCmp, const char *gridname)
	: m_config(uCmp, pCmp), m_gridname(gridname) {}

	/// CTOR (full)
	BiotProblem(const char* uCmp, int uorder, const char *pCmp, int porder, const char *gridname)
	: m_config(uCmp, uorder, pCmp, porder), m_gridname(gridname)  {}

	/// DTOR
	virtual ~BiotProblem() {}

	// get grid name
	const char* get_gridname() const { return m_gridname.c_str();}

	/// Add subset parameters.
	void add_subset_parameters(const BiotSubsetParameters &p)
	{ m_params.push_back(p); }

#ifdef WITH_JSON
	/// Allows adding descriptions.
	void add_subset_parameters(const char* &json_string)
	{
		std::stringstream ss;
		ss << json_string;

		nlohmann::json json;
		ss >> json;

		BiotSubsetParameters _p = json; // implicit use of from_json!
		m_params.push_back(_p);
	}
#endif



	/// Get characteristic time.
	virtual double get_char_time()
	{
		std::vector<BiotSubsetParameters>::const_iterator it= m_params.begin();
		double tchar = 0.0;

		for(it = m_params.begin(); it != m_params.end(); it++)  {
			tchar=std::max(tchar, DefaultCharTime(*it));
		}

		return tchar;
	}

	virtual double start_time()
	{ return 0.0; }

	virtual double end_time()
	{ return 2.0*get_char_time(); }


	/// Adding all elem discs to domain disc.
	virtual void add_elem_discs(SmartPtr<TDomainDisc> dd, bool bSteadyStateMechanics=true)
	{
		// Using factory to create elem discs.
		typedef BiotElemDiscFactory<TDomain> TBiotElemDiscFactory;
		typedef BiotElemDisc<TDomain> TBiotElemDisc;


		BiotDiscConfig conf = config();
		conf.m_bSteadyStateMechanics = bSteadyStateMechanics;

		BiotElemDiscFactory<TDomain> elemDiscFactory(conf);

		// Iterate over parameter sets.
		std::vector<BiotSubsetParameters>::iterator it;
		for(it = m_params.begin(); it != m_params.end(); it++)    {

			// Create elem discs and add
			SmartPtr<TBiotElemDisc> biot = elemDiscFactory.create_elem_discs(*it);

			dd->add(biot->pressure_disc().template cast_dynamic<TElemDisc>());
			dd->add(biot->displacement_disc().template cast_dynamic<TElemDisc>());

		}
	}

protected:
	//! Add stabilizationto domain disc.
	virtual void add_stab_discs(SmartPtr<TDomainDisc> dd, bool bSteadyStateMechanics=true)
	{
		UG_ASSERT(bSteadyStateMechanics == true, "ERROR: Only implemented for Mass matrix!");

		if (config().m_dStab <= 0.0) return;

		SmartPtr<TElemDisc> pDiscStab;
		std::vector<BiotSubsetParameters>::iterator it;
		for(it = m_params.begin(); it != m_params.end(); it++) {
			double gamma = it->get_lambda() + 2.0*it->get_mu();
		  	pDiscStab = make_sp(new TConvectionDiffusionStab(config().m_pCmp.c_str(), it->get_subsets().c_str(), config().m_dStab/gamma));
		  	dd->add(pDiscStab.template cast_dynamic<TElemDisc>());
		}
	}

public:
	//! Add stabilization to domain disc.
	virtual void add_uzawa_discs(SmartPtr<TDomainDisc> dd, bool bSteadyStateMechanics=true)
	{

		UG_ASSERT(bSteadyStateMechanics == true, "ERROR: Only implemented for Mass matrix!");

		SmartPtr<TConvectionDiffusionFV1> pDiscUzawa;
		std::vector<BiotSubsetParameters>::iterator it;
		for(it = m_params.begin(); it != m_params.end(); it++) {
			pDiscUzawa = make_sp(new TConvectionDiffusionFV1(config().m_pCmp.c_str(), it->get_subsets().c_str()));
			pDiscUzawa->set_mass_scale(it->get_beta());
			dd->add(pDiscUzawa.template cast_static<TElemDisc>());
	    }

	}

	//! This add all boundary conditions.
	virtual void add_boundary_conditions(SmartPtr<TDomainDisc> dd, bool bSteadyStateMechanics=true)
	{}

	//! Initial values
	virtual void interpolate_start_values(SmartPtr<TGridFunction> u, double t0)
	{}

	//! Post-processing (per time step)
	virtual bool post_processing(SmartPtr<TGridFunction> u, size_t step, double time) {return true; }



	const BiotDiscConfig& config() const { return m_config; }
	int get_porder() const {return config().m_pOrder;}
	int get_uorder() const {return config().m_uOrder;}

protected:
	BiotDiscConfig m_config;
	std::vector<BiotSubsetParameters> m_params;
	const std::string m_gridname;
};


#ifdef WITH_JSON
template <typename TDomain, typename TAlgebra>
void to_json(nlohmann::json &j, const BiotProblem<TDomain,TAlgebra> &p);

template <typename TDomain, typename TAlgebra>
void from_json(const nlohmann::json &j, BiotProblem<TDomain,TAlgebra> &p);
#endif


} // namespace Poroelasticity
} // namespace ug

#endif
