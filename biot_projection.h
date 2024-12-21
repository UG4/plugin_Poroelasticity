/*
 * Copyright (c) 2019-2022:  G-CSC, Goethe University Frankfurt
 * Copyright (c) 2022-2025:  MSQC, Goethe University Frankfurt
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
#pragma once

#ifndef __BIOT_PROJECTION_H__
#define __BIOT_PROJECTION_H__

// For 

// For 
#include "lib_algebra/operator/interface/linear_operator_inverse.h" // Solver


#include "lib_disc/time_disc/theta_time_step.h" 	// ThetaTimeStep
#include "../Limex/time_disc/time_integrator.hpp" 	// ConstStepLinearTimeIntegrator

#ifdef UG_JSON
#include <nlohmann/json.hpp>
#endif

namespace ug {
namespace Poroelasticity {



//! Creates consistent initial values.
/*!	Solves Au = gradp (for given p).
*/ 
template <typename TDomain, typename TAlgebra>
class BiotProjection
{
public:	
	typedef BiotProblem<TDomain, TAlgebra> TBiotProblem;
	typedef typename TBiotProblem::TDomainDisc TDomainDisc;

	using TSolver = ILinearOperatorInverse<typename TAlgebra::vector_type>;
	using TGridFunction = typename TBiotProblem::TGridFunction;
protected:
	using TTimeIntegrator = ConstStepLinearTimeIntegrator<TDomain,TAlgebra>;

public:	
	//! Constructor (TODO: construct from BiotProblem???)
	BiotProjection(SmartPtr<TDomainDisc> dd0, SmartPtr<TSolver> solver) 
	: m_solver(solver), m_dd0(dd0), time_integrator(SPNULL) {}

	//! Computes a consistent solution 
	void apply (SmartPtr<TGridFunction> sol)
	{
		// Create time integrator (if required)
		if (time_integrator.invalid())
		{
			using TTimeDisc = ThetaTimeStep<TAlgebra>;
			SmartPtr<TTimeDisc> time_disc = make_sp<TTimeDisc> (new TTimeDisc(m_dd0));
			time_integrator = make_sp<TTimeIntegrator> (new TTimeIntegrator(time_disc, m_solver));
		}

		// Perform one step.
		double tau0 = 1.0; // dummy tau
		time_integrator->set_time_step(tau0);
		time_integrator->set_precision_bound(1e-12);
		time_integrator->apply(sol, tau0, sol, 0.0);

	}

protected:

	SmartPtr<TSolver> m_solver; 
	SmartPtr<TDomainDisc> m_dd0;			// domain disc for initial values
	SmartPtr<TTimeIntegrator> time_integrator;
};


} // namespace Poroelasticity
} // namespace ug

#endif
