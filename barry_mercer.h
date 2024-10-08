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

#ifndef POROELASTICITY_BARRY_MERCER_H_
#define POROELASTICITY_BARRY_MERCER_H_

#include "lib_disc/io/vtkoutput.h" // VTKOutput
#include "lib_disc/function_spaces/interpolate.h" // Interpolate

#include "lib_disc/spatial_disc/elem_disc/dirac_source/lagrange_dirac_source.h"
#include "lib_disc/spatial_disc/constraints/dirichlet_boundary/lagrange_dirichlet_boundary.h"

#include "barry_mercer_data.h"
#include "biot_tools.h"

namespace ug {
namespace Poroelasticity{



//! Auxiliary class for compution errors as 'StdGlobPosData'.
template <class TGridFunction>
class BarryMercerErrorData {

	enum normTypes {L2NORM_P=0, L2NORM_UX, L2NORM_UY, H1SEMI_UX, H1SEMI_UY, _SIZE};

	double m_normErr[_SIZE];
	double m_normSol[_SIZE];
	double m_normRef[_SIZE];

protected:
	void ComputeNorms(TGridFunction &uref, double *norms)
	{
		int porder = 2;
		int uorder = 4;

		norms[L2NORM_P] = L2Norm(uref, "p", porder);
		norms[L2NORM_UX] = L2Norm(uref, "ux", uorder);
		norms[L2NORM_UY] = L2Norm(uref, "uy", uorder);

		norms[H1SEMI_UX] = H1SemiNorm(uref, "ux", uorder);
		norms[H1SEMI_UY] = H1SemiNorm(uref, "uy", uorder);
	};

	/*
	function PrintNorms(normDesc)
	for key, val in pairs(normDesc) do  UG_LOG(key.."\t"..val) end
	end

	function CompareNorms(normDesc, errDesc, refDesc)
	for key, val in pairs(normDesc) do  UG_LOG(key.."\t"..val.."\t"..errDesc[key].."\t("..errDesc[key]/refDesc[key]..")") end
	end
	*/
	void PrintNorms(double *norms)
	{
		for (int key=0; key<_SIZE; ++key)
		{ UG_LOG("KEY" << key << " \t" <<norms[key] << std::endl); }
	}

	void CompareNorms(double *normDesc,  double *errDesc, double *refDesc)
	{
			for (int key=0; key<_SIZE; ++key)
			{
				UG_LOG("KEY" << key << " \t"
						<< normDesc[key] << "\t"
						<< errDesc[key]   << "\t ("
						<< (errDesc[key]/refDesc[key]) << ")" << std::endl);
				//UG_LOG(key.."\t"..val.."\t"..errDesc[key].."\t("..errDesc[key]/refDesc[key]..")")
			}
	}
protected:
	SmartPtr<BarryMercerRefPressure::pos_data_type> m_spPressure;
	SmartPtr<BarryMercerRefDispX::pos_data_type> m_spDispX;
	SmartPtr<BarryMercerRefDispY::pos_data_type> m_spDispY;
public:
	void init (BarryMercerData& dimCoeffs)
	{
		m_spPressure = make_sp(new BarryMercerRefPressure(dimCoeffs));
		m_spDispX = make_sp(new BarryMercerRefDispX(dimCoeffs));
		m_spDispY = make_sp(new BarryMercerRefDispY(dimCoeffs));
	}

	void eval(BarryMercerData& dimCoeffs, TGridFunction &u, int step, double time)
	{
	int napprox = BarryMercerNondimensional::NAPPROX;
	const int dim = TGridFunction::dim;

	if (napprox<=0) return;
	 UG_LOG ("NAPPROX =" << napprox << std::endl);

	 // VTK output.
	 typedef VTKOutput<dim> TVTKOutput;
	 TVTKOutput vtk = TVTKOutput();

	 // Aux. subsets
	 const char* mandatory_subsets = "INNER,SINGULARITY,CORNERS,HORIZONTAL,VERTICAL";

	 // Aux. vector.
	 SmartPtr<TGridFunction> uref = u.clone();
	 uref->set(0.0);

	 const double charTime = dimCoeffs.tchar;
	 UG_LOG ("time =" << time << "("<< time/charTime <<")" << std::endl);

	 // Evaluate reference data.
	 this->init(dimCoeffs);
	 Interpolate(m_spPressure, uref, "p", mandatory_subsets, time);
	 Interpolate(m_spDispX, uref, "ux", mandatory_subsets, time);
	 Interpolate(m_spDispY, uref, "uy", mandatory_subsets, time);

	 // Print solution.
	 vtk.print("BarryMercer2D_Sol.vtu", u, step, time);
	 vtk.print("BarryMercer2D_Ref.vtu", *uref, step, time);

	 // Compute norms.
	 ComputeNorms(u, m_normSol);
	 ComputeNorms(*uref, m_normRef);

	 UG_LOG ("REFERENCE:");
	 PrintNorms(m_normRef);

	 // Compute errors.
	 VecScaleAdd((typename TGridFunction::vector_type&) *uref, 1.0, *uref, -1.0, u);
	 vtk.print("BarryMercer2D_Err.vtu", *uref, step, time);
	 ComputeNorms(*uref, m_normErr);

	 UG_LOG ("SOLUTION/ERROR:" << std::endl);
	 CompareNorms(m_normSol, m_normErr, m_normRef);



	 // More output.
	 UG_LOG("deltaP:\t" << time << "\t" << time/charTime  << "\t"
			 << m_normErr[L2NORM_P] << "\t"<< m_normSol[L2NORM_P]<<"\t"<< m_normRef[L2NORM_P]<< std::endl);

	 UG_LOG ("deltaU1A:\t"<<time<<"\t"<<time/charTime<<"\t"
			<< m_normErr[H1SEMI_UX] << "\t"<< m_normSol[H1SEMI_UX]<<"\t"<< m_normRef[H1SEMI_UX]<< std::endl);

	 UG_LOG ("deltaU2A:\t"<<time<<"\t"<<time/charTime<<"\t"
			<< m_normErr[H1SEMI_UY] << "\t"<< m_normSol[H1SEMI_UY]<<"\t"<< m_normRef[H1SEMI_UY] << std::endl);

	 UG_LOG ("deltaU1B:\t"<<time<<"\t"<<time/charTime<<"\t"
			   << m_normErr[L2NORM_UX] << "\t"<< m_normSol[L2NORM_UX]<<"\t"<< m_normRef[L2NORM_UX]<< std::endl);

	 UG_LOG ("deltaU2B:\t"<<time<<"\t"<<time/charTime<<"\t"
			   << m_normErr[L2NORM_UY] << "\t"<< m_normSol[L2NORM_UY]<<"\t"<< m_normRef[L2NORM_UY]<< std::endl);


	}

};







/// Implementation as a Biot problem
template <typename TDomain, typename TAlgebra>
class BarryMercerProblem
: public BiotProblem<TDomain,TAlgebra>
{
public:

	typedef BiotProblem<TDomain,TAlgebra> base_type;
	typedef DomainDiscretization<TDomain,TAlgebra> TDomainDisc;

private:
	typedef DiracSourceDisc<TDomain> TDiracSourceDisc;
	typedef DirichletBoundary<TDomain,TAlgebra> TDirichletBoundary;

public:

	BarryMercerProblem(const char* uCmp, const char *pCmp)
	: base_type(uCmp, pCmp, "../grids/barrymercer2D-tri.ugx"), m_a(1.0), m_b(1.0), nskip(1)
	{
		set_default_parameters();
	}

	BarryMercerProblem(const BiotDiscConfig& config)
	: base_type(config, "../grids/barrymercer2D-tri.ugx"), m_a(1.0), m_b(1.0), nskip(1)
	{
		set_default_parameters();
	}


	virtual ~BarryMercerProblem() {}


protected:

	void set_default_parameters()
	{
		 double E = 1e+5; 		// Young's elasticity modulus [Pa]
		 double nu = 0.4;       // Poisson's ratio  [1]

		 double lambda = (E*nu)/((1.0+nu)*(1.0-2.0*nu));
		 double mu = 0.5*E/(1+nu);

		 double kappa = 1e-5;  	// Permeability [m*m]
		 double muf = 1e-3;     // Pa*s    => Diff Coeff 1e-9
		 double alpha = 1.0;

		 //double Kcomp = E/(3*(1-2*nu)); // compression (or bulk) modulus)
		 double Kv = 2.0*E/(1+nu)*(1.0-nu)/(1.0-2.0*nu)  ; // uni-axial drained bulk modulus

		 double beta_uzawa=(alpha*alpha)/Kv*(2.0-2.0*nu);

		 base_type::m_params.resize(1);
		 base_type::m_params[0]=BiotSubsetParameters("INNER", alpha, kappa/muf, 0.0, lambda, mu, beta_uzawa);
	}

public:



	double start_time() override
	{ return 0.0; }

	double end_time() override
	{ return 2.0*M_PI*base_type::get_char_time(); }

	void add_elem_discs(SmartPtr<TDomainDisc> dd, bool bSteadyStateMechanics=true) override
	{
		// Add default Biot discs.
		base_type::add_elem_discs(dd, bSteadyStateMechanics);

		// Add point source.
		SmartPtr<TDiracSourceDisc> m_pointSourceDisc = make_sp(new TDiracSourceDisc("p", "SINGULARITY"));

	    double beta = base_type::m_params[0].get_kappa()*( base_type::m_params[0].get_lambda() + 2* base_type::m_params[0].get_mu());
	    SmartPtr<BarryMercerPointSource::pos_data_type> m_pointSourceFunc;
	    m_pointSourceFunc = make_sp(new BarryMercerPointSource(beta));

	    MathVector<2> point(0.25, 0.25);
	    m_pointSourceDisc->add_source(m_pointSourceFunc, point);
		dd->add(m_pointSourceDisc.template cast_static<typename TDiracSourceDisc::base_type>());

		// Add default stabilization.
		const BiotDiscConfig& discretization = base_type::config();
		if (discretization.m_uOrder ==  discretization.m_pOrder) {
			base_type::add_stab_discs(dd, bSteadyStateMechanics);
		}
	}


	//! Initial values
	void interpolate_start_values(SmartPtr<typename base_type::TGridFunction> u, double t0) override
	{ u->set(0.0); }

	//! Add all boundary conditions.
	virtual void add_boundary_conditions(SmartPtr<TDomainDisc> dd, bool bSteadyStateMechanics=true) override
	{
		SmartPtr<TDirichletBoundary> m_spDirichlet = make_sp(new TDirichletBoundary(false));
		m_spDirichlet->add(0.0, "p", "VERTICAL,HORIZONTAL,CORNERS");
		m_spDirichlet->add(0.0, "ux", "HORIZONTAL,CORNERS");
		m_spDirichlet->add(0.0, "uy", "VERTICAL,CORNERS");

		dd->add(m_spDirichlet.template cast_static<typename TDirichletBoundary::base_type> ());
	}


	/// Post-processing (per time step)
	bool post_processing(SmartPtr<typename base_type::TGridFunction> u, size_t step, double time) override
	{
		if ((step-1) % nskip !=0 ) return true; // Execute first n

		BarryMercerData dimData( base_type::get_char_time(), base_type::m_params[0].get_lambda(), base_type::m_params[0].get_mu());
		m_errData.eval(dimData, *u, step, time);
		return true;
	}

	void set_skip(size_t skip)
	{ nskip = skip; }

protected:

	/// Inverse of consolidation coefficient.
	double default_beta() const {
		return  base_type::m_params[0].get_kappa()*( base_type::m_params[0].get_lambda() + 2*base_type::m_params[0].get_mu());
	}

protected:
	double m_a;
	double m_b;

	size_t nskip;


	BarryMercerErrorData<typename base_type::TGridFunction> m_errData;

};
}
}

#endif /* POROELASTICITY_BARRY_MERCER_H_ */
