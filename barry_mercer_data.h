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

#ifndef POROELASTICITY_BARRY_MERCER_DATA_H_
#define POROELASTICITY_BARRY_MERCER_DATA_H_

#include "biot_tools.h"

namespace ug {
namespace Poroelasticity{


//! Non-dimensional solution.
class BarryMercerNondimensional {

public:
	double Pressure2D(double x, double y, double t) const;
	double VelX2D(double x, double y, double t) const;
	double VelY2D(double x, double y, double t) const;
	static const size_t NAPPROX;

protected:
	//! Compute coefficient from Eq. (24) in Barry & Mercer, ACME, 1999 (for $\omega=1).
	inline double FourierCoeff_P(int n, int q, double t_hat) const
	{

	  if ((n%4==0) || (q%4==0)) { return 0.0; }

	  const double lambda_n = n*m_PI;
	  const double lambda_q = q*m_PI;
	  const double _lambda_nq = lambda_n*lambda_n + lambda_q*lambda_q;

	  const double val1 = -2.0 * sin(lambda_n*X0) * sin(lambda_q*Y0);
	  const double val2 = (_lambda_nq*sin(t_hat) - cos(t_hat) + exp(-_lambda_nq*t_hat));
	  const double val3 = (1 + _lambda_nq*_lambda_nq);

	  return (val1*val2)/val3;
	}

	static const double m_PI;
	static const double X0;
	static const double Y0;

};






//! Dimensional coefficients for Barry-Mercer benchmark.
struct BarryMercerData
{
public:
	// BarryMercerData () : a(1.0), b(1.0), tchar(1.0) {}
	BarryMercerData (double tchar_, double lamdda_, double mu_)
	: a(1.0), b(1.0), tchar(tchar_),
	  lambda(lamdda_), mu(mu_) {}

	double a;
	double b;
	double tchar;

	double lambda;
	double mu;
};



//! Evaluate reference pressure.
class BarryMercerRefPressure
: public StdGlobPosData<BarryMercerRefPressure, number, 2, void>
{
public:
	//! Export base type
	typedef StdGlobPosData<BarryMercerRefPressure, number, 2, void> pos_data_type;

	//! CTOR
	BarryMercerRefPressure (const BarryMercerData &dimCoeffs)
	: m_nonDimData(), m_dimData(dimCoeffs) {}

	//! Define eval function.
	inline void evaluate(number& p, const MathVector<2>& x, number time, int si) const
	{
		p = m_nonDimData.Pressure2D(x[0]/m_dimData.a, x[1]/m_dimData.b, time/m_dimData.tchar);
		p *= (m_dimData.lambda + 2.0 *m_dimData.mu);
	}

protected:
	const BarryMercerNondimensional m_nonDimData;
	const BarryMercerData &m_dimData;
};


//! Evaluate reference pressure.
class BarryMercerRefDispX
: public StdGlobPosData<BarryMercerRefDispX, number, 2, void>
{
public:
	//! Export base type
	typedef StdGlobPosData<BarryMercerRefDispX, number, 2, void> pos_data_type;

	//! CTOR
	BarryMercerRefDispX (const BarryMercerData &dimData)
	: m_nonDimData(), m_dimData(dimData) {}

	//! Define eval function.
	inline void evaluate(number& p, const MathVector<2>& x, number time, int si) const
	{ p = m_nonDimData.VelX2D(x[0]/m_dimData.a, x[1]/m_dimData.b, time/m_dimData.tchar); }

protected:
	const BarryMercerNondimensional m_nonDimData;
	const BarryMercerData &m_dimData;
};

//! Evaluate reference pressure.
class BarryMercerRefDispY
: public StdGlobPosData<BarryMercerRefDispY, number, 2, void>
{
public:
	//! Export base type
	typedef StdGlobPosData<BarryMercerRefDispY, number, 2, void> pos_data_type;

	//! CTOR
	BarryMercerRefDispY(const BarryMercerData &dimData)
	: m_nonDimData(), m_dimData(dimData) {}

	//! Define eval function.
	inline void evaluate(number& p, const MathVector<2>& x, number time, int si) const
	{ p = m_nonDimData.VelY2D(x[0]/m_dimData.a, x[1]/m_dimData.b, time/m_dimData.tchar); }

protected:
	const BarryMercerNondimensional m_nonDimData;
	const BarryMercerData &m_dimData;
};



//! This defines a point source as 'StdGlobPosData'
class BarryMercerPointSource
: public StdGlobPosData<BarryMercerPointSource, number, 2, void>
{
public:
	//! Export base type
	typedef StdGlobPosData<BarryMercerPointSource, number, 2, void> pos_data_type;

	//! CTOR
	BarryMercerPointSource (const double consolidation)
	: m_nonDimData(), m_beta(consolidation) {}

	//! Define eval function.
	inline void evaluate(number& val, const MathVector<2>& x, number time, int si) const
	{
		double beta_ = get_beta();
		val = 2.0*beta_*sin(beta_*time);
	}

	inline double get_beta() const
	{return m_beta;}

protected:
	const BarryMercerNondimensional m_nonDimData;
	double m_beta;
};


}
}

#endif /* POROELASTICITY_BARRY_MERCER_H_ */
