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

#include "barry_mercer_data.h"
#include <boost/lexical_cast.hpp>

namespace ug {
namespace Poroelasticity{

const double BarryMercerNondimensional::X0 = 0.25;
const double BarryMercerNondimensional::Y0 = 0.25;
const double BarryMercerNondimensional::m_PI = ug::PI;
const size_t BarryMercerNondimensional::NAPPROX = 512;


//! double t_hat = m_beta*t;

/*! Pressure is normalized and should be multiplied by (BARRY_MERCER_DATA.LAMBDA+2.0*BARRY_MERCER_DATA.MU) */
double BarryMercerNondimensional::Pressure2D(double x, double y, double t_hat) const
{

  const int N = NAPPROX;

  //double x0 = X0;
  //double y0= Y0;

  double sinbt = sin(t_hat);
  double cosbt = cos(t_hat);

  double pp = 0.0;

  const int m = N;
  // ALT: for (int m=2; m<=N ; ++m) {
      for (int k=1; k <= m-1; ++k) {

    	 // const double lambda_k = k*m_PI;
    	  const double sinkx = sin(m_PI*k*x);

    	  for (int q=1; q <= m-1; ++q) {
    	 // ALT:  const int q = m-k;
    	  //const double lambda_q = q*m_PI;
    	  const double sinqy = sin(m_PI*q*y);

    	  const double _coeff_kq = FourierCoeff_P(k,q, t_hat);
    	  //  std::cout << k<< "," << q << ":" << boost::lexical_cast<string>(_coeff_kq) << std::endl;
       //-- print ("x= ("..x..", "..y.."\tn="..n.."q="..q..",m="..m.."), c=".._coeff_nq.."\t"..sinnx.."\t"..sinqy.."=> \tpp="..pp.."\tupdate=".._coeff_nq * sinnx * sinqy)
    	  pp += _coeff_kq * sinkx * sinqy;
    }

   // -- print("break: "..sinnx)
  }

//  -- print ("x= ("..t..","..t_hat..","..x..","..y..")"..pp)
  return -4.0*pp;
}

double BarryMercerNondimensional::VelX2D(double x, double y, double t_hat) const
{

	int N = NAPPROX;
//	double t_hat = BARRY_MERCER_DATA.BETA*t;

	double u = 0.0;

	for (int m=2; m<=N; ++m) {
		for (int k=1; k<=m-1; ++k) {


			//int n=k;
			double _lambda_k = m_PI*k;
			double _coskx = cos(_lambda_k*x)*_lambda_k;

			//-- for q=1,N do
			int q = m-k;
			double _lambda_kq = (m_PI*m_PI)*(k*k+ q*q);
			double _coeff_kq = FourierCoeff_P(k,q, t_hat);
			// -- < print ("x= ("..x..","..y..","..n..","..q.."), coeff="..(_coeff_nq/_lambda_nq)..", u="..u)
			u += _coeff_kq * _coskx * sin(m_PI*q*y)/_lambda_kq;
		}
	}

	return 4.0*u;
}


// local t_hat = m_beta*t;
double BarryMercerNondimensional::VelY2D(double x, double y, double t_hat) const
{

	int N = NAPPROX;
	double u = 0.0;

	for (int m=2; m<=N; ++m) {
		for (int k=1; k<=m-1; ++k) {
				//-- for q=1,N do
			int q = m-k;
			double _cosqx = cos(m_PI*q*y)*m_PI*q;
				//  -- for n=1,N do
			double _lambda_kq = (m_PI*m_PI)*(k*k+ q*q);
        	double _coeff_kq = FourierCoeff_P(k,q, t_hat);
				//  --  print ("x= ("..x..","..y..")"..pp..", c=".._coeff_nq)
			u += _coeff_kq * sin(m_PI*k*x)*_cosqx/_lambda_kq;
		}
	}

	return 4.0*u;
}

}
}





