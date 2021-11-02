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



#include <string>

#include "bridge/util.h"

// replace this with util_domain_dependent.h or util_algebra_dependent.h
// to speed up compilation time
#include "bridge/util_domain_algebra_dependent.h"

#include "common/ug_config.h"
#include "common/error.h"

// Plugin stuff
#include "biot_tools.h"
#include "barry_mercer.h"

using namespace std;
using namespace ug::bridge;

#ifdef UG_PARALLEL
#include "pcl/pcl_util.h"
#endif

namespace ug{
namespace Poroelasticity{

/** 
 *  \defgroup Poroelasticity Poroelasticity
 *  \ingroup plugins
 *  This is a small sample plugin.
 *  \{
 */

void PluginSaysHello()
{
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
	cout << "Hello, I'm your plugin on proc " <<
				pcl::ProcRank() << "." << endl;
	pcl::SynchronizeProcesses();
#else
	UG_LOG("Hello, I'm your personal plugin in serial environment!\n");
#endif
}

void CrashFct(const string& reason)
{
	UG_THROW("I Crash because: "<< reason);
}

void CrashFctFatal(const string& reason)
{
	UG_THROW("I Crash fatal because: "<< reason);
}

void PluginCrashes()
{
	try{
		CrashFct("Some funny reason");
	}
	catch(bad_alloc& err)
	{
		UG_LOG("Bad Alloc");
	}
	UG_CATCH_THROW("Something wrong in PluginCrashes");
}

void PluginCrashesFatal()
{
	try{
		CrashFctFatal("Some fatal reason");
	}
	catch(bad_alloc& err)
	{
		UG_LOG("Bad Alloc");
	}
	UG_CATCH_THROW("Something wrong in PluginCrashesFatal");
}

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts
 * of the plugin. All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
//	useful defines
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	{
		typedef BiotProblem<TDomain, TAlgebra> T;

		string name = string("BiotProblem").append(suffix);
		reg.add_class_<T>(name, grp)
		   .template add_constructor<void (*)(const char*,const char*,const char*)>("ucmp(s)#pcmp(s)")
		   .template add_constructor<void (*)(const char*,int,const char*,int,const char*)>("ucmp(s)#uorder#pcmp(s)#porder")
		   .add_method("get_uorder", &T::get_uorder)
		   .add_method("get_porder", &T::get_porder)
		   .add_method("get_gridname", &T::get_gridname)
		   .add_method("get_char_time", &T::get_char_time)
		   .add_method("start_time", &T::start_time)
		   .add_method("end_time", &T::end_time)
		   .add_method("add_elem_discs", &T::add_elem_discs)
		   //.add_method("add_stab_discs", &T::add_stab_discs)
		   .add_method("add_uzawa_discs", &T::add_uzawa_discs)
		   .add_method("interpolate_start_values", &T::interpolate_start_values)
		   .add_method("post_processing", &T::post_processing)
		 //  .add_method("add_uzawa_discs", &T::add_uzawa_discs)
		   .add_method("add_boundary_conditions", static_cast<void (T::*)(SmartPtr<typename T::TDomainDisc>, bool)> (&T::add_boundary_conditions))
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BiotProblem", tag);

	}



}

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
//	useful defines
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

	{
			typedef BiotElemDisc<TDomain> T;

			string name = string("BiotElemDisc").append(suffix);
			reg.add_class_<T>(name, grp)
			   .add_constructor()
			   .add_method("pressure_disc", &T::pressure_disc)
			   .add_method("displacement_disc", &T::displacement_disc)
			   .add_method("compression_linker", &T::compression_linker)
			   .add_method("divergence", &T::divergence)

			   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "BiotElemDisc", tag);

		}
	{
		typedef BiotElemDiscFactory<TDomain> T;
		typedef IElemDisc<TDomain> TElemDisc;

		string name = string("BiotElemDiscFactory").append(suffix);
		reg.add_class_<T>(name, grp)
			.template add_constructor<void (*)(const char*, int,
					const char*, int, bool)>("")
			.add_method("create_elem_discs", &T::create_elem_discs)
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "BiotElemDiscFactory", tag);

	}

}

/**
 * Function called for the registration of Dimension dependent parts
 * of the plugin. All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
//	useful defines
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts
 * of the plugin. All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
//	useful defines
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

}

/**
 * Function called for the registration of Domain and Algebra independent parts
 * of the plugin. All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{

	{
		typedef BiotSubsetParameters T;

		string name = string("BiotSubsetParameters");
		reg.add_class_<T>(name, grp)
		    .add_constructor()
			.add_constructor<void (*)(const char*,number, number, number, number, number, number)>("Subset(s)#...");
			// .set_construct_as_smart_pointer(true);

	}

	{
		// Bessel functions.
		reg.add_function("BesselJ0", &BesselJ0, grp);
		reg.add_function("BesselJ1", &BesselJ1, grp);
	}
}

}; // end Functionality



struct FunctionalityFor2D
{
	template <typename TDomain, typename TAlgebra>
	static void DomainAlgebra(Registry& reg, string grp)
	{
	//	useful defines
	const string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	const string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	{
			typedef BiotProblem<TDomain, TAlgebra> TBase;
			typedef BarryMercerProblem<TDomain, TAlgebra> T;

			string name = string("BarryMercerProblem").append(suffix);
			reg.add_class_<T,TBase>(name, grp)
				.template add_constructor<void (*)(const char*,const char*)>("ucmp(s)#pcmp(s)")
			 //  .template add_constructor<void (*)(const char*,int,const char*,int)>("ucmp(s)#uorder#pcmp(s)#porder")
			  // .add_method("add_elem_discs", &T::add_elem_discs)

			 //  .add_method("add_uzawa_discs", &T::add_stab_discs)
			 //  .add_method("add_boundary_conditions", static_cast<void (T::*)(SmartPtr<typename T::TDomainDisc>, bool)> (&T::add_boundary_conditions))
			   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "BarryMercerProblem", tag);

		}
	}
};

// end group sample_plugin
/// \}

}// end of namespace Poroelasticity


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_Poroelasticity(Registry* reg, string grp)
{
	grp.append("/Poroelasticity");
	typedef Poroelasticity::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(*reg,grp);
		// RegisterDimensionDependent<Functionality>(*reg,grp);
		// RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterAlgebraDependent<Functionality>(*reg,grp);
		//RegisterDomainAlgebraDependent<Functionality>(*reg,grp);

		// Register only for 2D/3D (cf. SmallStrainMechanics)
		RegisterDomain2d3dDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomain2dAlgebraDependent<Poroelasticity::FunctionalityFor2D>(*reg,grp);

	}
	UG_REGISTRY_CATCH_THROW(grp);
}

extern "C" UG_API void
FinalizeUGPlugin_Poroelasticity()
{
}


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////


}//	end of namespace ug
