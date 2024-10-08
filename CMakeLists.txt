# Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
# Author: Sebastian Reiter
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

################################################################################
# SAMPLE
#
# This script shows how easy it is to add a build script for your plugins.
# Simply specify the plugins name and a list of your source files.
################################################################################

# Name of your plugin and sources.
set(pluginName	Poroelasticity)
set(SOURCES	biot_tools.cpp barry_mercer_data.cpp biot_plugin.cpp)


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.6)



# JSON should be found by ugcore
if (NOT ${nlohmann_json_FOUND})
	find_package(nlohmann_json QUIET 3.2.0)
endif (NOT ${nlohmann_json_FOUND})

message(STATUS "JSON found=" ${nlohmann_json_FOUND} ${nlohmann_json_INCLUDES})

# Enable C++-14 (for Boost on Intel)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


project(UG_PLUGIN_${pluginName})


# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	#target_include_directories(${pluginName} )
	# and additional libraries
	add_library(${pluginName} SHARED ${SOURCES})
	link_directories(${UG_ROOT_PATH}/bin/plugins)
	target_link_libraries (${pluginName} PRIVATE ug4 ConvectionDiffusion SmallStrainMechanics )
endif(buildEmbeddedPlugins)
