# Copyright 2013-2017 Imperial College London
# Copyright 2013-2017 Andreas Schuh
# Copyright 2018 King's College London
# SPDX-License-Identifier: Apache-2.0
#
# This file has been modified from the original source.

################################################################################
# @file  BasisProject.cmake
# @brief Sets basic information about the MIRTK module and calls basis_project().
#
# This file defines basic information about a project by calling
# the basis_project() function. This basic information, also known as metadata,
# is used by CMake BASIS to setup the project. The dependencies to other modules
# have to be specified here such that the top-level IRTK project can analyze the
# inter-module dependencies, as well as dependencies on third-party libraries.
#
# @sa https://cmake-basis.github.io/standard/modules.html
#
# @ingroup BasisSettings
################################################################################

# Note: The #<*> dependency patterns are required by the basisproject tool and
#       should be kept on a separate line as last commented argument of the
#       corresponding options of the basis_project() command. The TEMPLATE
#       option and set argument are also required by this tool and should not
#       be changed manually. The argument is updated by basisproject --update.

basis_project (

  # ----------------------------------------------------------------------------
  # meta-data
  NAME        "SVR"
  VERSION     "1.0.0"
  SOVERSION   "0"
  PACKAGE     "MIRTK"
  AUTHORS     "Maria Deprez, Naomi Clarke"
  DESCRIPTION "SVR Reconstruction for MRI Medical Imaging."
  COPYRIGHT   "2024 King's College London"
  CONTACT     "Maria Deprez <maria.deprez@kcl.ac.uk>, Naomi Clarke <naomi.5.clarke@kcl.ac.uk>"
  TEMPLATE    "mirtk-module/1.0"

  # ----------------------------------------------------------------------------
  # dependencies
  DEPENDS
    MIRTK{Common,Numerics,Image,IO,Transformation,Registration}
    Boost{filesystem,program_options,timer}
  #<dependency>
  OPTIONAL_DEPENDS
    Python{Interp}
    TBB{tbb,malloc}
    #<optional-dependency>
  TOOLS_DEPENDS
    Python{Interp}
  TEST_DEPENDS
    Boost{unit_test_framework}
)
