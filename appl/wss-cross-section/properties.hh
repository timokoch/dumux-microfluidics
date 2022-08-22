// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup OnePTests
 * \brief The properties for the incompressible test
 */

#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_PROPERTIES_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/fvspatialparams1pconstant.hh>

#include "problem.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct OnePIncompressible { using InheritsFrom = std::tuple<OneP>; };
struct OnePIncompressibleBox { using InheritsFrom = std::tuple<OnePIncompressible, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePIncompressible>
{ using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePIncompressible> { using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePIncompressible>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FVPorousMediumFlowSpatialParamsOnePConstant<GridGeometry, Scalar>;
};

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePIncompressible> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePIncompressible>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar> >;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };

// define a TypeTag for a quad precision test
#if HAVE_QUADMATH
namespace TTag {
struct OnePIncompressibleTpfaQuad { using InheritsFrom = std::tuple<OnePIncompressibleTpfa>; };
} // end namespace TTag
template<class TypeTag>
struct Scalar<TypeTag, TTag::OnePIncompressibleTpfaQuad> { using type = Dune::Float128; };
#endif
}
#endif
