
#include "datastruct/scanner/DetCoord.hpp"
#include "types.h"

#pragma once

namespace yrt::pet::petsird
{
	// tolerance of 0.1 micron for these calculations:
	constexpr float EPSILON = 1e-4;

	// Apply transformation to a coordinate
	::petsird::Coordinate
	    transforms_coord(const ::petsird::RigidTransformation& transform,
	                     const ::petsird::Coordinate& coord);
	::petsird::Coordinate getCentroid(const ::petsird::BoxShape& box);

	// Returns:
	// - The DetCoord object for the scanner
	// - The original detector indices after the reshuffling
	// - Detectors per ring
	// - Number of rings
	std::tuple<std::shared_ptr<DetCoord>, std::vector<size_t>, size_t, size_t>
	    toDetCoord(std::vector<Vector3D>& points,
	               std::vector<Vector3D>& orientations);
	// Function that computes crystal depth, size in Z, and size in
	//  transaxial and gives you the orientation of the crystal
	// Returns:
	// - Crystal size in Z
	// - Crystal size in transaxial
	// - Crystal depth
	// - Orientation unit vector
	std::tuple<float, float, float, Vector3D>
	    getCrystalInfo(const ::petsird::BoxShape& box);
}  // namespace yrt::pet::petsird
