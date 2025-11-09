
#pragma once

#include "DetectorCorrespondenceMap.hpp"
#include "yrt-pet/datastruct/scanner/Scanner.hpp"
#include "petsird_helpers.h"
#include "petsird/types.h"

namespace yrt::petsird
{
	// tolerance of 0.1 micron for these calculations:
	constexpr float EPSILON = 1e-4;
	using TimeBlockCollection = std::vector<::petsird::TimeBlock>;

	// Apply transformation to a coordinate
	::petsird::Coordinate
	    transforms_coord(const ::petsird::RigidTransformation& transform,
	                     const ::petsird::Coordinate& coord);
	::petsird::Coordinate getCentroid(const ::petsird::BoxShape& box);

	// Returns:
	// - The Scanner
	// - The original detector indices after the reshuffling
	std::tuple<Scanner, DetectorCorrespondenceMap>
	    toScanner(const ::petsird::ScannerInformation& scannerInfo);
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

namespace petsird_helpers
{
	std::array<ExpandedDetectionBin, 2> expand_detection_bin_pair(
	    const ScannerInformation& scanner,
	    const std::array<TypeOfModule, 2>& type_of_module_pair,
	    const std::array<DetectionBin, 2>& detection_bin_pair);
}  // namespace petsird_helpers
