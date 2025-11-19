
#pragma once

#include "DetectorCorrespondenceMap.hpp"
#include "petsird/types.h"
#include "petsird_helpers.h"
#include "yrt-pet/datastruct/scanner/Scanner.hpp"

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
}  // namespace yrt::petsird

namespace petsird_helpers
{
	std::array<ExpandedDetectionBin, 2> expand_detection_bin_pair(
	    const ScannerInformation& scanner,
	    const std::array<TypeOfModule, 2>& type_of_module_pair,
	    const std::array<DetectionBin, 2>& detection_bin_pair);


	float get_detection_efficiency_from_pair(
	    const ScannerInformation& scanner,
	    const TypeOfModulePair& type_of_module_pair,
	    const std::array<DetectionBin, 2>& module_index_pair,
	    const std::array<uint32_t, 2>& element_index_pair);

}  // namespace petsird_helpers
