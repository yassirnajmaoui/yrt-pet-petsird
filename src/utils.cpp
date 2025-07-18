//
// Created by yn257 on 7/17/25.
//

#include "utils.hpp"
#include "DetectorCorrespondenceMap.hpp"
#include "petsird_helpers/geometry.h"

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

::petsird::Coordinate
    yrt::pet::petsird::getCentroid(const ::petsird::BoxShape& box)
{
	const size_t numCorners = box.corners.size();

	::petsird::Coordinate centroid{};
	for (size_t i = 0; i < numCorners; ++i)
	{
		for (int dim = 0; dim < 3; dim++)
		{
			centroid.c[dim] += box.corners[i].c[dim];
		}
	}
	for (int dim = 0; dim < 3; dim++)
	{
		// Divide the sum by number of corners to get the average
		centroid.c[dim] = centroid.c[dim] / static_cast<float>(numCorners);
	}
	return centroid;
}

petsird::Coordinate yrt::pet::petsird::transforms_coord(
    const ::petsird::RigidTransformation& transform,
    const ::petsird::Coordinate& coord)
{
	const xt::xarray<float> hom = xt::linalg::dot(
	    ::petsird_helpers::geometry::transform_to_mat44(transform),
	    ::petsird_helpers::geometry::coordinate_to_homogeneous(coord));
	return ::petsird_helpers::geometry::homogeneous_to_coordinate(hom);
}

std::tuple<Scanner, yrt::pet::petsird::DetectorCorrespondenceMap>
    yrt::pet::petsird::toScanner(
        const ::petsird::ScannerInformation& scannerInfo)
{
	struct IndexedPoint
	{
		Vector3D point;
		Vector3D orientation;
		DetectorCorrespondenceMap::DetectorKey originalKey;
	};

	DetectorCorrespondenceMap correspondenceMap;

	const ::petsird::ScannerGeometry& scannerGeom =
	    scannerInfo.scanner_geometry;

	// Get detecting elements
	::petsird::TypeOfModule numTypeOfModules =
	    scannerGeom.replicated_modules.size();
	size_t totalNumDets = 0;
	for (::petsird::TypeOfModule typeOfModule_i = 0;
	     typeOfModule_i < numTypeOfModules; typeOfModule_i++)
	{
		totalNumDets +=
		    petsird_helpers::get_num_det_els(scannerInfo, typeOfModule_i);
	}

	// PETSIRD "dimensions": type, module, detector
	// Reshuffled to: (DOI), ring, transaxial (YRT-PET dimensions)

	std::vector<IndexedPoint> indexedPoints;
	indexedPoints.resize(totalNumDets);

	// Crystal properties
	float crystalSize_z, crystalSize_trans, crystalDepth;

	// Flat index for the detectors in the flat unshuffled LUT
	det_id_t detId = 0;

	// Properties to measure
	bool isFirstCrystal = true;
	float minZ{};
	float maxZ{};
	float maxDistanceCenter{};

	for (::petsird::TypeOfModule typeOfModule_i = 0;
	     typeOfModule_i < numTypeOfModules; typeOfModule_i++)
	{
		// Number of modules of this type
		const auto& replicatedModuleType =
		    scannerGeom.replicated_modules[typeOfModule_i];
		size_t numModules = replicatedModuleType.NumberOfObjects();

		// Number of detectors in the modules of this type
		const auto& detectors = replicatedModuleType.object.detecting_elements;
		size_t numDetectors = detectors.NumberOfObjects();

		// Get the centroid of the detectors of this module type (reference)
		const auto& currentDetectorShape = detectors.object.shape;
		auto centroid = getCentroid(currentDetectorShape);

		Vector3D crystalOrientation_yrt;
		std::tie(crystalSize_z, crystalSize_trans, crystalDepth,
		         crystalOrientation_yrt) = getCrystalInfo(currentDetectorShape);
		::petsird::Coordinate crystalOrientation{{crystalOrientation_yrt.x,
		                                          crystalOrientation_yrt.y,
		                                          crystalOrientation_yrt.z}};

		for (uint32_t module_i = 0; module_i < numModules; module_i++)
		{
			auto moduleTransform = replicatedModuleType.transforms[module_i];

			// Get the detector orientations (We assume here that the detectors
			//  in a module are parallel)

			for (uint32_t detector_i = 0; detector_i < numDetectors;
			     detector_i++)
			{
				auto detectorTransform = detectors.transforms[detector_i];

				// Transform crystal centroid
				auto totalTransform =
				    petsird_helpers::geometry::mult_transforms(
				        {moduleTransform, detectorTransform});
				auto transformedDetectorCentroid =
				    transforms_coord(totalTransform, centroid);

				const float crystalX = transformedDetectorCentroid.c[0];
				const float crystalY = transformedDetectorCentroid.c[1];
				const float crystalZ = transformedDetectorCentroid.c[2];

				// Assign position to the detCoord
				indexedPoints[detId].point.x = crystalX;
				indexedPoints[detId].point.y = crystalY;
				indexedPoints[detId].point.z = crystalZ;

				// Remove rotation from the transform
				auto translationFreeMatrix = totalTransform.matrix;
				translationFreeMatrix[3] = 0;
				translationFreeMatrix[7] = 0;
				translationFreeMatrix[11] = 0;
				::petsird::RigidTransformation crystalRotation{
				    translationFreeMatrix};

				// Assign orientation to the detCoord
				auto rotatedCrystalOrientation =
				    transforms_coord(crystalRotation, crystalOrientation);
				indexedPoints[detId].orientation.x =
				    rotatedCrystalOrientation.c[0];
				indexedPoints[detId].orientation.y =
				    rotatedCrystalOrientation.c[1];
				indexedPoints[detId].orientation.z =
				    rotatedCrystalOrientation.c[2];

				float distanceCenter = std::sqrt(std::pow(crystalX, 2.0f) +
				                                 std::pow(crystalY, 2.0f) +
				                                 std::pow(crystalZ, 2.0f));

				if (isFirstCrystal)
				{
					minZ = crystalZ;
					maxZ = crystalZ;
					maxDistanceCenter = distanceCenter;
					isFirstCrystal = false;
				}
				else
				{
					minZ = std::min(crystalZ, minZ);
					maxZ = std::max(crystalZ, maxZ);
					maxDistanceCenter =
					    std::max(maxDistanceCenter, distanceCenter);
				}

				// Add the properties for the correspondence afterwards
				indexedPoints[detId].originalKey =
				    DetectorCorrespondenceMap::DetectorKey{
				        typeOfModule_i, module_i, detector_i};

				detId++;
			}
		}
	}

	float axialFOV = maxZ - minZ;

	// Step 1: Sort by Z
	std::sort(indexedPoints.begin(), indexedPoints.end(),
	          [](const IndexedPoint& a, const IndexedPoint& b)
	          { return a.point.z < b.point.z; });

	// Step 2: Group points into axial rings
	std::vector<std::vector<IndexedPoint>> rings;
	std::vector<IndexedPoint> currentRing;
	float lastZ = indexedPoints[0].point.z;

	for (const auto& pt : indexedPoints)
	{
		const float z = pt.point.z;
		if (std::abs(z - lastZ) > EPSILON)
		{
			rings.push_back(currentRing);
			currentRing.clear();
			lastZ = z;
		}
		currentRing.push_back(pt);
	}
	if (!currentRing.empty())
	{
		rings.push_back(currentRing);
	}

	// Check that all rings have the same number of detectors
	bool isFirstRing = true;
	size_t numRings = rings.size();
	size_t detsPerRing{};
	for (size_t ring_i = 0; ring_i < numRings; ring_i++)
	{
		if (isFirstRing)
		{
			detsPerRing = rings[ring_i].size();
			isFirstRing = false;
		}
		else
		{
			if (rings[ring_i].size() != detsPerRing)
			{
				throw std::runtime_error(
				    "Not all rings have the name number of detectors: ring_i=" +
				    std::to_string(ring_i));
			}
		}
	}

	// Sort detectors in each ring transaxially
	for (auto& ring : rings)
	{
		std::sort(ring.begin(), ring.end(),
		          [](const IndexedPoint& a, const IndexedPoint& b)
		          {
			          double thetaA = std::atan2(a.point.y, a.point.x);
			          double thetaB = std::atan2(b.point.y, b.point.x);
			          return thetaA < thetaB;
		          });
	}

	// Return DetCoord
	auto detCoord = std::make_shared<DetCoordOwned>();
	detCoord->allocate(totalNumDets);
	detId = 0;

	for (size_t ring_i = 0; ring_i < numRings; ring_i++)
	{
		detsPerRing = rings[ring_i].size();
		for (size_t det_i = 0; det_i < detsPerRing; det_i++)
		{
			const auto& indexedPoint = rings[ring_i][det_i];

			// Add to DetCoord object
			detCoord->setXpos(detId, indexedPoint.point.x);
			detCoord->setYpos(detId, indexedPoint.point.y);
			detCoord->setZpos(detId, indexedPoint.point.z);
			detCoord->setXorient(detId, indexedPoint.orientation.x);
			detCoord->setYorient(detId, indexedPoint.orientation.y);
			detCoord->setZorient(detId, indexedPoint.orientation.z);

			// Add to correspondence
			correspondenceMap.addMapping(indexedPoint.originalKey.type,
			                             indexedPoint.originalKey.module,
			                             indexedPoint.originalKey.det, detId);

			detId++;
		}
	}

	// Get Scanner properties
	// How to get the minimum angle difference and maximum ring difference ?
	Scanner scanner{scannerInfo.model_name,
	                axialFOV,
	                crystalSize_z,
	                crystalSize_trans,
	                crystalDepth,
	                maxDistanceCenter,
	                detsPerRing,
	                numRings,
	                /*Placeholder: */ 1,
	                /*Placeholder: */ numRings - 1,
	                /*Placeholder: */ 1,
	                1};
	scanner.setDetectorSetup(detCoord);

	return {scanner, correspondenceMap};
}

std::tuple<float, float, float, Vector3D>
    yrt::pet::petsird::getCrystalInfo(const ::petsird::BoxShape& box)
{
	// Get depth dimension
	const auto vertices = box.corners;

	bool isFirstCheck = true;
	float largestDistance{};
	Vector3D largestDistanceUnitVector{};
	std::tuple<int, int> largestDistance_ij;

	for (size_t i = 0; i < vertices.size(); ++i)
	{
		for (size_t j = i + 1; j < vertices.size(); ++j)
		{
			const auto& a = vertices[i];
			const auto& b = vertices[j];

			int same = 0;
			if (std::abs(a.c[0] - b.c[0]) < EPSILON)
			{
				same++;
			}
			if (std::abs(a.c[1] - b.c[1]) < EPSILON)
			{
				same++;
			}
			if (std::abs(a.c[2] - b.c[2]) < EPSILON)
			{
				same++;
			}

			// Only consider pairs that differ in exactly one dimension (i.e.,
			//  share an edge)
			if (same != 2)
			{
				continue;
			}

			Vector3D edgeVector = {b.c[0] - a.c[0], b.c[1] - a.c[1],
			                       b.c[2] - a.c[2]};
			const float distance = std::sqrt(edgeVector.x * edgeVector.x +
			                                 edgeVector.y * edgeVector.y +
			                                 edgeVector.z * edgeVector.z);

			if (isFirstCheck)
			{
				largestDistance = distance;
				largestDistance_ij = {i, j};
				largestDistanceUnitVector = edgeVector.getNormalized();
				isFirstCheck = false;
			}
			else
			{
				if (distance > largestDistance)
				{
					largestDistance = distance;
					largestDistance_ij = {i, j};
					largestDistanceUnitVector = edgeVector.getNormalized();
				}
			}
		}
	}

	if (largestDistanceUnitVector.z > EPSILON)
	{
		std::cerr << "Warning: The crystals given have an orientation with a Z "
		             "component: ["
		          << largestDistanceUnitVector.x << ", "
		          << largestDistanceUnitVector.y << ", "
		          << largestDistanceUnitVector.z << "]" << std::endl;
	}

	// Find the size in the Z dimension
	isFirstCheck = true;
	float zMin{}, zMax{};
	for (const auto& p : vertices)
	{
		if (isFirstCheck)
		{
			zMin = p.c[2];
			zMax = p.c[2];
			isFirstCheck = false;
		}
		else
		{
			zMin = std::min(zMin, p.c[2]);
			zMax = std::max(zMax, p.c[2]);
		}
	}
	float crystalSize_z = zMax - zMin;

	// Get the transaxial size (by first getting the transaxial direction)

	// This assumes the Z direction is not the largest unit vector orientation
	constexpr Vector3D zDir{0, 0, 1};

	const Vector3D thirdDir =
	    zDir.crossProduct(largestDistanceUnitVector).getNormalized();

	float minProj{}, maxProj{};
	isFirstCheck = true;

	for (const auto& vertex : vertices)
	{
		Vector3D v{vertex.c[0], vertex.c[1], vertex.c[2]};
		float proj = v.scalProd(thirdDir);

		if (isFirstCheck)
		{
			minProj = proj;
			maxProj = proj;
			isFirstCheck = false;
		}
		else
		{
			minProj = std::min(minProj, proj);
			maxProj = std::max(maxProj, proj);
		}
	}
	float crystalSize_trans = maxProj - minProj;

	return {crystalSize_z, crystalSize_trans, largestDistance,
	        largestDistanceUnitVector};
}

std::array<petsird::ExpandedDetectionBin, 2>
    petsird_helpers::expand_detection_bin_pair(
        const ScannerInformation& scanner,
        const std::array<TypeOfModule, 2>& type_of_module_pair,
        const std::array<DetectionBin, 2>& detection_bin_pair)
{

	assert(type_of_module < scanner.scanner_geometry.replicated_modules.size());

	std::array<ExpandedDetectionBin, 2> result;

	for (int det_i = 0; det_i < 2; det_i++)
	{
		const auto& rep_module =
		    scanner.scanner_geometry
		        .replicated_modules[type_of_module_pair[det_i]];
		const auto& energy_bin_edges =
		    scanner.event_energy_bin_edges[type_of_module_pair[det_i]];

		const uint32_t num_en = energy_bin_edges.NumberOfBins();
		const uint32_t num_el_per_module =
		    rep_module.object.detecting_elements.transforms.size();

		const auto& bin = detection_bin_pair[det_i];
		const auto det = bin / num_en;
		result[det_i] = {det / num_el_per_module, det % num_el_per_module,
		                 bin % num_en};
	}
	return result;
}
