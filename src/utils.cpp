//
// Created by yn257 on 7/17/25.
//

#include "utils.hpp"
#include "petsird_helpers/geometry.h"

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

::petsird::Coordinate
    yrt::pet::petsird::getCentroid(const ::petsird::BoxShape& box)
{
	size_t numCorners = box.corners.size();
	::petsird::Coordinate centroid;
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
	xt::xarray<float> hom = xt::linalg::dot(
	    ::petsird_helpers::geometry::transform_to_mat44(transform),
	    ::petsird_helpers::geometry::coordinate_to_homogeneous(coord));
	return ::petsird_helpers::geometry::homogeneous_to_coordinate(hom);
}

std::tuple<std::shared_ptr<DetCoord>, std::vector<size_t>, size_t, size_t>
    yrt::pet::petsird::toDetCoord(std::vector<Vector3D>& points,
                                  std::vector<Vector3D>& orientations)
{
	struct IndexedPoint
	{
		Vector3D point;
		size_t originalIndex;
	};

	std::vector<IndexedPoint> indexedPoints;
	indexedPoints.resize(points.size());
	for (size_t i = 0; i < points.size(); ++i)
	{
		indexedPoints[i] = {points[i], i};
	}

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
		float z = pt.point.z;
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
	size_t detsInRing{};
	for (size_t ring_i = 0; ring_i < numRings; ring_i++)
	{
		if (isFirstRing)
		{
			detsInRing = rings[ring_i].size();
			isFirstRing = false;
		}
		else
		{
			if (rings[ring_i].size() != detsInRing)
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
	detCoord->allocate(points.size());
	std::vector<size_t> originalIndices(points.size());
	size_t detectorId = 0;

	for (size_t ring_i = 0; ring_i < numRings; ring_i++)
	{
		detsInRing = rings[ring_i].size();
		for (size_t det_i = 0; det_i < detsInRing; det_i++)
		{
			detCoord->setXpos(detectorId, rings[ring_i][det_i].point.x);
			detCoord->setYpos(detectorId, rings[ring_i][det_i].point.y);
			detCoord->setZpos(detectorId, rings[ring_i][det_i].point.z);

			size_t originalIndex = rings[ring_i][det_i].originalIndex;
			Vector3D originalOrientation = orientations[originalIndex];
			originalIndices[detectorId] = originalIndex;

			detCoord->setXorient(detectorId, originalOrientation.x);
			detCoord->setYorient(detectorId, originalOrientation.y);
			detCoord->setZorient(detectorId, originalOrientation.z);

			detectorId++;
		}
	}

	return {std::move(detCoord), originalIndices, detsInRing, numRings};
}

std::tuple<float, float, float, Vector3D>
    yrt::pet::petsird::getCrystalInfo(const ::petsird::BoxShape& box)
{
	// Get depth dimension
	auto vertices = box.corners;

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
			float distance = std::sqrt(edgeVector.x * edgeVector.x +
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
	Vector3D zDir{0, 0, 1};
	Vector3D thirdDir =
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
