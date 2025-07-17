#include "datastruct/scanner/Scanner.hpp"
#include "utils.hpp"
#include "utils/Utilities.hpp"

#include "H5Cpp.h"
#include "binary/protocols.h"
#include "hdf5.h"
#include "hdf5/protocols.h"
#include "petsird_helpers/create.h"
#include "petsird_helpers/geometry.h"
#include "protocols.h"
#include "types.h"

#include "CLI11.hpp"
#include <string>


int main(int argc, char** argv)
{
	CLI::App app{"YourApp description"};

	// Variables to hold parsed values
	std::string input_fname;
	int numSubsets = 0;
	int numIterations = 0;
	std::string psfKernel_fname;
	std::string outScannerLUT_fname;
	std::string outScannerJSON_fname;
	std::string outSensImage_fname;
	std::string outImage_fname;

	// Add options
	app.add_option("-i,--input", input_fname, "Input PETSIRD file")
	    ->required()
	    ->check(CLI::ExistingFile);

	app.add_option("--num-subsets", numSubsets, "Number of subsets");

	app.add_option("--num-iterations", numIterations, "Number of iterations");

	app.add_option("--psf", psfKernel_fname, "PSF kernel file")
	    ->check(CLI::ExistingFile);

	app.add_option("--out-scanner-lut", outScannerLUT_fname,
	               "Output Scanner LUT file");
	app.add_option("--out-scanner-json", outScannerJSON_fname,
	               "Output Scanner JSON file");

	app.add_option("--out-sens", outSensImage_fname,
	               "Output sensitivity image file");
	app.add_option("--out-recon", outImage_fname,
	               "Output reconstructed image file");

	CLI11_PARSE(app, argc, argv);

	// Your logic here
	std::cout << "Input PETSIRD file: " << input_fname << std::endl;

	/*
	 * Assumtions made by this program:
	 * - There is no DOI (This needs to be fixed!)
	 * - The axial dimension of the scanner is the Z dimension
	 * - The axial dimension in the definition of the crystal box shape is Z
	 *    (The crystal is aligned with the Z axis)
	 * - Crystals are always perpendicular to the Z axis
	 * - All crystals have the same dimensions
	 * */

	// TODO: Make this script compatible with HDF5 files

	// Read PETSIRD FILE
	auto reader = petsird::binary::PETSIRDReader(input_fname);

	// Read the header and get the scanner
	petsird::Header header;
	reader.ReadHeader(header);
	const petsird::ScannerInformation& scannerInfo = header.scanner;
	const petsird::ScannerGeometry& scannerGeom = scannerInfo.scanner_geometry;

	// Get detecting elements
	petsird::TypeOfModule numTypeOfModules =
	    scannerGeom.replicated_modules.size();
	size_t totalNumDets = 0;
	for (petsird::TypeOfModule typeOfModule_i = 0;
	     typeOfModule_i < numTypeOfModules; typeOfModule_i++)
	{
		totalNumDets +=
		    petsird_helpers::get_num_det_els(scannerInfo, typeOfModule_i);
	}

	// PETSIRD "dimensions": type, module, detector
	// Should be reshuffled to: DOI, ring, transaxial (YRT-PET dimensions)

	std::vector<Vector3D> crystalPositions;
	std::vector<Vector3D> crystalOrientations;
	crystalPositions.resize(totalNumDets);
	crystalOrientations.resize(totalNumDets);

	// Crystal properties
	float crystalSize_z, crystalSize_trans, crystalDepth;

	// flat index for the detectors in the YRT-PET LUT
	det_id_t detectorIdx = 0;

	// Properties to measure
	bool isFirstCrystal = true;
	float minZ{};
	float maxZ{};
	float maxDistanceCenter{};

	for (petsird::TypeOfModule typeOfModule_i = 0;
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
		auto centroid = yrt::pet::petsird::getCentroid(currentDetectorShape);

		Vector3D crystalOrientation_yrt;
		std::tie(crystalSize_z, crystalSize_trans, crystalDepth,
		         crystalOrientation_yrt) =
		    yrt::pet::petsird::getCrystalInfo(currentDetectorShape);
		petsird::Coordinate crystalOrientation{{crystalOrientation_yrt.x,
		                                        crystalOrientation_yrt.y,
		                                        crystalOrientation_yrt.z}};

		for (size_t module_i = 0; module_i < numModules; module_i++)
		{
			auto moduleTransform = replicatedModuleType.transforms[module_i];

			// Get the detector orientations (We assume here that the detectors
			// in a module are parallel)

			for (size_t detector_i = 0; detector_i < numDetectors; detector_i++)
			{
				auto detectorTransform = detectors.transforms[detector_i];

				// Transform crystal centroid
				auto totalTransform =
				    petsird_helpers::geometry::mult_transforms(
				        {moduleTransform, detectorTransform});
				auto transformedDetectorCentroid =
				    yrt::pet::petsird::transforms_coord(totalTransform,
				                                        centroid);

				float crystalX = transformedDetectorCentroid.c[0];
				float crystalY = transformedDetectorCentroid.c[1];
				float crystalZ = transformedDetectorCentroid.c[2];

				// Assign position to the detCoord
				crystalPositions[detectorIdx].x = crystalX;
				crystalPositions[detectorIdx].y = crystalY;
				crystalPositions[detectorIdx].z = crystalZ;

				// Remove rotation from the transform
				auto noTranslationMatrix = totalTransform.matrix;
				noTranslationMatrix[3] = 0;
				noTranslationMatrix[7] = 0;
				noTranslationMatrix[11] = 0;
				petsird::RigidTransformation crystalRotation{
				    noTranslationMatrix};
				// Assign orientation to the detCoord
				auto rotatedCrystalOrientation =
				    yrt::pet::petsird::transforms_coord(crystalRotation,
				                                        crystalOrientation);
				crystalOrientations[detectorIdx].x =
				    rotatedCrystalOrientation.c[0];
				crystalOrientations[detectorIdx].y =
				    rotatedCrystalOrientation.c[1];
				crystalOrientations[detectorIdx].z =
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
					distanceCenter =
					    std::max(maxDistanceCenter, distanceCenter);
				}

				// TODO: Compute:
				//  - Get the rotation part of
				//  - Crystal size (z, trans, depth)
				//  - Dets per ring
				//  - Number of rings
				//  - Number of DOI layers
				detectorIdx++;
			}
		}
	}

	float axialFOV = maxZ - minZ;


	// Prepare detCoord
	auto [detCoord, detOriginalIndices, detsPerRing, numRings] =
	    yrt::pet::petsird::toDetCoord(crystalPositions, crystalOrientations);

	// TODO: reshuffle detectors

	if (!outScannerLUT_fname.empty())
	{
		std::cout << "Output scanner LUT file: " << outScannerLUT_fname
		          << std::endl;
		detCoord->writeToFile(outScannerLUT_fname);
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
	                numRings - 1,
	                1,
	                1};
	scanner.setDetectorSetup(detCoord);

	return 0;
}
