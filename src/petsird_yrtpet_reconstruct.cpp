#include "datastruct/scanner/Scanner.hpp"
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
	std::cout << "Input PETSIRD file: " << input_fname << "\n"
	          << "Number of subsets: " << numSubsets << "\n"
	          << "Number of iterations: " << numIterations << "\n"
	          << "PSF file: " << psfKernel_fname << "\n"
	          << "Output sensitivity image: " << outSensImage_fname << "\n"
	          << "Output reconstructed image: " << outImage_fname << std::endl;

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

	// Prepare detCoord
	DetCoordOwned detCoord;
	detCoord.allocate(totalNumDets);

	// PETSIRD "dimensions": type, module, detector
	// Should be reshuffled to: DOI, ring, transaxial (YRT-PET dimensions)

	// flat index for the detectors in the YRT-PET LUT
	det_id_t detectorIdx = 0;

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
		// TODO: hide this in a helper function in the beginning of this file
		const auto& currentDetectorShape = detectors.object.shape;
		size_t numCorners = currentDetectorShape.corners.size();
		petsird::Coordinate centroid;
		for (size_t i = 0; i < numCorners; ++i)
		{
			for (int dim = 0; dim < 3; dim++)
			{
				centroid.c[dim] += currentDetectorShape.corners[i].c[dim];
			}
		}
		for (int dim = 0; dim < 3; dim++)
		{
			// Divide the sum by number of corners to get the average
			centroid.c[dim] = centroid.c[dim] / static_cast<float>(numCorners);
		}

		for (size_t module_i = 0; module_i < numModules; module_i++)
		{
			auto moduleTransform = replicatedModuleType.transforms[module_i];
			for (size_t detector_i = 0; detector_i < numDetectors; detector_i++)
			{
				auto detectorTransform = detectors.transforms[detector_i];

				// Transform crystal centroid
				auto transformedDetectorCentroid =
				    petsird_helpers::geometry::mult_transforms_coord(
				        {moduleTransform, detectorTransform}, centroid);

				// Assign to the detCoord
				detCoord.setXpos(detectorIdx, transformedDetectorCentroid.c[0]);
				detCoord.setYpos(detectorIdx, transformedDetectorCentroid.c[1]);
				detCoord.setZpos(detectorIdx, transformedDetectorCentroid.c[2]);
				// TODO: Add orientation too
				detectorIdx++;
			}
		}
	}

	if (!outScannerLUT_fname.empty())
	{
		detCoord.writeToFile(outScannerLUT_fname);
	}

	std::cout << Util::compiledWithCuda() << std::endl;

	return 0;
}
