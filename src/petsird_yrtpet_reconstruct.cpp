#include "yrt-pet/datastruct/projection/ListMode.hpp"
#include "yrt-pet/datastruct/scanner/Scanner.hpp"
#include "yrt-pet/utils/ReconstructionUtils.hpp"
#include "yrt-pet/utils/Utilities.hpp"

#include "PETSIRDListMode.hpp"
#include "PETSIRDNorm.hpp"
#include "utils.hpp"

#include "hdf5.h"
#include "petsird/binary/protocols.h"
#include "petsird/hdf5/protocols.h"
#include "petsird/protocols.h"
#include "petsird/types.h"
#include "petsird_helpers/create.h"
#include "petsird_helpers/geometry.h"

#include "CLI11.hpp"
#include <string>


int main(int argc, char** argv)
{
	CLI::App app{"PETSIRD reconstruction executable using YRT-PET"};

	// Variables to hold parsed values
	bool useTOF;
	bool useGPU;
	bool useNorm;
	std::string input_fname;
	int numSubsets = 0;
	int numIterations = 0;
	int numThreads = -1;
	std::string imageParams_fname;
	std::string psfKernel_fname;
	std::string attImage_fname;
	std::string outScannerLUT_fname;
	std::string outScannerJSON_fname;
	std::string outSensImage_fname;
	std::string sensImage_fname;
	std::string outImage_fname;

	// Add options
	app.add_option("-i,--input", input_fname, "Input PETSIRD file")
	    ->required()
	    ->check(CLI::ExistingFile);

	if (yrt::util::compiledWithCuda())
	{
		app.add_flag("--gpu", useGPU, "Use GPU acceleration");
	}

	app.add_option("--num_threads", numThreads, "Number of threads to use");

	app.add_option("--num_subsets", numSubsets, "Number of subsets")
	    ->default_val(1);

	app.add_option("--num_iterations", numIterations, "Number of iterations")
	    ->default_val(10);

	app.add_option("-p, --params", imageParams_fname, "Image parameters file")
	    ->required()
	    ->check(CLI::ExistingFile);

	app.add_option("--psf", psfKernel_fname, "PSF kernel file")
	    ->check(CLI::ExistingFile);

	app.add_option("--att", attImage_fname, "Attenuation image file")
	    ->check(CLI::ExistingFile);

	app.add_flag("--norm", useNorm, "Apply normalisation correction");

	app.add_flag("--tof", useTOF, "Use TOF information");

	app.add_option("--out_scanner_lut", outScannerLUT_fname,
	               "Output scanner LUT file");
	// app.add_option("--out-scanner-json", outScannerJSON_fname,
	//                "Output scanner JSON file");

	app.add_option("--sens", sensImage_fname,
	               "Pre-existing sensitivity image filename");
	app.add_option("--out_sens", outSensImage_fname,
	               "Output sensitivity image file");
	app.add_option("-o, --out", outImage_fname,
	               "Output reconstructed image file")
	    ->required();

	CLI11_PARSE(app, argc, argv);

	std::cout << "Input PETSIRD file: " << input_fname << std::endl;

	/*
	 * Assumptions made by this program:
	 * - There is no DOI (This needs to be fixed!)
	 * - The axial dimension of the scanner is the Z dimension
	 * - The axial dimension in the definition of the crystal box shape is Z
	 *    (The crystal is aligned with the Z axis)
	 * - Crystals are always perpendicular to the Z axis (no helmets yet)
	 * - All crystals have the same dimensions
	 * */

	// TODO:
	//  - Make this script compatible with HDF5 files

	yrt::globals::setNumThreads(numThreads);

	// Read PETSIRD FILE
	auto reader = petsird::binary::PETSIRDReader(input_fname);

	// Read the header and get the scanner
	petsird::Header header;
	reader.ReadHeader(header);
	const petsird::ScannerInformation& scannerInfo = header.scanner;

	// Prepare detCoord
	auto [scanner, correspondenceMap] = yrt::petsird::toScanner(scannerInfo);

	if (!outScannerLUT_fname.empty())
	{
		std::cout << "Output scanner LUT file: " << outScannerLUT_fname
		          << std::endl;
		auto detSetup = scanner.getDetectorSetup();
		detSetup->writeToFile(outScannerLUT_fname);
	}

	// TODO: Save the scanner's JSON file

	// ListMode l = yrt::petsird::PETSIRDListMode();
	//  Read the header and get the scanner
	yrt::petsird::TimeBlockCollection timeBlocks;
	timeBlocks.reserve(1ull << 31);

	const bool readingTimeBlock = reader.ReadTimeBlocks(timeBlocks);
	if (!readingTimeBlock)
	{
		throw std::runtime_error("Error while reading time blocks");
	}

	auto lm = std::make_unique<yrt::petsird::PETSIRDListMode>(
	    scanner, scannerInfo, correspondenceMap, timeBlocks, useTOF);

	// Initialize reconstruction
	auto osem = yrt::util::createOSEM(scanner, useGPU);
	osem->setListModeEnabled(true);

	// Read image parameters
	yrt::ImageParams params{imageParams_fname};
	osem->setImageParams(params);

	if (!psfKernel_fname.empty())
	{
		osem->addImagePSF(psfKernel_fname);
	}

	std::unique_ptr<yrt::petsird::PETSIRDNorm> norm;
	if (useNorm)
	{
		norm = std::make_unique<yrt::petsird::PETSIRDNorm>(scanner, scannerInfo,
		                                                   correspondenceMap);
		osem->setSensitivityHistogram(norm.get());
	}

	std::unique_ptr<yrt::Image> attImage;
	if (!attImage_fname.empty())
	{
		attImage = std::make_unique<yrt::ImageOwned>(attImage_fname);
		osem->setAttenuationImage(attImage.get());
	}

	std::vector<std::unique_ptr<yrt::Image>> sensImages;
	if (!sensImage_fname.empty())
	{
		sensImages.push_back(
		    std::make_unique<yrt::ImageOwned>(sensImage_fname));
	}
	else
	{
		osem->generateSensitivityImages(sensImages, outSensImage_fname);
	}

	osem->setSensitivityImages(sensImages);

	osem->setDataInput(lm.get());

	osem->num_MLEM_iterations = numIterations;
	osem->num_OSEM_subsets = numSubsets;

	if (useTOF)
	{
		float tofResolution_ps =
		    scannerInfo.tof_resolution[0][0] * 2.0f / yrt::SPEED_OF_LIGHT_MM_PS;
		osem->addTOF(tofResolution_ps, 5);
	}

	osem->reconstruct(outImage_fname);

	std::cout << "Done." << std::endl;

	return 0;
}
