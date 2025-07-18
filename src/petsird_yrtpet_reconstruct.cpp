#include "datastruct/projection/ListMode.hpp"
#include "datastruct/scanner/Scanner.hpp"
#include "utils/ReconstructionUtils.hpp"
#include "utils/Utilities.hpp"

#include "PETSIRDListMode.hpp"
#include "utils.hpp"

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
	CLI::App app{"PETSIRD reconstruction executable using YRT-PET"};

	// Variables to hold parsed values
	bool useGPU;
	std::string input_fname;
	int numSubsets = 0;
	int numIterations = 0;
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

	if (Util::compiledWithCuda())
	{
		app.add_flag("--gpu", useGPU, "Use GPU acceleration");
	}

	app.add_option("--num-subsets", numSubsets, "Number of subsets")
	    ->default_val(1);

	app.add_option("--num-iterations", numIterations, "Number of iterations")
	    ->required();

	app.add_option("-p, --params", imageParams_fname, "Image parameters file")
	    ->required()
	    ->check(CLI::ExistingFile);

	app.add_option("--psf", psfKernel_fname, "PSF kernel file")
	    ->check(CLI::ExistingFile);

	app.add_option("--att", attImage_fname, "Attenuation image file")
	    ->check(CLI::ExistingFile);

	app.add_option("--out-scanner-lut", outScannerLUT_fname,
	               "Output scanner LUT file");
	// app.add_option("--out-scanner-json", outScannerJSON_fname,
	//                "Output scanner JSON file");

	app.add_option("--sens", sensImage_fname,
	               "Pre-existing sensitivity image filename");
	app.add_option("--out-sens", outSensImage_fname,
	               "Output sensitivity image file");
	app.add_option("--out-recon", outImage_fname,
	               "Output reconstructed image file")
	    ->required();

	CLI11_PARSE(app, argc, argv);

	// Your logic here
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
	//  - Read TOF resolution from the PETSIRD file and add it in "addTOF"
	//  - Convert TOF indices into values in picoseconds and use them for recon
	//  - Add possibility to provide an attenuation image
	//  - Add normalisation


	// Read PETSIRD FILE
	auto reader = petsird::binary::PETSIRDReader(input_fname);

	// Read the header and get the scanner
	petsird::Header header;
	reader.ReadHeader(header);
	const petsird::ScannerInformation& scannerInfo = header.scanner;

	// Prepare detCoord
	auto [scanner, correspondenceMap] =
	    yrt::pet::petsird::toScanner(scannerInfo);

	if (!outScannerLUT_fname.empty())
	{
		std::cout << "Output scanner LUT file: " << outScannerLUT_fname
		          << std::endl;
		auto detSetup = scanner.getDetectorSetup();
		detSetup->writeToFile(outScannerLUT_fname);
	}

	// TODO: Save the scanner's JSON file

	// ListMode l = yrt::pet::petsird::PETSIRDListMode();
	//  Read the header and get the scanner
	yrt::pet::petsird::TimeBlockCollection timeBlocks;
	timeBlocks.reserve(50ull << 10);

	const bool readingTimeBlock = reader.ReadTimeBlocks(timeBlocks);
	if (!readingTimeBlock)
	{
		throw std::runtime_error("Error while reading time blocks");
	}

	yrt::pet::petsird::PETSIRDListMode lm(scanner, scannerInfo,
	                                      correspondenceMap, timeBlocks);

	// Initialize reconstruction
	auto osem = Util::createOSEM(scanner, useGPU);
	osem->setListModeEnabled(true);

	// Read image parameters
	ImageParams params{imageParams_fname};
	osem->setImageParams(params);

	if (!psfKernel_fname.empty())
	{
		osem->addImagePSF(psfKernel_fname);
	}

	std::unique_ptr<Image> attImage;
	if (!attImage_fname.empty())
	{
		attImage = std::make_unique<ImageOwned>(attImage_fname);
		osem->setAttenuationImage(attImage.get());
	}

	std::vector<std::unique_ptr<Image>> sensImages;
	if (!sensImage_fname.empty())
	{
		sensImages.push_back(std::make_unique<ImageOwned>(sensImage_fname));
	}
	else
	{
		osem->generateSensitivityImages(sensImages, outSensImage_fname);
	}

	osem->setSensitivityImages(sensImages);

	osem->setDataInput(&lm);

	osem->num_MLEM_iterations = numIterations;
	osem->num_OSEM_subsets = numSubsets;

	osem->reconstruct(outImage_fname);

	std::cout << "Done." << std::endl;

	return 0;
}
