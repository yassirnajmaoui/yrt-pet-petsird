
#include "datastruct/scanner/Scanner.hpp"
#include "petsird_helpers.h"

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
	std::string outSensImage_fname;
	std::string outImage_fname;

	// Add options
	app.add_option("-i,--input", input_fname, "Input PETSIRD file")
	    ->required()
	    ->check(CLI::ExistingFile);

	app.add_option("-s,--subsets", numSubsets, "Number of subsets")->required();

	app.add_option("-n,--iterations", numIterations, "Number of iterations")
	    ->required();

	app.add_option("--psf", psfKernel_fname, "PSF kernel file")
	    ->required()
	    ->check(CLI::ExistingFile);

	app.add_option("--out-sens", outSensImage_fname,
	               "Output sensitivity image file")
	    ->required();

	app.add_option("--out-recon", outImage_fname,
	               "Output reconstructed image file")
	    ->required();

	CLI11_PARSE(app, argc, argv);

	// Your logic here
	std::cout << "Input PETSIRD file: " << input_fname << "\n"
	          << "Number of subsets: " << numSubsets << "\n"
	          << "Number of iterations: " << numIterations << "\n"
	          << "PSF file: " << psfKernel_fname << "\n"
	          << "Output Sensitivity image: " << outSensImage_fname << "\n"
	          << "Output Reconstructed image: " << outImage_fname << "\n";


	return 0;
}
