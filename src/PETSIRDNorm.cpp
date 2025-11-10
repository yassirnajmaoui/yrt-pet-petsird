#include "PETSIRDNorm.hpp"

#include "petsird_helpers.h"
#include "yrt-pet/datastruct/scanner/Scanner.hpp"

namespace yrt::petsird
{
	PETSIRDNorm::PETSIRDNorm(
	    const Scanner& pr_scanner,
	    const ::petsird::ScannerInformation& pr_scannerInfo,
	    const DetectorCorrespondenceMap& pr_correspondence)
	    : yrt::Histogram3D(pr_scanner),
	      mr_correspondence(pr_correspondence),
	      mr_scannerInfo(pr_scannerInfo)
	{
	}

	float PETSIRDNorm::getProjectionValue(bin_t binId) const
	{
		const det_pair_t detPair = getDetectorPair(binId);
		// Higher number -> more sensitive
		// TODO: implement this. But need to somehow have "histogram bins" be
		//  the "detection_bin" that petsird defines

		auto d1_tuple = mr_correspondence.getDetectorFromFlatIndex(detPair.d1);
		auto d2_tuple = mr_correspondence.getDetectorFromFlatIndex(detPair.d2);

		::petsird::ExpandedDetectionBin d1_expandedBin{};
		d1_expandedBin.module_index = std::get<1>(d1_tuple);
		d1_expandedBin.element_index = std::get<2>(d1_tuple);
		// TODO: Support energy level
		auto d1_bin = petsird_helpers::make_detection_bin(
		    mr_scannerInfo, std::get<0>(d1_tuple), d1_expandedBin);

		::petsird::ExpandedDetectionBin d2_expandedBin{};
		d2_expandedBin.module_index = std::get<1>(d2_tuple);
		d2_expandedBin.element_index = std::get<2>(d2_tuple);
		// TODO: Support energy level
		auto d2_bin = petsird_helpers::make_detection_bin(
		    mr_scannerInfo, std::get<0>(d2_tuple), d2_expandedBin);

		float sensitivity = petsird_helpers::get_detection_efficiency(
		    mr_scannerInfo, {std::get<0>(d1_tuple), std::get<0>(d2_tuple)},
		    d1_bin, d2_bin);

		if (sensitivity != 1.0 && sensitivity != 0)
		{
			std::cout << "Sensitivity: " << sensitivity << std::endl;
		}

		return sensitivity;
	}

	void PETSIRDNorm::setProjectionValue(bin_t binId, float val)
	{
		(void)binId;
		(void)val;
		throw std::logic_error("Unimplemented");
	}

	void PETSIRDNorm::incrementProjection(bin_t binId, float val)
	{
		(void)binId;
		(void)val;
		throw std::logic_error("Unimplemented");
	}

	void PETSIRDNorm::clearProjections(float p_value)
	{
		(void)p_value;
		throw std::logic_error("Unimplemented");
	}


}  // namespace yrt::petsird
