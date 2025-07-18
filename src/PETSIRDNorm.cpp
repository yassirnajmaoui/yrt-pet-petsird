#include "PETSIRDNorm.hpp"

#include "petsird_helpers.h"

namespace yrt::pet::petsird
{
	PETSIRDNorm::PETSIRDNorm(
	    const Scanner& pr_scanner,
	    const ::petsird::ScannerInformation& pr_scannerInfo,
	    const DetectorCorrespondenceMap& pr_correspondence)
	    : Histogram3D(pr_scanner),
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


}  // namespace yrt::pet::petsird
