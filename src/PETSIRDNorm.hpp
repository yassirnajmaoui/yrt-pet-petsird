#pragma once

#include "DetectorCorrespondenceMap.hpp"
#include "datastruct/projection/Histogram3D.hpp"
#include "protocols.h"

namespace yrt::pet::petsird
{
	class PETSIRDNorm final : public Histogram3D
	{
	public:
		PETSIRDNorm(const Scanner& pr_scanner,
		            const ::petsird::ScannerInformation& pr_scannerInfo,
		            const DetectorCorrespondenceMap& pr_correspondence);

		float getProjectionValue(bin_t binId) const override;

		// N/A
		void setProjectionValue(bin_t binId, float val) override;
		void incrementProjection(bin_t binId, float val) override;
		void clearProjections(float p_value) override;


	private:
		const DetectorCorrespondenceMap mr_correspondence;
		const ::petsird::ScannerInformation& mr_scannerInfo;
	};
}  // namespace yrt::pet::petsird
