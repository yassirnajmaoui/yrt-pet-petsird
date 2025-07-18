
#pragma once

#include "datastruct/projection/ListMode.hpp"
#include "utils.hpp"

namespace yrt::pet::petsird
{
	class PETSIRDListMode : public ListMode
	{
	public:
		PETSIRDListMode(const Scanner& pr_scanner,
		                const ::petsird::ScannerInformation& pr_scannerInfo,
		                const DetectorCorrespondenceMap& pr_correspondence,
		                const TimeBlockCollection& pr_timeBlocks);

		// Appends the events in the given time blocks into the list of events
		void readTimeBlocks(const TimeBlockCollection& timeBlocks);

	private:
		const DetectorCorrespondenceMap mr_correspondence;
		const ::petsird::ScannerInformation& mr_scannerInfo;

		std::vector<timestamp_t> timestamps;  // in ms
		std::vector<det_id_t> d1s;            // index in the YRT-PET LUT
		std::vector<det_id_t> d2s;            // index in the YRT-PET LUT
		std::vector<float> tofs;              // in ps
		                                      // TODO: Motion
	};
}  // namespace yrt::pet::petsird
