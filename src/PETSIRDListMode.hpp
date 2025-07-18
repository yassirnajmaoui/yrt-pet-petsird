
#pragma once

#include "DetectorCorrespondenceMap.hpp"
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

		det_id_t getDetector1(bin_t id) const override;
		det_id_t getDetector2(bin_t id) const override;
		det_pair_t getDetectorPair(bin_t id) const override;
		size_t count() const override;
		timestamp_t getTimestamp(bin_t id) const override;

		bool hasTOF() const override;
		float getTOFValue(bin_t id) const override;

	private:
		const DetectorCorrespondenceMap mr_correspondence;
		const ::petsird::ScannerInformation& mr_scannerInfo;

		std::vector<timestamp_t> m_timestamps;  // in ms
		std::vector<det_id_t> m_d0s;            // index in the YRT-PET LUT
		std::vector<det_id_t> m_d1s;            // index in the YRT-PET LUT
		std::vector<float> m_tofs;              // in ps
		                                        // TODO: Motion
	};
}  // namespace yrt::pet::petsird
