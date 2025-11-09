#include "PETSIRDListMode.hpp"

#include <petsird_helpers.h>

namespace yrt::petsird
{
	PETSIRDListMode::PETSIRDListMode(
	    const Scanner& pr_scanner,
	    const ::petsird::ScannerInformation& pr_scannerInfo,
	    const DetectorCorrespondenceMap& pr_correspondence,
	    const TimeBlockCollection& pr_timeBlocks)
	    : ListMode(pr_scanner),
	      mr_correspondence(pr_correspondence),
	      mr_scannerInfo(pr_scannerInfo)
	{
		readTimeBlocks(pr_timeBlocks);
	}

	void PETSIRDListMode::readTimeBlocks(const TimeBlockCollection& timeBlocks)
	{
		// TODO: Increase capacity of the std::vectors to increase performance

		// const size_t numEventsEstimate = ???
		// m_timestamps.reserve(numEventsEstimate);
		// m_d0s.reserve(numEventsEstimate);
		// m_d1s.reserve(numEventsEstimate);

		const size_t numTimeBlocks = timeBlocks.size();
		timestamp_t currentTime{};

		for (size_t timeBlock_i = 0; timeBlock_i < numTimeBlocks; timeBlock_i++)
		{
			const auto& timeBlock = timeBlocks[timeBlock_i];
			if (std::holds_alternative<::petsird::EventTimeBlock>(timeBlock))
			{
				const auto& eventTimeBlock =
				    std::get<::petsird::EventTimeBlock>(timeBlock);
				currentTime = eventTimeBlock.time_interval.start;

				// Here we only accumulate prompt events
				const auto& promptEvents = eventTimeBlock.prompt_events;

				const size_t numTypesOfModules = promptEvents.size();

				for (::petsird::TypeOfModule mtype0 = 0;
				     mtype0 < numTypesOfModules; mtype0++)
				{
					const auto& promptEvents_mtype0 = promptEvents[mtype0];

					if (promptEvents_mtype0.size() != numTypesOfModules)
					{
						throw std::runtime_error(
						    "File is not properly formed: The number of module "
						    "types is not consistent in the list-mode events.");
					}

					for (::petsird::TypeOfModule mtype1 = 0;
					     mtype1 < numTypesOfModules; mtype1++)
					{
						const auto& promptEvents_mtype01 =
						    promptEvents_mtype0[mtype1];
						for (const auto& promptEvent : promptEvents_mtype01)
						{
							auto [d0_expanded, d1_expanded] =
							    petsird_helpers::expand_detection_bin_pair(
							        mr_scannerInfo, {mtype0, mtype1},
							        promptEvent.detection_bins);
							det_id_t d0flatIdx = mr_correspondence.getFlatIndex(
							    mtype0, d0_expanded.module_index,
							    d0_expanded.element_index);
							det_id_t d1flatIdx = mr_correspondence.getFlatIndex(
							    mtype1, d1_expanded.module_index,
							    d1_expanded.element_index);

							promptEvent.tof_idx;  // TODO: Store TOF value in ps

							// Add to the cumulative vectors
							m_timestamps.emplace_back(currentTime);
							m_d0s.emplace_back(d0flatIdx);
							m_d1s.emplace_back(d1flatIdx);
						}
					}
				}
			}
		}
	}

	det_id_t PETSIRDListMode::getDetector1(bin_t id) const
	{
		return m_d0s[id];
	}

	det_id_t PETSIRDListMode::getDetector2(bin_t id) const
	{
		return m_d1s[id];
	}

	det_pair_t PETSIRDListMode::getDetectorPair(bin_t id) const
	{
		return {m_d0s[id], m_d1s[id]};
	}

	size_t PETSIRDListMode::count() const
	{
		return m_d0s.size();
	}

	timestamp_t PETSIRDListMode::getTimestamp(bin_t id) const
	{
		return m_timestamps[id];
	}

	bool PETSIRDListMode::hasTOF() const
	{
		// TODO: Figure out how to read TOF value in ps from idx
		return false;
	}

	float PETSIRDListMode::getTOFValue(bin_t id) const
	{
		(void) id;
		return 0;  // Will be ignored because of hasTOF
	}


}  // namespace yrt::pet::petsird