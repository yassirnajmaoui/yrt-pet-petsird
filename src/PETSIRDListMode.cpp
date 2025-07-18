#include "PETSIRDListMode.hpp"

#include <petsird_helpers.h>


namespace yrt::pet::petsird
{
	PETSIRDListMode::PETSIRDListMode(
	    const Scanner& pr_scanner,
	    const ::petsird::ScannerInformation& pr_scannerInfo,
	    const DetectorCorrespondenceMap& pr_correspondence,
	    const TimeBlockCollection& pr_timeBlocks)
	    : ListMode(pr_scanner),
	      mr_scannerInfo(pr_scannerInfo),
	      mr_correspondence(pr_correspondence)
	{
		readTimeBlocks(pr_timeBlocks);
	}

	void PETSIRDListMode::readTimeBlocks(const TimeBlockCollection& timeBlocks)
	{
		// TODO: Increase capacity of the std::vectors to increase performance

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
						}
					}
				}
			}
		}
	}
}  // namespace yrt::pet::petsird