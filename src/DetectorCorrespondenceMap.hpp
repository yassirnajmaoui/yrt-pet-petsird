#pragma once

#include "utils/Types.hpp"

#include <stdexcept>
#include <unordered_map>

namespace yrt::pet::petsird
{
	class DetectorCorrespondenceMap
	{
	public:
		struct DetectorKey
		{
			uint32_t type;
			uint32_t module;
			uint32_t det;

			// Needed for unordered_map
			bool operator==(const DetectorKey& other) const
			{
				return type == other.type && module == other.module &&
				       det == other.det;
			}
		};

		// Set the correspondence
		void addMapping(uint32_t type, uint32_t module, uint32_t det,
		                det_id_t value);

		// Get the correspondence
		det_id_t getFlatIndex(uint32_t type, uint32_t module,
		                      uint32_t det) const;

		// Get (type, module, det) from value
		std::tuple<uint32_t, uint32_t, uint32_t>
		    getDetectorFromFlatIndex(det_id_t value) const;

		// Check if a detector is present
		bool contains(uint32_t type, uint32_t module, uint32_t det) const;

	private:
		// Custom hash function
		struct DetectorKeyHash
		{
			std::size_t operator()(const DetectorKey& k) const
			{
				return std::hash<uint32_t>()(k.type) ^
				       (std::hash<uint32_t>()(k.module) << 1) ^
				       (std::hash<uint32_t>()(k.det) << 2);
			}
		};

		std::unordered_map<DetectorKey, det_id_t, DetectorKeyHash> map;
		std::unordered_map<det_id_t, DetectorKey> reverseMap;
	};
}  // namespace yrt::pet::petsird
