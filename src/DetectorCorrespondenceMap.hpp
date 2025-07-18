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
		                det_id_t value)
		{
			const DetectorKey key{type, module, det};
			map[key] = value;
			reverseMap[value] = key;
		}

		// Get the correspondence
		det_id_t getFlatIndex(uint32_t type, uint32_t module,
		                      uint32_t det) const
		{
			const DetectorKey key{type, module, det};
			const auto it = map.find(key);
			if (it == map.end())
			{
				throw std::out_of_range("Detector not found in map.");
			}
			return it->second;
		}

		// Get (type, module, det) from value
		std::tuple<uint32_t, uint32_t, uint32_t>
		    getDetectorFromFlatIndex(det_id_t value) const
		{
			const auto it = reverseMap.find(value);
			if (it == reverseMap.end())
				throw std::out_of_range("Value not found.");
			const auto& key = it->second;
			return std::make_tuple(key.type, key.module, key.det);
		}

		// Check if a detector is present
		bool contains(uint32_t type, uint32_t module, uint32_t det) const
		{
			return map.find({type, module, det}) != map.end();
		}

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
