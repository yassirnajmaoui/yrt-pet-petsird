#include "DetectorCorrespondenceMap.hpp"

void yrt::petsird::DetectorCorrespondenceMap::addMapping(uint32_t type,
                                                              uint32_t module,
                                                              uint32_t det,
                                                              yrt::det_id_t value)
{
	const DetectorKey key{type, module, det};
	map[key] = value;
	reverseMap[value] = key;
}

yrt::det_id_t yrt::petsird::DetectorCorrespondenceMap::getFlatIndex(
    uint32_t type, uint32_t module, uint32_t det) const
{
	const DetectorKey key{type, module, det};
	const auto it = map.find(key);
	if (it == map.end())
	{
		throw std::out_of_range("Detector not found in map.");
	}
	return it->second;
}

std::tuple<uint32_t, uint32_t, uint32_t>
    yrt::petsird::DetectorCorrespondenceMap::getDetectorFromFlatIndex(
        yrt::det_id_t value) const
{
	const auto it = reverseMap.find(value);
	if (it == reverseMap.end())
		throw std::out_of_range("Value not found.");
	const auto& key = it->second;
	return std::make_tuple(key.type, key.module, key.det);
}

bool yrt::petsird::DetectorCorrespondenceMap::contains(uint32_t type,
                                                            uint32_t module,
                                                            uint32_t det) const
{
	return map.find({type, module, det}) != map.end();
}
