#ifndef _AMU_HPP
#define _AMU_HPP

#include <meas_vec.hpp>

#include <VKVKRepresentation.hpp>

void computeAmu(const RegoType& rego,
		const jack_t<VKVKRepFiniteVol>& rep,
		const jack_t<VKVKRepInfiniteVol>& infVolRep);

#endif
