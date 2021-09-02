#ifndef _CACHEDDATA_HPP
#define _CACHEDDATA_HPP

#include "meas_vec.hpp"

#include <externer.hpp>
#include <params.hpp>

constexpr char cachedAveCorrPath[]="cachedAveCorr.dat";

EXTERN map<AveId,djvec_t> aveCorrCache;
EXTERN bool hasTostoreCachedAveCorr INIT_EXTERN_TO(=false);

bool loadCachedAveCorr();

void storeCachedAveCorr();

#endif
