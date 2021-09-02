#ifndef _DATA_HPP
#define _DATA_HPP

#include <array>
#include <string>

#include <meas_vec.hpp>

#include <externer.hpp>
#include <params.hpp>

using namespace std;

djvec_t getAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const size_t iMes);
djvec_t getTMAve(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb);

djvec_t getAveForRego(const size_t iSourceMin,const size_t iSourceMax,const size_t iGammaComb,const RegoType& rego);

#endif
