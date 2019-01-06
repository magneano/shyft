/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once
#include "core_serialization.h"
namespace shyft::core {

/** @brief method stack parameters that works across methods
 *
 * Most parameters are connected to each method used in the stack
 * This parameter structure contains the one that controls
 * the internal flow between the method-layers and/or composition.
 * 
 */
struct mstack_parameter {
    mstack_parameter()=default;
    explicit mstack_parameter( double r_direct_response_fraction):reservoir_direct_response_fraction{r_direct_response_fraction}{}
    double reservoir_direct_response_fraction{1.0};///< range 0..1,unit-less, semantic:precipitation, in any form that hits reservoirs goes to direct response(no routing)  
    x_serialize_decl();
};
    
}
x_serialize_export_key(shyft::core::mstack_parameter);
