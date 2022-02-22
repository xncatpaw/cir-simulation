/**
 * @file simu_type.hpp
 * @author bx.zh 
 * @brief defines the type of simulation.
 * @date 2022-02-10
 * 
 * @copyright Copyleft (c) all wrongs reserved.
 */
#pragma once

namespace cir
{
    enum SimuType {PRECISE, SECOND_ORDER=2, THIRD_ORDER=3};
    enum SimuScheme {IMP_3=3, IMP_4=4, EXP=5};
}