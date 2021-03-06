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
    enum Scheme {CMB=-1, EUL=2, IMP_3=3, IMP_4=4, EXP=5, TG=6, QE=7, EXT=9};
    static const Scheme ALL_SCHEME[] = {CMB, EUL, IMP_3, IMP_4, EXP, TG, QE, EXT};
    static const Scheme ALL_CIR_SCHEME[] = {EUL, IMP_3, IMP_4, EXP, QE, EXT};
    // enum CIRScheme {EUL=2, IMP_3=3, IMP_4=4, EXP=5, TG=6, QE=7, EXT=9};
    // enum HESScheme {EUL=2, CMB=3};
}