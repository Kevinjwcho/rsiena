/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ThreeWayPerceivedReciprocityEffect.h   *** DEPRECATED ORPHAN ***
 *
 * This file is no longer compiled into the package (removed from
 * AllEffects.h, EffectFactory.cpp, sources.list, src/Makevars.win).
 *
 * Rationale:
 *   The original intra-slice statistic
 *       s_i = sum_k Y[i]_{i,k} * Y[i]_{k,i}
 *   is identically zero in this three-way fork because
 *   initializeFRAN.r:2289-2296 marks perceiver i's own row of slice Y[i]
 *   as structural-zero (value 10) — perceiver's self-reported ties live
 *   exclusively in Y[self], not in the perceived slice.  So the
 *   statistic above always evaluates to 0 and recovery is impossible.
 *
 *   The corrected cross-network statistic
 *       s = sum_i sum_j Y[self]_{i,j} * Y[i]_{j,i}
 *   is now implemented via GenericNetworkEffect +
 *   PerceiverRestrictedInTieFunction; see
 *   src/model/effects/generic/PerceiverRestrictedInTieFunction.{h,cpp}
 *   and the EffectFactory.cpp "percRecip" branch.
 *
 * Safe to delete this file (and the matching .cpp) once you're sure no
 * downstream tooling still references the class name.
 *****************************************************************************/

#ifndef THREEWAYPERCEIVEDRECIPROCITYEFFECT_H_
#define THREEWAYPERCEIVEDRECIPROCITYEFFECT_H_

// Intentionally empty — see notice above.

#endif /*THREEWAYPERCEIVEDRECIPROCITYEFFECT_H_*/
