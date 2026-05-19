/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ThreeWayPerceiverEffect.h
 *
 * Description: This file contains the definition of the
 * ThreeWayPerceiverEffect class.
 *
 * Three-way SAOM: perceiver-attribute effects operating on a perceived slice
 * Y[i].  The perceiver (slice owner) index i is parsed at initialization from
 * the variable name string "Y[n]" (1-based n → 0-based C++ index).
 *
 * Supported modes:
 *   TW_PERC_EGO       percX         contribution = (v_i - mean_v)
 *   TW_PERC_ALT_SAME  percAltSame   contribution = sim(v_i, v_k)
 *   TW_PERC_SEND_SAME percSenderSame contribution = sim(v_i, v_j)
 *
 * where i = perceiver, j = ego (sender in the slice), k = alter (receiver).
 * "sim" is the mean-centered similarity: I(v_a == v_b) - mu_sim.
 *****************************************************************************/

#ifndef THREEWAYPERCEIVEREFFECT_H_
#define THREEWAYPERCEIVEREFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Mode enum
// ----------------------------------------------------------------------------

enum ThreeWayPerceiverMode
{
	TW_PERC_EGO,        ///< percX:          v_i (mean-centered)
	TW_PERC_ALT_SAME,   ///< percAltSame:    sim(v_i, v_alter)
	TW_PERC_SEND_SAME   ///< percSenderSame: sim(v_i, v_ego)
};


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Perceiver-attribute network effects for three-way SAOM perceived slices.
 */
class ThreeWayPerceiverEffect : public CovariateDependentNetworkEffect
{
public:
	ThreeWayPerceiverEffect(const EffectInfo * pEffectInfo,
		ThreeWayPerceiverMode mode);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	ThreeWayPerceiverMode lmode;

	/// 0-based actor index of the slice owner (perceiver i)
	int lperceiverIndex;
};

}

#endif /*THREEWAYPERCEIVEREFFECT_H_*/
