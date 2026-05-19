/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ThreeWayPerceiverSenderAgreeEffect.h
 *
 * Description: This file contains the definition of the
 * ThreeWayPerceiverSenderAgreeEffect class.
 *
 * Three-way SAOM: percSenderAgree effect (Eq 3.12 in the paper).
 *
 * Statistic for slice Y[i]:
 *   s(x) = sum_{j,k} x_{ijk} * x_{jjk} * (I(v_i == v_j) - mu_sim)
 *
 * where
 *   i   = perceiver (slice owner, 0-based index lperceiverIndex)
 *   j   = ego (sender whose micro-step is evaluated)
 *   k   = alter (receiver of the potential tie)
 *   x_{ijk}  = tie in perceived slice Y[i]  (the current network)
 *   x_{jjk}  = j's self-reported tie to k,  read from Y[self]
 *
 * The perceiver index is parsed at initialization from the variable name
 * "Y[n]" (1-based n on the R side → 0-based index here).
 * The self-reported network is located by replacing "[n]" with "[self]".
 *****************************************************************************/

#ifndef THREEWAYPERCEIVERENDERAGREEEFFECT_H_
#define THREEWAYPERCEIVERENDERAGREEEFFECT_H_

#include "CovariateDependentNetworkEffect.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;

// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Perceiver-sender-agreement network effect for three-way SAOM perceived slices.
 */
class ThreeWayPerceiverSenderAgreeEffect : public CovariateDependentNetworkEffect
{
public:
	explicit ThreeWayPerceiverSenderAgreeEffect(
		const EffectInfo * pEffectInfo);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);

private:
	/// 0-based actor index of the slice owner (perceiver i)
	int lperceiverIndex;

	/// Pointer to the current Y[self] network (not owned; valid within period)
	const Network * lpSelfNetwork;
};

}

#endif /*THREEWAYPERCEIVERENDERAGREEEFFECT_H_*/
