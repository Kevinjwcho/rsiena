/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ThreeWayPerceiverSenderAgreeEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ThreeWayPerceiverSenderAgreeEffect class.
 *****************************************************************************/

#include <string>
#include <stdexcept>

#include "ThreeWayPerceiverSenderAgreeEffect.h"
#include "model/EffectInfo.h"
#include "model/State.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo  the effect descriptor
 */
ThreeWayPerceiverSenderAgreeEffect::ThreeWayPerceiverSenderAgreeEffect(
		const EffectInfo * pEffectInfo) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lperceiverIndex = -1;
	this->lpSelfNetwork   = 0;
}


/**
 * Initializes this effect.
 *
 * Calls the parent initializer (which resolves interactionName1 to a
 * covariate), then:
 *   (a) parses the perceiver index from the variable name "Y[n]", and
 *   (b) looks up the self-reported companion network whose name is formed
 *       by replacing "[n]" with "[self]" (e.g., "Y[self]").
 *
 * @param[in] pData    the observed data
 * @param[in] pState   state of dependent variables at the start of the period
 * @param[in] period   the current period
 * @param[in] pCache   cache object for speed-up calculations
 */
void ThreeWayPerceiverSenderAgreeEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentNetworkEffect::initialize(pData, pState, period, pCache);

	// Parse perceiver index from variable name "Y[n]" (1-based → 0-based).
	string varName = this->pEffectInfo()->variableName();
	string::size_type lb = varName.find('[');
	string::size_type rb = varName.find(']');

	if (lb == string::npos || rb == string::npos || rb <= lb + 1)
	{
		throw logic_error(
			"ThreeWayPerceiverSenderAgreeEffect: cannot parse perceiver index "
			"from variable name '" + varName + "'. Expected format 'Y[n]'.");
	}

	// "Y[shared]" is the canonical (renamed) copy of Y[1] produced by the
	// R-side sharedCov mechanism.  Treat it as perceiver 0 (0-based).
	// All other slices carry their 1-based numeric index, e.g. "Y[2]" → 1.
	string innerStr = varName.substr(lb + 1, rb - lb - 1);
	if (innerStr == "shared")
	{
		this->lperceiverIndex = 0;
	}
	else
	{
		this->lperceiverIndex = stoi(innerStr) - 1;   // 1-based → 0-based
	}

	// Derive and look up the self-reported network name, e.g. "Y[self]".
	string selfName = varName.substr(0, lb) + "[self]";
	this->lpSelfNetwork = pState->pNetwork(selfName);

	if (!this->lpSelfNetwork)
	{
		throw logic_error(
			"ThreeWayPerceiverSenderAgreeEffect: self-reported network '"
			+ selfName + "' not found in state. "
			"Ensure Y[self] is included in the sienaData object.");
	}
}


/**
 * Calculates the contribution of flipping the tie (ego → alter) in Y[i]
 * to the statistic.
 *
 * Contribution = x_{ego,alter}(Y[self]) * (I(v_perceiver == v_ego) - mu_sim)
 *
 * ego        = j  (sender whose decision is evaluated)
 * alter      = k  (receiver of the potential tie)
 * perceiver  = i  (slice owner, stored in lperceiverIndex)
 */
double ThreeWayPerceiverSenderAgreeEffect::calculateContribution(int alter) const
{
	// x_{jjk}: j's self-reported tie to k.
	double selfTie =
		static_cast<double>(this->lpSelfNetwork->tieValue(this->ego(), alter));

	// (I(v_i == v_j) - mu_sim): centered similarity between perceiver and ego.
	double sim =
		this->actor_similarity(this->lperceiverIndex, this->ego());

	return selfTie * sim;
}


/**
 * Contribution of the existing tie (ego → alter) to the network statistic,
 * with missingness guards applied.
 */
double ThreeWayPerceiverSenderAgreeEffect::tieStatistic(int alter)
{
	if (this->missing(this->lperceiverIndex) || this->missing(this->ego()))
	{
		return 0;
	}
	return this->calculateContribution(alter);
}

}
