/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ThreeWayPerceiverEffect.cpp
 *
 * Description: This file contains the implementation of the
 * ThreeWayPerceiverEffect class.
 *****************************************************************************/

#include <string>
#include <stdexcept>

#include "ThreeWayPerceiverEffect.h"
#include "model/EffectInfo.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo  the effect descriptor
 * @param[in] mode         which perceiver-attribute statistic to compute
 */
ThreeWayPerceiverEffect::ThreeWayPerceiverEffect(
		const EffectInfo * pEffectInfo,
		ThreeWayPerceiverMode mode) :
	CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lmode = mode;
	this->lperceiverIndex = -1;
}


/**
 * Initializes this effect.
 *
 * Calls the parent initializer (which resolves interactionName1 to a
 * covariate), then parses the perceiver index from the variable name.
 * The variable name must have the form "Y[n]" where n is the 1-based
 * perceiver actor number (as assigned on the R side).
 *
 * @param[in] pData    the observed data
 * @param[in] pState   state of dependent variables at the start of the period
 * @param[in] period   the current period
 * @param[in] pCache   cache object for speed-up calculations
 */
void ThreeWayPerceiverEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	CovariateDependentNetworkEffect::initialize(pData, pState, period, pCache);

	// Parse "n" from variable name of the form "Y[n]" and convert to 0-based.
	string varName = this->pEffectInfo()->variableName();
	string::size_type lb = varName.find('[');
	string::size_type rb = varName.find(']');

	if (lb == string::npos || rb == string::npos || rb <= lb + 1)
	{
		throw logic_error(
			"ThreeWayPerceiverEffect: cannot parse perceiver index from "
			"variable name '" + varName + "'. Expected format 'Y[n]'.");
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
}


/**
 * Calculates the contribution of a potential tie flip (ego → alter) to
 * the statistic, assuming the perceiver's slice is being updated.
 *
 * ego   = j  (the sender whose decision is evaluated)
 * alter = k  (the receiver of the potential tie)
 * perceiver = i  (the slice owner, stored in lperceiverIndex)
 */
double ThreeWayPerceiverEffect::calculateContribution(int alter) const
{
	switch (this->lmode)
	{
	case TW_PERC_EGO:
		// percX:  v_i - mean(v)
		// Independent of alter and ego; same for every tie in this slice.
		return this->value(this->lperceiverIndex);

	case TW_PERC_ALT_SAME:
		// percAltSame:  I(v_i == v_k) - mu_sim
		return this->actor_similarity(this->lperceiverIndex, alter);

	case TW_PERC_SEND_SAME:
		// percSenderSame:  I(v_i == v_j) - mu_sim
		// Constant across k for a given ego j.
		return this->actor_similarity(this->lperceiverIndex, this->ego());

	default:
		return 0;
	}
}


/**
 * Contribution of the existing tie (ego → alter) to the network statistic,
 * with missingness guards applied.
 */
double ThreeWayPerceiverEffect::tieStatistic(int alter)
{
	// Always guard on the perceiver's covariate being observed.
	if (this->missing(this->lperceiverIndex))
	{
		return 0;
	}

	// Additional guard depending on which actor's value enters the formula.
	switch (this->lmode)
	{
	case TW_PERC_ALT_SAME:
		if (this->missing(alter))
		{
			return 0;
		}
		break;

	case TW_PERC_SEND_SAME:
		if (this->missing(this->ego()))
		{
			return 0;
		}
		break;

	default:
		break;
	}

	return this->calculateContribution(alter);
}

}
