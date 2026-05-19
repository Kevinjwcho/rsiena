/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PerceiverRestrictedInTieFunction.h
 *
 * Description:
 *   Cross-network "perceiver-restricted in-tie" alter function used to
 *   implement the three-way perceived reciprocity statistic
 *
 *       s = Σ_i Σ_j  Y[self]_{i, j} · Y[i]_{j, i}
 *
 *   on the OBJECTIVE of each perceived slice Y[i] (Direction 1).
 *
 *   Layout:
 *     * Outcome network (the function's "outer" effect variable):
 *           "Y[k]" — a perceived slice, k = 1..K
 *       (Y[self] does not evolve in the three-way fork, so the
 *       statistic can't live on Y[self]'s objective.)
 *     * Cross network (the function's referenced network via
 *       interaction1):  "Y[self]"
 *
 *   For a tie flip (ego, alter) in Y[k]:
 *       if alter == k - 1   (alter == slice's perceiver, 0-based)
 *           return Y[self]_{alter, ego}    = perceiver's outgoing tie to ego
 *       else
 *           return 0
 *
 *   Aggregated across the K slices via the existing crprod-family share
 *   machinery (interaction1 = Y[self] is share-marked for shareParameters
 *   on threeway), all K rows are tied to one β and the moment vector
 *   sums to the user's target statistic above.
 *
 *   The class is a thin extension of OneModeNetworkAlterFunction
 *   (which already plumbs the cache for the interaction1 network).
 *   The perceiver index is parsed once at construction from the
 *   *outcome* variable name (passed in as a second argument), not from
 *   the function's own network name.
 *****************************************************************************/

#ifndef PERCEIVERRESTRICTEDINTIEFUNCTION_H_
#define PERCEIVERRESTRICTEDINTIEFUNCTION_H_

#include "OneModeNetworkAlterFunction.h"

namespace siena
{

class PerceiverRestrictedInTieFunction: public OneModeNetworkAlterFunction
{
public:
	PerceiverRestrictedInTieFunction(std::string interactionName,
		std::string outcomeName);

	virtual double value(int alter) const;

private:
	/// 0-based actor index of the slice owner (perceiver i), parsed once
	/// at construction from the outcome variable name "Y[k]" (-> k - 1).
	/// -1 indicates a parse failure; the function then returns 0 always.
	int lperceiverIndex;
};

}

#endif /* PERCEIVERRESTRICTEDINTIEFUNCTION_H_ */
