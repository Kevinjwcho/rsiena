/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: PerceiverRestrictedInTieFunction.cpp
 *
 * Description: Implementation of PerceiverRestrictedInTieFunction.
 *****************************************************************************/

#include <string>
#include <stdexcept>

#include "PerceiverRestrictedInTieFunction.h"
#include "model/tables/NetworkCache.h"

using namespace std;

namespace siena
{

PerceiverRestrictedInTieFunction::PerceiverRestrictedInTieFunction(
		string interactionName, string outcomeName) :
	OneModeNetworkAlterFunction(interactionName),
	lperceiverIndex(-1)
{
	// Parse "Y[k]" from the OUTCOME variable name -> k - 1.  The
	// canonical "Y[shared]" name is already remapped to "Y[1]" before
	// reaching the C++ side (initializeFRAN.r:316), so we only ever see
	// numeric slice tags.  Anything else leaves lperceiverIndex at -1
	// and value() returns 0 unconditionally (safe no-op).
	string::size_type lb = outcomeName.find('[');
	string::size_type rb = outcomeName.find(']');
	if (lb != string::npos && rb != string::npos && rb > lb + 1)
	{
		string innerStr = outcomeName.substr(lb + 1, rb - lb - 1);
		if (innerStr != "shared")
		{
			try
			{
				lperceiverIndex = stoi(innerStr) - 1;   // 1-based -> 0-based
			}
			catch (const exception &)
			{
				lperceiverIndex = -1;
			}
		}
	}
}

/**
 * Returns Y[self]_{alter, ego} ONLY when alter equals this slice's
 * perceiver index (parsed from the outcome variable name).  Otherwise
 * returns 0.  This implements the perceiver-restriction: only tie flips
 * whose alter IS the perceiver of the current perceived slice contribute,
 * and they contribute the perceiver's own back-tie value from Y[self].
 */
double PerceiverRestrictedInTieFunction::value(int alter) const
{
	if (alter != this->lperceiverIndex)
	{
		return 0.0;
	}
	// inTieValue(alter) on the interaction1 network (Y[self]) returns
	// Y[self]_{alter, ego}.  Since alter == perceiverIndex here, this is
	// exactly Y[self]_{perceiver, ego} = perceiver's outgoing tie to the
	// actor doing the toggle in Y[i] — the term the user's statistic
	// asks for.
	return this->pNetworkCache()->inTieValue(alter);
}

}
