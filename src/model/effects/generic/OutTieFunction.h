/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: OutTieFunction.h
 *
 * Description: This file contains the definition of the
 * OutTieFunction class.
 *****************************************************************************/


#ifndef OUTTIEFUNCTION_H_
#define OUTTIEFUNCTION_H_

#include "NetworkAlterFunction.h"

namespace siena
{

/**
 * Defines a function that returns the values of ties from the ego
 * to alters in a network with the given name.
 */
class OutTieFunction: public NetworkAlterFunction
{
public:
	OutTieFunction(std::string networkName);
	OutTieFunction(std::string networkName, const bool simulatedState);
				
	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);		
	virtual void initialize(const Data * pData,
		State * pState, State * pSimulatedState, int period, Cache * pCache);

	virtual double value(int alter) const;
};

/**
 * Centered version of OutTieFunction: returns the tie value minus a constant
 * fixed at the mean tie value of the network at its first observation. Used
 * for cross-network effects whose predictor is another dependent network, so
 * that the effect is not a second density effect (three-way Group 2).
 */
class CenteredOutTieFunction: public OutTieFunction
{
public:
	CenteredOutTieFunction(std::string networkName);
	virtual void initialize(const Data * pData,
		State * pState, int period, Cache * pCache);
	virtual void initialize(const Data * pData,
		State * pState, State * pSimulatedState, int period, Cache * pCache);
	virtual double value(int alter) const;
	double center() const { return this->lcenter; }
private:
	void computeCenter(const Data * pData);
	double lcenter {};
};

}

#endif /* OUTTIEFUNCTION_H_ */
