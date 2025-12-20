#pragma once
#include "Option.h"
#include <vector>

// Base class for Asian options, derived from Option
class AsianOption : public Option {   //Asian Option class derived from Option
private:
	std::vector<double> _timeSteps; //Timesteps of the AsianOption : (t1,t2,....,tm)
    double _strike;

public:
    AsianOption(double expiry, const std::vector<double>& monitoringTimes); //Constructor that takes the expiry date and monitoring times
    ~AsianOption(); //Destructor


    std::vector<double> getTimeSteps() const override;                      // Getter for time steps
    double payoffPath(const std::vector<double>&) const override;           // Payoff from a full price path
    bool isAsianOption() const override;                                    // Returns true

    double getStrike() const override;                                      // Getter for strike
    virtual OptionType getOptionType() const = 0;                           // Call or Put
    virtual double payoff(double) const = 0;                                // Payoff from average price
};
