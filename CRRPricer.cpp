#include <iostream>
#include "EuropeanVanillaOption.h"
#include "AmericanOption.h"

#include "BinaryTree.h"
#include "CRRPricer.h"

#include <cmath>
// Standard CRR Constructor
CRRPricer::CRRPricer(Option* o, int depth, double S0, double U, double D, double R) : _option(o), _depth(depth), _asset_price(S0), _up(U), _down(D),  _interest_rate (R) {
    if (!(_down<_interest_rate && _interest_rate<_up)) {throw std::invalid_argument("Arbitrage");} // Arbitrage check: The model admits no arbitrage if D < R < U
    if(o->isAsianOption()) {throw std::invalid_argument("Error : Asian option !");} // Constraint: CRR does not support Asian options (Path Dependent)

    _tree.setDepth(depth);
    _exerciseTree.setDepth(depth);
}
// Overloaded constructor
CRRPricer::CRRPricer(Option* option, int depth, double asset_price, double r, double volatility) : _asset_price(asset_price), _option(option), _depth(depth) {
    if (option->isAsianOption()) {throw std::invalid_argument("Error: Asian Option not supported");}
    else {
        // Calculate U, D, R based on Black-Scholes volatility and risk-free rate
        calculateParamsBS(depth, r, volatility);
        // Verify arbitrage condition after calculation
        if (!(_down<_interest_rate && _interest_rate<_up)) {throw std::invalid_argument("Arbitrage");} 

        _tree.setDepth(depth);
        _exerciseTree.setDepth(depth);
    }
}
CRRPricer::~CRRPricer() { }

void CRRPricer::DisplayTree() {
    _tree.display();
}

void CRRPricer::OptionDisplayTree() {
    _H.display();
}

void CRRPricer::compute() { //Fills the H tree with the option prices
    _H.setDepth(_depth);

    if (_depth==1) { //Case where depth = 1
        _H.setNode(0,0,_option->payoff(_tree.getNode(0,0)));
    }

    else{
        double q = (_interest_rate - _down)/(_up - _down);
        // Build the Stock Price Tree | S(n, i) = S0 * (1+U)^i * (1+D)^(n-i)
        _tree.setNode(0, 0, _asset_price);
        for (int n = 1; n <= _depth; n++) {

            for (int i = 0; i <= n; i++) {
                double S = _asset_price * pow(1 + _up, i) * pow(1 + _down, n - i);
                _tree.setNode(n, i, S);
            }
        }
        // Initialize Terminal Option Values at Maturity (n = N) | H(N, i) = payoff(S(N, i))
        for(int i = 0; i <= _depth; i++) {
            double S = _tree.getNode(_depth, i); 
            _H.setNode(_depth, i, _option->payoff(S));
            _exerciseTree.setNode(_depth, i, true); // forced to exercise at maturity
        }
        // Backward Induction
        if (_option->isAmericanOption()) {
            for (int i = _depth - 1; i >= 0; i--) {//go back up in the tree in order to calculate option's value at each node
                for (int j = 0; j <= i; j++) {
                    double continuous_value = (q * _H.getNode(i+1, j+1) + (1-q) * _H.getNode(i+1, j)) / (1+ _interest_rate);
                    double spot_price_at_node = _tree.getNode(i, j);
                    double intrinsic_value = _option->payoff(spot_price_at_node); 
                    _H.setNode(i, j, std::max(continuous_value, intrinsic_value));

                    // Set the exercise condtion for American options
                    _exerciseTree.setNode(i, j, continuous_value <= intrinsic_value);
                }
            }
        }
        // European options
        else {
            for (int i = _depth - 1; i >= 0; i--) {//go back up the tree by calculating the value of the option at each node
                for (int j = 0; j <= i; j++) {
                    double continuous_value = (q * _H.getNode(i+1, j+1) + (1-q) * _H.getNode(i+1, j)) / (1+ _interest_rate);
                    _H.setNode(i, j, continuous_value);
                    _exerciseTree.setNode(i, j, false);
                }
            }
        }
    }
}
    
double CRRPricer::operator()(bool closed_form) {
    if (!closed_form) { //Closed_form is false so we use compute() and return H(0,0)
        compute();
        return _H.getNode(0, 0);
    } 
    
    else { //Closed_form is true so we use the formula
        if (_option->isAmericanOption()) {
			throw std::invalid_argument("Closed formula not available for American options");
		}
        compute();
        double q = (_interest_rate - _down)/(_up - _down);
            double sum = 0;
            int N = _depth;
            
            for (int i = 0; i <= N; i++) {
                double combination_coeff = tgamma(N + 1.0) / 
                                        (tgamma(i + 1.0) * tgamma(N - i + 1.0));
                double probability = combination_coeff * pow(q, i) * pow(1.0 - q, N - i);
                double payoff = _option->payoff(_tree.getNode(N, i)); 
                sum += probability * payoff;
            }
            // Discount back to time 0
            return sum / pow(1.0 + _interest_rate, N);
        }
}

// Computes U, D, and R to approximate Black-Scholes dynamics
void CRRPricer::calculateParamsBS(int depth, double r, double volatility) {
    double h = _option->getExpiry() / depth;
    _up = std::exp((r + (std::pow(volatility, 2) / 2)) * h + volatility * std::sqrt(h)) - 1;
    _down = std::exp((r + (std::pow(volatility, 2) / 2)) * h - volatility * std::sqrt(h)) - 1;
    _interest_rate = std::exp(r * h) - 1;
}

//Indicate the exercise condition for american options true/false
bool CRRPricer::getExercise(int n, int i) {
	 return _exerciseTree.getNode(n, i); 
}