
#ifndef DSEPP_GRADIENDDESCENT_H
#define DSEPP_GRADIENDDESCENT_H

#include <functional>
#include <algorithm>
#include "../types.h"
#include "Extra_functions.h"

class C_Gradiend_Descent{
public:
    vector<double> Parameters;
    vector<double> Derivatives;
    std::function<double(t_dArray1D)> Cost_function;
    const double delta = 0.01;
    const double eps = 0.00001;
    t_dArray1D alpha;
    int NumberParameters;

    C_Gradiend_Descent(){
        alpha.resize(1);
        Parameters.resize(1);
        Derivatives.resize(1);
    }

    C_Gradiend_Descent(int __NumberParameters){
        NumberParameters= __NumberParameters;
        alpha.resize(NumberParameters,1.0);
        Parameters.resize(NumberParameters);
        Derivatives.resize(NumberParameters);
        setInitialParameters();
    }

    void setAlpha(t_dArray1D __alpha){
        alpha = __alpha;
    }

    void setCostFunction(std::function<double(t_dArray1D)> & __Cost_function){
        Cost_function = __Cost_function;
    }

    void setInitialParameters(){
        std::generate(Parameters.begin(),Parameters.end(), []{return 1.0;} );
    }

    void setInitialParameters(t_dArray1D _parameters){
        Parameters = _parameters;
    }

    void setDerivatives(){
        if(Cost_function){
            int p_ctr =0;
            for (auto parameter_i: Parameters){
                std::function<double(double)>  Wraped_Cost_function = std::bind(&C_Gradiend_Descent::CostFunctionWrapper, this,
                                                                                Parameters , std::placeholders::_1, p_ctr);
                Derivatives[p_ctr] = derivative(Wraped_Cost_function, parameter_i, delta);
                cout << "Derivative = " << Derivatives[p_ctr] << " at " << parameter_i << endl;
                p_ctr++;
            }
        } else { std::cout << "Cost funсtion is not set" << endl; assert(1);}
    }

    void changeAlphas(){
        int p_ctr =0;
        for (auto deriv: Derivatives){
            alpha[p_ctr] = 0.01/(fabs(deriv) + 0.01);
            p_ctr++;
        }
    }

    void changeParameters(){
        int p_ctr =0;
        for (auto parameter_i: Parameters){
            Parameters[p_ctr] = parameter_i - alpha[p_ctr]*Derivatives[p_ctr];
            p_ctr++;
        }
    }

    void minimizeCostFunction(){
        double cost_value = 1.0;
        while (fabs(cost_value) > eps){
            setDerivatives();
            changeAlphas();
            changeParameters();
            printParameters();
            cost_value = Cost_function(Parameters);
            cout << "Cost function = " << cost_value << endl;
        }
        printParameters();
    }


    void printParameters(){
        int p_ctr =0;
        for (auto parameter_i: Parameters){
            cout << "Parameter_" << p_ctr << " = " << parameter_i << endl;
            p_ctr++;
        }
    }

    double CostFunctionWrapper(t_dArray1D & Parameters, double WrapParameter, int WrapIndex){
        t_dArray1D TempParameters = Parameters;
        TempParameters[WrapIndex] = WrapParameter;
        return Cost_function(TempParameters);
    }
};

#endif //DSEPP_GRADIENDDESCENT_H
