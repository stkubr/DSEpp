//
// Created by stkubr on 06.04.15.
//

#ifndef _DSEPP_INTEGRATOR_PATH_HPP_
#define _DSEPP_INTEGRATOR_PATH_HPP_

template <typename T_out,typename T_contour ,typename T_in> class C_Integrator_Path{
protected:
    int NumAmps;
    std::function <T_in(T_in, T_in)> * IntegrationWeightFunction;

    C_Integrator_Path(int _NumAps):NumAmps(_NumAps){}

    void setIntegrationWeightFunction(std::function<T_in(T_in, T_in)> *_IntegrationWeightFunction){
        IntegrationWeightFunction = _IntegrationWeightFunction;
    }

public:
    static C_Integrator_Path * createIntegrator(int _NumPoints, std::function <T_in(T_in, T_in)> * _IntegrationWeightFunction ){
        C_Integrator_Path * p;
        p = new C_Integrator_Path(_NumPoints);
        p -> setIntegrationWeightFunction(_IntegrationWeightFunction);
        return p;
    }

    vector<T_out> getResult(T_contour &PathStorage, T_in &Point)
    {
        vector<T_out> result(NumAmps,0);
        vector<T_in> sum(NumAmps,0);
        T_in sumNorm = 0;
        T_in IntegrationWeight;

        for (int j=0;j< PathStorage[0].size();j++){
            IntegrationWeight = (*IntegrationWeightFunction)(PathStorage[0][j], Point) * PathStorage[1][j];
            for(int i = 0; i < NumAmps; i++) sum[i] += IntegrationWeight * PathStorage[i+2][j];
            sumNorm += IntegrationWeight;
        }

        for(int i = 0; i < NumAmps; i++) result[i] = sum[i]/sumNorm;
        return result;
    }
};

#endif //_DSEPP_INTEGRATOR_PATH_HPP_
