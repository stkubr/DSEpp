#pragma once

#include <iostream>
#include <source/types.h>

class C_Integration_Mockup
{
	private :
	int intNumIntegrationPoints, intFunctionType;
	double dDownLimit, dUpLimit;
	double dIntegrationResult;
	double (C_Integration_Mockup::*TestingFunction) (double);
	
	
	public:
	C_Integration_Mockup(){
		intNumIntegrationPoints=1000;
		dDownLimit=-1.0;
		dUpLimit=1.0;
	}
	
	void setFunctionType(int __intFunctionType){
		intFunctionType=__intFunctionType;
		switch (intFunctionType)
		{
			case 0: 
				TestingFunction=&C_Integration_Mockup::FunctionSine3x;
				break;
			case 1: 
				TestingFunction=&C_Integration_Mockup::FunctionCheb;
				break;
			case 2: 
				TestingFunction=&C_Integration_Mockup::FunctionPower;
				break;
			case 3: 
				TestingFunction=&C_Integration_Mockup::FunctionChebCheb;
				break;
			case 4: 
				TestingFunction=&C_Integration_Mockup::FunctionPowerCheb;
				break;
			default:
				TestingFunction=&C_Integration_Mockup::FunctionLine;
				std::cout << "Function Line selected" << std::endl;
		}
	}
	
	t_dMatrix IntegrandMockup(double x){
		t_dMatrix result(1,1);
		result(0,0)=(this->*TestingFunction)(x);
		return result;
	}
	
	double FunctionLine(double x){
		return x;
	}
	
	double FunctionSine3x(double x){
		return sin(3.0*x);
	}
	
	double FunctionCheb(double x){
		return sqrt(1-x*x);
	}
	
	double FunctionPower(double x){
		return cosh(sinh(x));
	}
	
	double FunctionChebCheb(double x){
		return 1.0;
	}
	
	double FunctionPowerCheb(double x){
		return cosh(sinh(x))/sqrt(1-x*x);
	}
	
	int getNumIntegrationPoints(){
		return intNumIntegrationPoints;
	}
	
	double getDownLimit(){
		return dDownLimit;
	}
	
	double getUpLimit(){
		return dUpLimit;
	}
};
