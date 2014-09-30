#pragma once

class C_AbsVertex: public C_AbsDiagram{
	public:
	void info(){
		std::cout << "AbsVertex" << std::endl;
	}
	
	virtual void Initialization() = 0;
	//virtual dcx getDressingAt(dcx point) = 0;
	
};

