#pragma once

class C_AbsVertex: public C_AbsDiagram{
	protected:
	C_AbstractKernel * Kernel;
	C_DedicMem_Abs * Memory_abs;
	
	public:
	void info(){
		cout << "AbsVertex" << endl;
	}
	
	virtual void Initialization() = 0;
	//virtual dcx getDressingAt(dcx point) = 0;
	
};
/*
enum VertexState_ID { Initialized=0, Dressed, CauchyDressed, VertexState_ID_End };

class C_VertexState {
	private:
	bool flag_
};
*/

