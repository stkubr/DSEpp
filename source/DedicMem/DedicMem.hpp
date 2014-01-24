/*
 * DedicMem.hpp
 *
 *  Created on: Jan 21, 2014
 *      Author: stkubr
 */

#ifndef DEDICMEM_HPP_
#define DEDICMEM_HPP_

#include <vector>
#include <fstream>
#include "../eigen/Eigen/Eigenvalues"
#include "../types.h"

class C_DedicMem_Abs{
	
	public:
	virtual void info() = 0;
	virtual ~C_DedicMem_Abs() {}
    friend class C_Memory_Manager;
};

class C_DedicMem_Quark: public C_DedicMem_Abs {
	public:
	t_cmplxArray4D S_grid;
	t_cmplxArray3D S_cont;
	
	public:
	void resizeGrid(int num_storages, int amps_num, int exter_num, int inter_num);
	void resizeContour(int num_storages, int amps_num, int cont_num);
	
	void RemoveGrid();

	void info();
};

class C_DedicMem_Kernel: public C_DedicMem_Abs {
	
	private:
	std::vector<t_cmplxMatrix2D > K_Matrix_Storage;
	
	public:
	t_cmplxArray4D VertexDressings;
	void info();

	void ResizeKstorage(int i);
	
	void SetKmatrixAt(int i ,t_cmplxMatrix2D * _K);
	
	void EraseKstorage();
	
	t_cmplxMatrix2D * GetKmatrixAt(int i);
	
};

class C_DedicMem_BSA: public C_DedicMem_Abs {
	
	public:
	std::vector< std::vector <t_cmplxDirac> > AmpStorage;
	Eigen::MatrixXcf EVMatrix;
	t_cmplxArray2D CauchyContour;
	t_cmplxArray3D CauchyGrid;
	
	public:
	
	// AmpStorage Manipulation
//----------------------------------------------------------------------
	void resizeAmpStorage(int amps, int points);
	
	void setAmpStorage(int amp, int point, t_cmplxDirac * Amp);
	
	t_cmplxDirac getAmpStorage(int amp, int point);
	
	void clearAmpStorage();
//----------------------------------------------------------------------


	//EVMatrix Manipulation
//----------------------------------------------------------------------
	void ResizeEVMatrix(int num_rad, int num_angle, int num_amplitudes, int num_chebs);
//----------------------------------------------------------------------
	
	//BSE_onCauchy Contour ang Grid Manipulation
//----------------------------------------------------------------------
	void resizeBSEContour(int num_amplitudes, int num_points);
	
	void resizeBSEGrid(int num_amplitudes, int num_ex, int num_in);
	
	void clearCauchyGrid();
	
//----------------------------------------------------------------------
	
	void info();
};

class C_DedicMem_1_LoopDiagram: public C_DedicMem_Abs{
	public:
	t_cmplxArray1D Path_Photon,Path_Pion_p,Path_Pion_m,Path_Quark_p_m,Path_Quark_p_p,Path_Quark_m_m;
	std::vector <t_cmplxMatrix> Photon_Stg,Pion_p_Stg,Pion_m_Stg,Quark_p_p_Stg,Quark_m_m_Stg,Quark_p_m_Stg;
	
	void ErasePathes();
	
	void info();
};

class C_DedicMem_2_LoopDiagram: public C_DedicMem_Abs{
	public:
	t_cmplxArray2D Pathes;
	std::vector< std::vector<t_cmplxMatrix> > Storages;
	
	void SetNumPathesAndStorages(int num_pathes, int num_storages);
	
	void ErasePathesAndStorages();
	
	void info();
};

class C_DedicMem_2_Loop_Int: public C_DedicMem_Abs{
	public:
	t_dArray3D PointsAndWieghts;
	std::vector<std::vector<int> > Counters;
	
	void ResizePointsAndWieghts(int total_points, int dim);
	
	void ErasePointsAndWieghts();
	
	void info();
};


class C_Memory_Manager{
	public:
	void CopyMemoryFrom(C_DedicMem_Quark * quark_memory, C_DedicMem_Kernel * kernel_memory);
	
	void CopyMemoryFrom(ifstream * DressingStream, C_DedicMem_Kernel * kernel_memory);

};


	class C_MemoryFactory{
	public:
		virtual C_DedicMem_Abs* CreateMemory()=0;
		virtual ~C_MemoryFactory() {}
	};

	class C_MemoryFactory_Quark: public C_MemoryFactory{
		public:
		C_DedicMem_Abs* CreateMemory() {
			return new C_DedicMem_Quark;
	    }
	};

	class C_MemoryFactory_Kernel: public C_MemoryFactory{
		public:
		C_DedicMem_Abs* CreateMemory() {
			return new C_DedicMem_Kernel;
	    }
	};

	class C_MemoryFactory_BSA: public C_MemoryFactory{
		public:
		C_DedicMem_Abs* CreateMemory() {
			return new C_DedicMem_BSA;
	    }
	};

	class C_MemoryFactory_1_LoopDiagram: public C_MemoryFactory{
		public:
		C_DedicMem_Abs* CreateMemory() {
			return new C_DedicMem_1_LoopDiagram;
	    }
	};

	class C_MemoryFactory_2_LoopDiagram: public C_MemoryFactory{
		public:
		C_DedicMem_Abs* CreateMemory() {
			return new C_DedicMem_2_LoopDiagram;
	    }
	};

	class C_MemoryFactory_2_Loop_Int: public C_MemoryFactory{
		public:
		C_DedicMem_Abs* CreateMemory() {
			return new C_DedicMem_2_Loop_Int;
	    }
	};

	extern C_MemoryFactory_Quark * DedicMemFactory_Quark;
	extern C_MemoryFactory_Kernel * DedicMemFactory_Kernel;
	/*C_MemoryFactory_BSA * DedicMemFactory_BSA = new C_MemoryFactory_BSA;
	C_MemoryFactory_1_LoopDiagram * DedicMemFactory_1_LoopDiagram = new C_MemoryFactory_1_LoopDiagram;
	C_MemoryFactory_2_LoopDiagram * DedicMemFactory_2_LoopDiagram = new C_MemoryFactory_2_LoopDiagram;
	C_MemoryFactory_2_Loop_Int * DedicMemFactory_2_Loop_Int = new C_MemoryFactory_2_Loop_Int;*/
	extern C_Memory_Manager * MemoryManager;

#endif /* DEDICMEM_HPP_ */
