/*
 * DedicMem.cpp
 *
 *  Created on: Jan 24, 2014
 *      Author: stkubr
 */

#include "DedicMem.hpp"


	void C_DedicMem_Quark::resizeGrid(int amps_num, int exter_num, int inter_num){
		S_grid.resize(amps_num,t_cmplxArray2D(exter_num,t_cmplxArray1D(inter_num)));
	}
	void C_DedicMem_Quark::resizeContour(int amps_num, int cont_num){
		S_cont.resize(amps_num, t_cmplxArray1D(cont_num));
	}

	void C_DedicMem_Quark::RemoveGrid(){
		S_grid.clear();
		std::cout << "Grid storage has been erased."<< std::endl;
	}

	void C_DedicMem_Quark::info(){
		std::cout << "Dedicated Memory for Quark" << std::endl;
	}




	t_cmplxArray4D VertexDressings;
	void C_DedicMem_Kernel::info(){
		std::cout << "Dedicated Memory for Kernel" << std::endl;
	}

	void C_DedicMem_Kernel::ResizeKstorage(int i){
		K_Matrix_Storage.resize(i);
	}

	void C_DedicMem_Kernel::SetKmatrixAt(int i ,t_cmplxMatrix2D * _K){
		K_Matrix_Storage[i]=(*_K);
	}

	void C_DedicMem_Kernel::EraseKstorage(){
			K_Matrix_Storage.clear();
	}

	t_cmplxMatrix2D * C_DedicMem_Kernel::GetKmatrixAt(int i){
		return &K_Matrix_Storage[i];
	}



	// AmpStorage Manipulation
//----------------------------------------------------------------------
	void C_DedicMem_BSA::resizeAmpStorage(int amps, int points){
		AmpStorage.resize(amps,std::vector <t_cmplxDirac>(points));
	}

	void C_DedicMem_BSA::setAmpStorage(int amp, int point, t_cmplxDirac * Amp){
		AmpStorage[amp][point]=(*Amp);
	}

	t_cmplxDirac C_DedicMem_BSA::getAmpStorage(int amp, int point){
		return AmpStorage[amp][point];
	}

	void C_DedicMem_BSA::clearAmpStorage(){
		AmpStorage.clear();
	}
//----------------------------------------------------------------------


	//EVMatrix Manipulation
//----------------------------------------------------------------------
	void C_DedicMem_BSA::ResizeEVMatrix(int num_rad, int num_angle, int num_amplitudes, int num_chebs)
	{
		EVMatrix = Eigen::MatrixXcf::Ones(num_rad*num_angle*num_amplitudes*num_chebs,num_rad*num_angle*num_amplitudes*num_chebs);
	}
//----------------------------------------------------------------------

	//BSE_onCauchy Contour ang Grid Manipulation
//----------------------------------------------------------------------
	void C_DedicMem_BSA::resizeBSEContour(int num_amplitudes, int num_points){
		CauchyContour.resize(num_amplitudes,t_cmplxArray1D(num_points));
	}

	void C_DedicMem_BSA::resizeBSEGrid(int num_amplitudes, int num_ex, int num_in){
		CauchyGrid.resize(num_amplitudes, t_cmplxArray2D(num_ex, t_cmplxArray1D(num_in)));
	}

	void C_DedicMem_BSA::clearCauchyGrid(){
		CauchyGrid.clear();
	}

//----------------------------------------------------------------------

	void C_DedicMem_BSA::info(){
		std::cout << "Dedicated Memory for BSE" << std::endl;
	}


	void C_DedicMem_1_LoopDiagram::ErasePathes()
	{
		Path_Photon.clear();
		Path_Pion_p.clear();
		Path_Pion_m.clear();
		Path_Quark_p_m.clear();
		Path_Quark_p_p.clear();
		Path_Quark_m_m.clear();
	}

	void C_DedicMem_1_LoopDiagram::info(){
		std::cout << "Dedicated Memory for 1 Loop Diagram" << std::endl;
	}



	void C_DedicMem_2_LoopDiagram::SetNumPathesAndStorages(int num_pathes, int num_storages){
		Pathes.resize(num_pathes);
		Storages.resize(num_storages);
	}

	void C_DedicMem_2_LoopDiagram::ErasePathesAndStorages()
	{
		Pathes.clear();
		Storages.clear();
	}

	void C_DedicMem_2_LoopDiagram::info(){
		std::cout << "Dedicated Memory for 1 Loop Diagram" << std::endl;
	}



	void C_DedicMem_2_Loop_Int::ResizePointsAndWieghts(int total_points, int dim){
		PointsAndWieghts.resize(2,t_dArray2D(total_points,t_dArray1D(dim)));
		Counters.resize(3,std::vector<int>(total_points));
	}

	void C_DedicMem_2_Loop_Int::ErasePointsAndWieghts(){
		PointsAndWieghts.clear();
	}

	void C_DedicMem_2_Loop_Int::info(){
		std::cout << "Dedicated Memory for 1 Loop Diagram" << std::endl;
	}


	void C_Memory_Manager::CopyMemoryFrom(C_DedicMem_Quark * quark_memory, C_DedicMem_Kernel * kernel_memory){
		int j,num_steps;//,sumj;sumj=0;
		t_cmplx temp;
		num_steps=(quark_memory->S_cont)[0].size();
		t_cmplxArray2D * KernelQuarkStorage;
		KernelQuarkStorage=&(kernel_memory->VertexDressings[0][0]);
		//std::cout << "after copy" << "  " << (kernel_memory->VertexDressings.size()) << std::endl;
		(*KernelQuarkStorage).resize(4,t_cmplxArray1D(4*num_steps));
		for (j=1;j<=num_steps;j++)
		{
			(*KernelQuarkStorage)[0][j-1]=0.0;
			(*KernelQuarkStorage)[1][j-1]=(quark_memory->S_cont)[0][j-1];
			(*KernelQuarkStorage)[2][j-1]=(quark_memory->S_cont)[1][j-1];
			(*KernelQuarkStorage)[3][j-1]=(quark_memory->S_cont)[3][j-1];
		}
	}

	void C_Memory_Manager::CopyMemoryFrom(ifstream * DressingStream, C_DedicMem_Kernel * kernel_memory){
		int j,num_steps;//,sumj;sumj=0;
		t_cmplx temp;
		//num_steps=kernel_memory->VertexDressings[0][0][0].size();
		t_cmplxArray2D * KernelQuarkStorage;

		//KernelQuarkStorage=&(kernel_memory->VertexDressings[1][0]);
		//std::cout << "after copy" << "  " << (kernel_memory->VertexDressings.size()) << std::endl;
		//(*KernelQuarkStorage).resize(3,dcxArray1D(num_steps));

		int num_P, num_p;
		if ((*DressingStream).is_open()){
			(*DressingStream) >> num_P >> num_p;
			//kernel_memory->VertexDressings[1].resize(num_P,dcxArray2D(7,dcxArray1D(num_P*num_p)));
			kernel_memory->VertexDressings[1].resize(num_P+1);//,dcxArray2D(7,dcxArray1D(10*10)));
			std::cout << num_P <<"  "<< num_p <<"  "<< kernel_memory->VertexDressings[1].size() << std::endl;
			for (int i = 0; i < num_P+1; i++){
				KernelQuarkStorage=&(kernel_memory->VertexDressings[1][i]);
				(*KernelQuarkStorage).resize(7,t_cmplxArray1D(num_p));
				//std::cout << kernel_memory->VertexDressings.size() << "  " << kernel_memory->VertexDressings[1].size() << "  " << kernel_memory->VertexDressings[1][0].size() << "  " << kernel_memory->VertexDressings[1][0][0].size() << std::endl;
				for (j=0;j<num_p;j++)
				{
					for (int k = 0; k < 7; k++) (*DressingStream) >> (*KernelQuarkStorage)[k][j];
					//for (int k = 0; k < 7; k++) std::cout << (*KernelQuarkStorage)[k][j];
					//std::cout << std::endl;
					//(*DressingStream) >> (*KernelQuarkStorage)[0][j] >> (*KernelQuarkStorage)[1][j] >> (*KernelQuarkStorage)[2][j];
					//std::cout << num_steps << "  " << (*KernelQuarkStorage)[0][j] << "  " << (*KernelQuarkStorage)[1][j] << "  " << (*KernelQuarkStorage)[2][j] << "  " << j << std::endl;
				}
				//cin.get();
			}
		}
		else std::cout << "Cant open file!" << std::endl;
	}


C_MemoryFactory_Quark * DedicMemFactory_Quark = new C_MemoryFactory_Quark;
C_MemoryFactory_Kernel * DedicMemFactory_Kernel = new C_MemoryFactory_Kernel;
/*C_MemoryFactory_BSA * DedicMemFactory_BSA = new C_MemoryFactory_BSA;
C_MemoryFactory_1_LoopDiagram * DedicMemFactory_1_LoopDiagram = new C_MemoryFactory_1_LoopDiagram;
C_MemoryFactory_2_LoopDiagram * DedicMemFactory_2_LoopDiagram = new C_MemoryFactory_2_LoopDiagram;
C_MemoryFactory_2_Loop_Int * DedicMemFactory_2_Loop_Int = new C_MemoryFactory_2_Loop_Int;*/
C_Memory_Manager * MemoryManager= new C_Memory_Manager;

