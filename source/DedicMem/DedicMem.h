#pragma once

class C_Memory_Manager;

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
	void resizeGrid(int num_storages, int amps_num, int exter_num, int inter_num){
		S_grid.resize(num_storages,t_cmplxArray3D(amps_num,t_cmplxArray2D(exter_num,t_cmplxArray1D(inter_num))));
	}
	void resizeContour(int num_storages, int amps_num, int cont_num){
		S_cont.resize(num_storages, t_cmplxArray2D(amps_num, t_cmplxArray1D(cont_num)));
	}
	
	void RemoveGrid(){
		S_grid.clear();
		std::cout << "Grid storage has been erased."<< std::endl;
	}

	void info(){
		std::cout << "Dedicated Memory for Quark" << std::endl;
	}	
};

class C_DedicMem_Kernel: public C_DedicMem_Abs {
	
	private:
	std::vector<t_cmplxMatrix2D > K_Matrix_Storage;
	
	public:
	t_cmplxArray4D VertexDressings;
	void info(){
		std::cout << "Dedicated Memory for Kernel" << std::endl;
	}

	void ResizeKstorage(int i){
		K_Matrix_Storage.resize(i);
	}
	
	void SetKmatrixAt(int i ,t_cmplxMatrix2D * _K){
		K_Matrix_Storage[i]=(*_K);
	}
	
	void EraseKstorage(){
			K_Matrix_Storage.clear();
	}
	
	t_cmplxMatrix2D * GetKmatrixAt(int i){
		return &K_Matrix_Storage[i];
	}
	
	
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
	void resizeAmpStorage(int amps, int points){
		AmpStorage.resize(amps,std::vector <t_cmplxDirac>(points));
	}
	
	void setAmpStorage(int amp, int point, t_cmplxDirac * Amp){
		AmpStorage[amp][point]=(*Amp);
	}
	
	t_cmplxDirac getAmpStorage(int amp, int point){
		return AmpStorage[amp][point];
	}
	
	void clearAmpStorage(){
		AmpStorage.clear();
	}
//----------------------------------------------------------------------


	//EVMatrix Manipulation
//----------------------------------------------------------------------
	void ResizeEVMatrix(int num_rad, int num_angle, int num_amplitudes, int num_chebs)
	{
		EVMatrix = Eigen::MatrixXcf::Ones(num_rad*num_angle*num_amplitudes*num_chebs,num_rad*num_angle*num_amplitudes*num_chebs);
	}
//----------------------------------------------------------------------
	
	//BSE_onCauchy Contour ang Grid Manipulation
//----------------------------------------------------------------------
	void resizeBSEContour(int num_amplitudes, int num_points){
		CauchyContour.resize(num_amplitudes,t_cmplxArray1D(num_points));
	}
	
	void resizeBSEGrid(int num_amplitudes, int num_ex, int num_in){
		CauchyGrid.resize(num_amplitudes, t_cmplxArray2D(num_ex, t_cmplxArray1D(num_in)));
	}
	
	void clearCauchyGrid(){
		CauchyGrid.clear();
	}
	
//----------------------------------------------------------------------
	
	void info(){
		std::cout << "Dedicated Memory for BSE" << std::endl;
	}
};

class C_DedicMem_1_LoopDiagram: public C_DedicMem_Abs{
	public:
	t_cmplxArray1D Path_Photon,Path_Pion_p,Path_Pion_m,Path_Quark_p_m,Path_Quark_p_p,Path_Quark_m_m;
	std::vector <t_cmplxMatrix> Photon_Stg,Pion_p_Stg,Pion_m_Stg,Quark_p_p_Stg,Quark_m_m_Stg,Quark_p_m_Stg;
	
	void ErasePathes()
	{
		Path_Photon.clear();
		Path_Pion_p.clear();
		Path_Pion_m.clear();
		Path_Quark_p_m.clear();
		Path_Quark_p_p.clear();
		Path_Quark_m_m.clear();
	}
	
	void info(){
		std::cout << "Dedicated Memory for 1 Loop Diagram" << std::endl;
	}
};

class C_DedicMem_2_LoopDiagram: public C_DedicMem_Abs{
	public:
	t_cmplxArray2D Pathes;
	std::vector< std::vector<t_cmplxMatrix> > Storages;
	
	void SetNumPathesAndStorages(int num_pathes, int num_storages){
		Pathes.resize(num_pathes);
		Storages.resize(num_storages);
	}
	
	void ErasePathesAndStorages()
	{
		Pathes.clear();
		Storages.clear();
	}
	
	void info(){
		std::cout << "Dedicated Memory for 1 Loop Diagram" << std::endl;
	}
};

class C_DedicMem_2_Loop_Int: public C_DedicMem_Abs{
	public:
	t_dArray3D PointsAndWieghts;
	std::vector<std::vector<int> > Counters;
	
	void ResizePointsAndWieghts(int total_points, int dim){
		PointsAndWieghts.resize(2,t_dArray2D(total_points,t_dArray1D(dim)));
		Counters.resize(3,std::vector<int>(total_points));
	}
	
	void ErasePointsAndWieghts(){
		PointsAndWieghts.clear();
	}
	
	void info(){
		std::cout << "Dedicated Memory for 1 Loop Diagram" << std::endl;
	}
};

class C_Memory_Manager{
	public:
	void CopyMemoryFrom(C_DedicMem_Quark * quark_memory, C_DedicMem_Kernel * kernel_memory){
		int j,num_steps;//,sumj;sumj=0;
		t_cmplx temp;
		num_steps=(quark_memory->S_cont)[0][0].size();
		t_cmplxArray2D * KernelQuarkStorage;
		KernelQuarkStorage=&(kernel_memory->VertexDressings[0][0]);
		//std::cout << "after copy" << "  " << (kernel_memory->VertexDressings.size()) << std::endl;
		(*KernelQuarkStorage).resize(4,t_cmplxArray1D(4*num_steps));
		for (j=1;j<=num_steps;j++) 
		{
			(*KernelQuarkStorage)[0][j-1]=0.0;
			(*KernelQuarkStorage)[1][j-1]=(quark_memory->S_cont)[0][0][j-1];
			(*KernelQuarkStorage)[2][j-1]=(quark_memory->S_cont)[0][1][j-1];
			(*KernelQuarkStorage)[3][j-1]=(quark_memory->S_cont)[0][3][j-1];
			
			(*KernelQuarkStorage)[0][j-1 + num_steps]=0.0;
			(*KernelQuarkStorage)[1][j-1 + num_steps]=(quark_memory->S_cont)[1][0][j-1];
			(*KernelQuarkStorage)[2][j-1 + num_steps]=(quark_memory->S_cont)[1][1][j-1];
			(*KernelQuarkStorage)[3][j-1 + num_steps]=(quark_memory->S_cont)[1][3][j-1];
			
			(*KernelQuarkStorage)[0][j-1 + 2*num_steps]=0.0;
			(*KernelQuarkStorage)[1][j-1 + 2*num_steps]=(quark_memory->S_cont)[2][0][j-1];
			(*KernelQuarkStorage)[2][j-1 + 2*num_steps]=(quark_memory->S_cont)[2][1][j-1];
			(*KernelQuarkStorage)[3][j-1 + 2*num_steps]=(quark_memory->S_cont)[2][3][j-1];
			//std::cout << -1.0*(quark_memory->S_cont)[2][3][j-1] << "  " << temp*(quark_memory->S_cont)[2][1][j-1] << std::endl;
			
			(*KernelQuarkStorage)[0][j-1 + 3*num_steps]=0.0;
			(*KernelQuarkStorage)[1][j-1 + 3*num_steps]=(quark_memory->S_cont)[3][0][j-1];
			(*KernelQuarkStorage)[2][j-1 + 3*num_steps]=(quark_memory->S_cont)[3][1][j-1];
			(*KernelQuarkStorage)[3][j-1 + 3*num_steps]=(quark_memory->S_cont)[3][3][j-1];
		}
	}
	
	void CopyMemoryFrom(ifstream * DressingStream, C_DedicMem_Kernel * kernel_memory){
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

C_MemoryFactory_Quark * DedicMemFactory_Quark = new C_MemoryFactory_Quark;
C_MemoryFactory_Kernel * DedicMemFactory_Kernel = new C_MemoryFactory_Kernel;
C_MemoryFactory_BSA * DedicMemFactory_BSA = new C_MemoryFactory_BSA;
C_MemoryFactory_1_LoopDiagram * DedicMemFactory_1_LoopDiagram = new C_MemoryFactory_1_LoopDiagram;
C_MemoryFactory_2_LoopDiagram * DedicMemFactory_2_LoopDiagram = new C_MemoryFactory_2_LoopDiagram;
C_MemoryFactory_2_Loop_Int * DedicMemFactory_2_Loop_Int = new C_MemoryFactory_2_Loop_Int;
C_Memory_Manager * MemoryManager= new C_Memory_Manager;
