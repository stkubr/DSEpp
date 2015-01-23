/*
 * Quark.cpp
 *
 *  Created on: Jan 23, 2014
 *      Author: stkubr
 */

#include "Quark.hpp"

C_Quark::C_Quark(){
	Memory_abs=DedicMemFactory_Quark->CreateMemory();
	Memory=(C_DedicMem_Quark*)Memory_abs;
	flag_dressed=false;
	flag_renormalization =false;
	num_amplitudes=2;
    num_IntegDimentions = 2;
	kinematicFactor=1.0/(16.0*pi*pi*pi);
	threadloc_integr_inx.resize(omp_get_num_threads());
	threadloc_p_momenta_inx.resize(omp_get_num_threads());
}

/*t_cmplx getTensorExpression(t_cmplxVector& p){
	return p*p;
}*/

void C_Quark::ReadParameters(string & _ParamPath){
	params.ReadParameters(_ParamPath);
}

t_cmplx C_Quark::getDressingFactor(){
	return Z2;
}

void C_Quark::setContourApex(double M2){
	params.M2_contour=M2;
}

// Set parameters to initial values
//----------------------------------------------------------------------
void C_Quark::InitialState(){
	A_renorm=1.0;
	B_renorm=0.0;
	B_mu=0.0;
	Z2=1.0;
	eps=1.0;
	check_res=0.0;
	flag_dressed=false;
	flag_renormalization =false;

    CauchyIntegratonWeight_lambda = [](t_cmplx path_var, t_cmplx coordinate) { return 1.0/(path_var - coordinate); };

	ResizeMemory();
}

// Resize all storages (internal and external), also create side objects like (Integrators, Kernels and etc.)
//----------------------------------------------------------------------
void C_Quark::InitializateIntegrators(){
    Integrator_momentum =C_Integrator_Line<t_cmplxMatrix,C_OneLoopIntegrator,double>::createIntegrator(
			params.num_prop_steps, params.LimDk, params.LimUk, num_amplitudes, qgausleg_lin_ID);
	Integrator_angle_Z =C_Integrator_Line<t_cmplxMatrix,C_OneLoopIntegrator,double>::createIntegrator(
			params.num_angle_steps, params.LimDk, params.LimUk, num_amplitudes, qgauscheb_ID);
	Integrator_momentum_short =C_Integrator_Line<t_cmplxMatrix,C_Quark,double>::createIntegrator(
			params.num_cutoff_steps, params.LimDk, params.LimUk, num_amplitudes, qgausleg_sym_ID);
    Integrator_cauchy=C_Integrator_Path<t_cmplx,t_cmplxArray2D,t_cmplx>::createIntegrator(num_amplitudes, &CauchyIntegratonWeight_lambda);

	Integrator_momentum->getNodes(&zz_rad,&w_rad);
	Integrator_momentum_short->getNodes(&zz_line,&w_line);
	Integrator_angle_Z->getNodes(&zz_angle,&w_angle);
}

void C_Quark::ResizeMemory(){
	Memory->resizeContour(num_amplitudes+2, 2*params.num_prop_steps + 2*params.num_cutoff_steps);
	Memory->resizeGrid(num_amplitudes+1, 2*params.num_prop_steps + 2*params.num_cutoff_steps, params.num_prop_steps*params.num_angle_steps);
}

// Allocator copy of "this" (used in parallel sections)
//----------------------------------------------------------------------
C_Quark * C_Quark::MakeCopy(){
	return new C_Quark(*this);
}

// Set Cauchy contour
//----------------------------------------------------------------------
void C_Quark::setContour(){
	t_cmplx last_point;
	t_cmplxArray2D tempContour;
	last_point = zz_rad[zz_rad.size()-1]*zz_rad[zz_rad.size()-1] - params.M2_contour/4.0 + params.LimUk*params.LimDk;
	Geometry::C_ParabolaContour contour(ii*sqrt(params.M2_contour)/2.0,
										ii*sqrt(params.M2_contour),
										last_point);
	contour.setParabolaContour(zz_rad,w_rad,zz_line,w_line);
	tempContour = contour.getParabolaContour();
	Memory->S_cont[0]=tempContour[0];
	Memory->S_cont[1]=tempContour[1];

	setInitialAandB();

	ofstream countur_data;
	countur_data.open ("Data_files/countur_data.dat");
	for (int j = 0; j < Memory->S_cont[0].size(); j++){
		countur_data << real(Memory->S_cont[0][j]) <<"  "<< imag(Memory->S_cont[0][j]) << std::endl;
	}
	countur_data.close();
	std::cout << "Contour has been set! points - " << Memory->S_cont[0].size() << std::endl;
}

// Set Cauchy contour
//----------------------------------------------------------------------
void C_Quark::setGrid() {
	// set (p - k)^2 grid
	for (int i = 0; i < Memory->S_cont[0].size(); i++) {
		for (int j = 0; j < params.num_prop_steps; j++) {
			for (int k = 0; k < params.num_angle_steps; k++) {
				Memory->S_grid[0][i][params.num_angle_steps * j + k] =
						/* p^2 */	Memory->S_cont[0][i]
						/* k^2 */   + zz_rad[j + 1] * zz_rad[j + 1]
					   /* (pk) */   - 2.0* sqrt(Memory->S_cont[0][i]* zz_rad[j + 1]* zz_rad[j + 1])* zz_angle[k + 1];
			}
		}
	}
	//draw grid
	ofstream grid_data;
	grid_data.open("Data_files/grid_data.dat");
	int i = 0;
	for (int j = 0; j < Memory->S_grid[0][0].size(); j++) {
		grid_data << real(Memory->S_grid[0][i][j]) << "  "
				  << imag(Memory->S_grid[0][i][j]) << std::endl;
	}
	grid_data.close();
	std::cout << "Grid has been set! points - "
			<< Memory->S_grid[0][0].size() * Memory->S_grid[0].size() << std::endl << std::endl;
}

// Initial guesses for A and B
//----------------------------------------------------------------------
void C_Quark::setInitialAandB(){
	auto initA = [] (t_cmplx z) -> t_cmplx {return 1.0 + exp(-z*z); };
	auto initB = [] (t_cmplx z) -> t_cmplx {return exp(-z*z); };

	for (int i = 0; i < Memory->S_cont[0].size(); i++){
		Memory->S_cont[2][i]=(initA(Memory->S_cont[0][1]));
		Memory->S_cont[3][i]=(initB(Memory->S_cont[0][1]));
	}
}

// Evaluate Cauchy integral on contour, at certain point
//----------------------------------------------------------------------
t_cmplxArray1D C_Quark::getResultCauchyAt(t_cmplx coordin){
    return Integrator_cauchy->getResult(Memory->S_cont, coordin);
}

// Evaluate Cauchy integral on contour, obtain Propagator on grid
//----------------------------------------------------------------------
void C_Quark::CalcPropGrid(){
#pragma omp parallel
{//begin of parallel pragma
		t_cmplx coordin;
		C_Quark * quark_copy;
		quark_copy=MakeCopy();
		t_cmplxArray1D S_temp_storage(num_amplitudes);
#pragma omp for
		for (int i = 0; i < Memory->S_cont[0].size(); i++){
			for (int j = 0; j < Memory->S_grid[0][0].size(); j++){
				coordin=Memory->S_grid[0][i][j];
				if (real(coordin)<params.LimUk*params.LimUk*params.EffectiveCutoff){
					S_temp_storage= quark_copy->getResultCauchyAt(coordin);
					Memory->S_grid[1][i][j]=(S_temp_storage[0]);
					Memory->S_grid[2][i][j]=(S_temp_storage[1]);
				} else {
					Memory->S_grid[1][i][j]=Z2*(1.0*params.HeavyLight -(params.HeavyLight-1.0)*real(coordin)) ;
					Memory->S_grid[2][i][j]=Z2*(params.m0);
				}
			}
		}
		delete quark_copy;
}//end of parallel pragma
}

// Evaluate DSE integral on grid, obtain quark propagator on contour
//----------------------------------------------------------------------
void C_Quark::CalcPropCont(){
	int num_contour=Memory->S_cont[0].size();
	std::function<t_cmplxMatrix(t_cmplxArray1D)>  bound_member_fn = std::bind(&C_Quark::Integrand_numerical,
			this,
			std::placeholders::_1);

#pragma omp parallel
{// start of parallel
		t_cmplxMatrix Temp_return(num_amplitudes,1),Bare_term(num_amplitudes,1);
		Bare_term(0, 0) = Z2;
		Bare_term(1, 0) = params.m0 - B_renorm;
#pragma omp for
		for (int i = 0; i < num_contour/2; i++){ // Iterating over upper part only
			threadloc_p_momenta_inx[omp_get_thread_num()] = i;
			threadloc_integr_inx[omp_get_thread_num()]=0;

			Temp_return = Bare_term + MultiDimInt2D_wo_nested(&bound_member_fn, num_amplitudes, 1);

			Memory->S_cont[2][i] = Temp_return(0, 0);
			Memory->S_cont[3][i] = Temp_return(1, 0);

			// Propagator is symmetric, the lower part is just mirrored conjugation
			Memory->S_cont[2][num_contour - 1 - i] = conj(Temp_return(0, 0));
			Memory->S_cont[3][num_contour - 1 - i] = conj(Temp_return(1, 0));
		}
}// end of parallel
}

// Set k and p vectors for the Numerical Integrand
//----------------------------------------------------------------------
void C_Quark::setKinematic(t_cmplxVector& k, t_cmplxVector& p, t_cmplx x, t_cmplx y, t_cmplx z){
	k.SetP4(0.0,0.0,sqrt(1.0-z*z)*y,y*z);
	p.SetP4(0.0,0.0,0.0,sqrt(x));
}

// Numerical form of the Integrand (suitable for general Kernel)
//----------------------------------------------------------------------
t_cmplxMatrix C_Quark::Integrand_numerical (t_cmplxArray1D integVariables){
	t_cmplx _A,_B,_B2,y2,pk,epsilon;

	int index_p = threadloc_p_momenta_inx[omp_get_thread_num()];
	int integration_inx = threadloc_integr_inx[omp_get_thread_num()];

	t_cmplx x = Memory->S_cont[0][index_p];
	t_cmplx y = integVariables[0];
	t_cmplx z = integVariables[1];


	t_cmplxVector k,p;
	setKinematic(k,p,x,y,z);

	t_cmplxMatrix result(num_amplitudes,1);
	t_cmplxDirac S(4,4),Proj1(4,4),Proj2(4,4);
	t_cmplxVector pion_momenta;

	pion_momenta=(p-k/2.0);
	pk=Memory->S_grid[0][index_p][integration_inx];
	_A=Memory->S_grid[1][index_p][integration_inx];
	_B=Memory->S_grid[2][index_p][integration_inx];
	epsilon=y*y*y/(pk*_A*_A+_B*_B);

	S=-ii*((p-k)*Y)*_A + I*_B;

	Proj1=-ii*(p*Y)/(p*p);
	Proj2=I;

    result(0, 0) = kinematicFactor * epsilon * Kernel->TraceKernelWithoutStoring(Proj1, S, k, p, pion_momenta, true);
    result(1, 0) = kinematicFactor * epsilon * Kernel->TraceKernelWithoutStoring(Proj2, S, k, p, pion_momenta, false);

	threadloc_integr_inx[omp_get_thread_num()]++;
	return result;
}

// Calculation A,B,M,Sigma_V,Sigma_S at point in contour
//----------------------------------------------------------------------
t_cmplxArray1D C_Quark::getPropAt(t_cmplx q){
	t_cmplxArray1D tempvector;
	t_cmplxArray1D storage(5);
	tempvector= getResultCauchyAt(q);

	// A
	storage[0]=tempvector[0];
	// B
	storage[1]=tempvector[1];
	// M
	storage[2]=storage[1]/storage[0];
	// Sigma_V
	storage[3]=(storage[0]/(q*storage[0]*storage[0]+storage[1]*storage[1]));
	// Sigma_S
	storage[4]=(storage[1]/(q*storage[0]*storage[0]+storage[1]*storage[1]));
	return storage;
}

// Calculate consequently Grid and Contour until converge
//----------------------------------------------------------------------
void C_Quark::PropSetAndCheck(){
	if(params.ReCalcProp){
		// First few steps without normalization to speed up convergence
		for (int i = 0; i < 2 ; i++){
			CalcPropGrid();
			CalcPropCont();
			write_Prop_re(20);
		}
		// Switching on normalization
		flag_renormalization =true;
		check_res=1.0;
		eps=1.0;
		while (eps>params.Accuracy){
			CalcPropGrid();
			CalcPropCont();
			Renormalize();
			PropCheck(100);
			write_Prop_re(20);
		}
	}
}

void C_Quark::Renormalize(){
	if (flag_renormalization){
		A_renorm=real(getResultCauchyAt(params.mu * params.mu)[0])-Z2*1.0;
		Z2=1.0 - (A_renorm);

		B_mu=real(getResultCauchyAt(params.mu * params.mu)[1]);
		B_renorm+=real(getResultCauchyAt(params.mu * params.mu)[1])-params.m0;
	}
}

// Check convergence
//----------------------------------------------------------------------
void C_Quark::PropCheck(int s){
	double Pu,Pd,x,scale,dp;
	t_cmplxArray1D storage;
	double res=0;
	scale = s;
	Pd=0.01;
	Pu=params.LimUk*0.8;
	dp=pow(10,(log10(Pu/Pd)/scale));
	std::vector<double> A_check(scale),B_check(scale);
	x=Pd;
	for (int i = 0; i < scale; i++){
		storage= getResultCauchyAt(x * x);
		A_check[i]=real(storage[0]);
		B_check[i]=real(storage[1]);
		res+=A_check[i]*B_check[i];
		x*=dp;
	}
	eps=fabs(res - check_res)/fabs(res);
	std::cout << "Z2 - " <<"  "<< Z2 <<"  "<< "A_mu" <<"  "<< real(getResultCauchyAt(params.mu * params.mu)[0]) <<"  "<<
			"m_renorm - " <<"  "<< B_mu <<"  "<< "Accuracy - " <<"  "<< eps << std::endl;
	check_res=res;
}

// Initialization (Dressing) of the Propagator
//----------------------------------------------------------------------
void C_Quark::DressPropagator(){
	if (flag_dressed==false)
	{
		PrintLine('-');
		std::cout << "Start Dressing for " << name << " with quark mass -" <<"  "<< params.m0 <<"  "<< "and contour -" <<"  "<< params.M2_contour << std::endl;
		PrintLine('-');
		setContour();
		setGrid();
		//LoadPropCountour();
		PropSetAndCheck();
		Memory->RemoveGrid();
		write_Prop_re(100);
		//CheckExtrapolation();
	}
}

// Draw Propagator at real line
//----------------------------------------------------------------------
void C_Quark::write_Prop_re(int s){
	double Pu,Pd,x,scale,dp;
	t_cmplxArray1D storage;
	scale = s;
	Pd=0.01;
	Pu=params.LimUk;
	dp=pow(10,(log10(Pu/Pd)/scale));
	x=Pd;
	std::cout.precision(7);
	ofstream data_Prop_re;
	data_Prop_re.open ("Data_files/data_Prop_re_new.dat");
	std::cout << std::endl;
	std::cout << "P^2" <<"         "<< "A(p)" <<"        "<< "B(p)"  << std::endl;
	for (int i = 0; i < scale; i++){
		storage= getResultCauchyAt(x * x);
		data_Prop_re << x*x <<"  "<< real(storage[0]) <<"  "<< real(storage[1])  << std::endl;
		std::cout << x*x <<"  "<< real(storage[0]) <<"  "<< real(storage[1])  << std::endl;
		x*=dp;
	}
	std::cout << std::endl;
	data_Prop_re.close();
}

// Write Propagator to file
//----------------------------------------------------------------------
void C_Quark::SavePropCountour(){
	ofstream SavePropStream;
	SavePropStream.open(SavePropPath);
	(SavePropStream) << "Z2" << '\t' << Z2 << std::endl;
	for (int j=1;j<=Memory->S_cont[0].size();j++){
		(SavePropStream) << Memory->S_cont[0][j-1] << '\t' <<
				Memory->S_cont[2][j-1] << '\t' <<
				Memory->S_cont[3][j-1]  << std::endl;
	}
	SavePropStream.close();
	std::cout << "Propagator was saved." << std::endl;
}

// Read Propagator from file
//----------------------------------------------------------------------
void C_Quark::LoadPropCountour(){
	string dummy;
	if (!params.ReCalcProp){
		ifstream PropContourStream;
		PropContourStream.open(SavePropPath);
		PropContourStream >> dummy >> Z2;
		if (PropContourStream.is_open()){
			for (int j=1;j<=Memory->S_cont[0].size();j++) {
				PropContourStream >> dummy >> Memory->S_cont[2][j-1] >> Memory->S_cont[3][j-1];
			} std::cout << "Propagator was loaded. And WILL NOT be recalculated" << std::endl;
		} else std::cout << "Cant open file!" << std::endl;
		PropContourStream.close();
	} else std::cout << "Propagator was NOT loaded. But WILL be recalculated." << std::endl;
}

// Export Propagator to file
// (exports all what is needed to perform Cauchy integration outside of this library: contour, weights, etc.)
//----------------------------------------------------------------------
void C_Quark::ExportPropagator(){
	int j;
	t_cmplx z_i,dz_i,temp;
	ofstream PropContourStream;
	PropContourStream.open ("Data_files/PropContourExportData_new.dat");
	std::cout << "Writing to file - \"Data_files/PropContourExportData_new.dat\"  " << std::endl;
	PropContourStream.precision(10);
	PropContourStream << 0.0001 << '\t' << 2000 << '\t' << Z2  << '\t' << std::endl;
	for (j=1;j<= Memory->S_cont[0].size() ;j++)
	{
		t_cmplx temp=1.0;
		z_i=(Memory->S_cont)[0][j-1];
		dz_i=temp*(Memory->S_cont)[1][j-1];
		PropContourStream << real(z_i) << '\t' << imag(z_i) << '\t'
				<< real(dz_i) << '\t' << imag(dz_i) << '\t' << real(temp) << '\t'
				<< real((Memory->S_cont)[2][j-1]) << '\t' << imag((Memory->S_cont)[2][j-1]) << '\t'
				<< real((Memory->S_cont)[3][j-1]) << '\t' << imag((Memory->S_cont)[3][j-1]) << std::endl;
	}
	PropContourStream.close();
	std::cout << "The quark propagator has beed exported" << std::endl;
}

// Saves Quark's A and B function on provided "Path" in provided "AmplutudeStorage"
//----------------------------------------------------------------------
void C_Quark::setPropagatorOnPath(std::vector<t_cmplxMatrix> (*AmplutudeStorage),t_cmplxArray1D (*Path)){
	std::cout << std::endl;
	std::cout << "Quark on Path extraction..." << std::endl;
	int num_points=(*Path).size();
	(*AmplutudeStorage).resize(num_points);
	t_cmplxArray1D temp_vector(5);
	t_cmplxMatrix temp_amp(2,1);
	std::cout << "points on the path - " << num_points << std::endl;
	for (int i = 0; i < num_points; i++)
	{
		temp_vector=getPropAt((*Path)[i]);
		temp_amp(0,0)=temp_vector[3];
		temp_amp(1,0)=temp_vector[4];
		(*AmplutudeStorage)[i]=temp_amp;
	}
	std::cout << "Quark on Path extraction finished." << std::endl;
}


// Gets sum A and B at 100 points. Used for Integration Test.
//----------------------------------------------------------------------
t_dArray1D C_Quark::GetTotalSum(){
	double Pu,Pd,x,scale,dp;
	t_cmplxArray1D storage;
	scale = 100;
	Pd=0.01;
	Pu=params.LimUk;
	dp=pow(10,(log10(Pu/Pd)/scale));
	t_dArray1D temp_result(2,0);
	x=Pd;
	for (int i = 0; i < scale; i++)
	{
		storage= getResultCauchyAt(x * x);
		temp_result[0]+=real(storage[0]);
		temp_result[1]+=real(storage[1]);
		x*=dp;
	}
	std::cout.precision(12);
	std::cout << temp_result[0] <<" "<< temp_result[1] << std::endl;
	return temp_result;
}


void C_Quark::CheckExtrapolation(){
	t_cmplx last_point;
		t_cmplxArray2D tempContour;
		t_cmplxArray2D tempABStorage;
		double M2=16.0;
		last_point = zz_rad[zz_rad.size()-1]*zz_rad[zz_rad.size()-1] - M2/4.0 + params.LimUk*params.LimDk;
		Geometry::C_ParabolaContour contour(ii*sqrt(M2)/2.0,
											ii*sqrt(M2),
											last_point);
		contour.setParabolaContour(zz_rad,w_rad,zz_line,w_line);
		tempContour = contour.getParabolaContour();

		for (int i = 0; i < tempContour.size(); i++) {
			tempABStorage.push_back(getResultCauchyAt(tempContour[0][i]));
		}

		std::cout << "Ready for Extrapolation?(enter)" << std::endl;
		std::cin.get();
		std::cout << "Yes Captain!" << std::endl;

		Memory->S_cont[0] = tempContour[0];
		Memory->S_cont[1] = tempContour[1];
		Memory->S_cont[2] = tempABStorage[0];
		Memory->S_cont[3] = tempABStorage[1];

		write_Prop_re(100);
}