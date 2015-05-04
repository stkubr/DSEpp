//
// Created by stkubr on 21.04.15.
//

#ifndef SHOW_PRECISION
#define SHOW_PRECISION 6
#endif

#ifndef DSEPP_BSE_TWOBODY_H
#define DSEPP_BSE_TWOBODY_H

#include <source/NumLibs/OneLoopIntegrator.hpp>
#include <source/NumLibs/Integrator.hpp>
#include <source/NumLibs/Integrator_Path.hpp>
#include <source/NumLibs/Support_functions.h>
#include <iomanip>
#include "BSE.h"
#include "BoundState_parameters.h"

class C_BSE_TwoBody: public C_BSE, public C_OneLoopIntegrator<t_cmplxMatrix, double, t_cmplxArray1D> {
protected:
    /// External Objects
    C_AbstractKernel * Kernel;
    C_Propagator * Parton_P;
    C_Propagator * Parton_M;

    C_DedicMem_BSE * Memory;

    C_BoundState_parameters params;

    int num_amplitudes;

    std::vector<std::vector<t_cmplxDirac>> threadloc_Projectors;

    t_cmplxArray2D threadloc_WeightCoeff;

    t_cmplxMatrix AMP,BUFFER_AMP,BUFFER_F_ex;

    t_dArray1D zz_rad, zz_cheb, zz_angleY, w_rad, w_cheb, w_angleY;

    bool flag_off_shell,flag_amp_desciption,flag_precalculation;

    // Constructor
    C_BSE_TwoBody() {
        Memory= static_cast<C_DedicMem_BSE *>(DedicMemFactory_BSE->CreateMemory());
        flag_amp_desciption=true;
        flag_precalculation=false;
        numIntegDimentions=3;
        params.setDefault();
    }

    virtual ~C_BSE_TwoBody(){
        delete Integrator_momentum;
        delete Integrator_angle_Z;
        delete Integrator_angle_Y;
        delete Memory;
    }

    void SetWeightCoeff(){
        int num_thread = omp_get_thread_num();
        for (int i = 0; i < num_amplitudes; i++) {
            threadloc_WeightCoeff[num_thread][i]=(threadloc_Projectors[num_thread][i]|threadloc_Projectors[num_thread][i]).Tr();
        }
    }

    void OrthogonalityCheck(){
        int num_thread = omp_get_thread_num();
        for (int i = 0; i < num_amplitudes; i++) {
            for (int j = 0; j < num_amplitudes; j++) {
                t_cmplx product;
                product = (threadloc_Projectors[num_thread][i]|threadloc_Projectors[num_thread][j]).Tr();
                std::cout << i+1 << "  " << j+1 << "  " << product << "  || ";
            }
            std::cout << std::endl << std::endl;
        }
    }

    void Initialization(){
        Integrator_momentum = C_Integrator_Line<t_cmplxMatrix, double>::createIntegrator(
                params.NumRadial, params.LimDk, params.LimUk, num_amplitudes, qgausleg_log_ID);
        Integrator_angle_Z=C_Integrator_Line<t_cmplxMatrix, double>::createIntegrator(
                params.NumCheb_nod1, params.LimDk, params.LimUk, num_amplitudes, qgauscheb_ID);
        Integrator_angle_Y=C_Integrator_Line<t_cmplxMatrix, double>::createIntegrator(
                params.NumAngleY, -1.0, 1.0, num_amplitudes, qgausleg_sym_ID);

        integrand_args.resize(3);
        ResizeALL();
        setInitialAMP();
    }

    void ResizeALL(){
        Integrator_momentum->getNodes(zz_rad,w_rad);
        Integrator_angle_Z->getNodes(zz_cheb,w_cheb);
        Integrator_angle_Y->getNodes(zz_angleY,w_angleY);

        BUFFER_AMP.Resize(params.NumRadial+1,params.Cheb_order*num_amplitudes+1);
        BUFFER_F_ex.Resize(params.NumRadial+1,(params.NumCheb_nod2)*num_amplitudes+1);
        AMP.Resize(params.NumRadial+1,params.Cheb_order*num_amplitudes+1);
    }

    void setInitialAMP(){
        for (int i = 1; i <= params.NumRadial; i++){
            for (int j = 1; j <=params.Cheb_order*num_amplitudes ; j++){
                AMP(i,j)=1.0;
            }
        }
    }

   C_Kinematics_1loop SetIntMomenta(t_cmplx x, t_cmplx y, t_cmplx z){
       C_Kinematics_1loop Momenta;
       Momenta.SetVectors_k(x,y,z);
       Momenta.P = C_BSE_TwoBody::threadloc_Momenta[omp_get_thread_num()].P;
       Momenta.p = C_BSE_TwoBody::threadloc_Momenta[omp_get_thread_num()].p;
       Momenta.SetVestors_k_for_S(params.zetta_part,Momenta.k);
       Momenta.SetVectors_q();
       return Momenta;
   }

    void PreCalculation(){
        Memory->resizeAmpStorage(num_amplitudes,params.NumRadial*params.NumCheb_nod1*params.NumAngleY);
        int int_ctr=0;
        for (int i = 1; i < zz_rad.size(); i++){
            for (int j = 1; j < zz_cheb.size(); j++){
                for (int k = 1; k < zz_angleY.size(); k++){
                    C_Kinematics_1loop momenta = SetIntMomenta(sqrt(zz_rad[i]),zz_angleY[k],zz_cheb[j]);
                    std::vector<t_cmplxDirac> WaveFunctions = SetWaveFunctions(momenta);
                    for (int amp_ctr = 0; amp_ctr < num_amplitudes; amp_ctr++){ Memory->setAmpStorage(amp_ctr, int_ctr, &WaveFunctions[amp_ctr]); }
                    int_ctr++;
                }
            }
        }
        flag_precalculation=true;
        //cin.get();
    }

    std::vector<t_cmplxDirac> SetWaveFunctions(C_Kinematics_1loop & Momenta){
        std::vector<t_cmplxDirac> WaveFunctions(num_amplitudes);
        if(!flag_precalculation){
            std::vector<t_cmplxDirac> Amplitudes(num_amplitudes);
            SetDiracStructures(Momenta.k,Momenta.P, Amplitudes);
            t_cmplxDirac S_p, S_m;
            setPropagators(Momenta.k_p, Momenta.k_m, S_p, S_m);
            for (int i = 0; i < num_amplitudes; i++){ WaveFunctions[i]=S_p*Amplitudes[i]*S_m;}
        } else {
            for (int i = 0; i < num_amplitudes; i++){ WaveFunctions[i]=Memory->AmpStorage[i][threadloc_Integ_ctr[omp_get_thread_num()]];}
        }
        return WaveFunctions;
    }

    t_cmplxDirac SetFullWaveFunction(std::vector<t_cmplxDirac> & WaveFunctions, t_cmplxArray1D & U_amp){
        t_cmplxDirac FullWaveFunction = WaveFunctions[0]; // start of the magic
        FullWaveFunction.Zero(); // end of the magic
        for (int i = 0; i < num_amplitudes; i++){ FullWaveFunction+=WaveFunctions[i]*U_amp[i]; }
        return FullWaveFunction;
    }

    virtual t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){
        t_cmplxMatrix dummy; std::cout << "Error - virtual call" << std::endl; assert(false); return dummy;};
    virtual void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> & DiracStructure){
        std::cout << "Error - virtual call" << std::endl; assert(false);};

    void preparePropagators(){
        Parton_P->dressPropagator();
        Parton_M->dressPropagator();
    }

    virtual void setPropagators(t_cmplxVector & K_plus, t_cmplxVector  & K_minus,
                        t_cmplxDirac & S_p, t_cmplxDirac & S_m)
    {
        t_cmplxArray1D quark_temp_sigma;
        quark_temp_sigma= Parton_P->PropagatorAtPoint(K_plus * K_plus);
        S_p=(-1.0*ii*(K_plus*Z)*quark_temp_sigma[3] + I*quark_temp_sigma[4]);

        quark_temp_sigma= Parton_M->PropagatorAtPoint(K_minus * K_minus);
        S_m=(-1.0*ii*(K_minus*Z)*(quark_temp_sigma[3]) + I*(quark_temp_sigma[4]));
    }

    t_cmplxArray1D SetDressing_normal(t_cmplx z){
        t_cmplxArray1D U_amp(num_amplitudes,0);
            for (int i = 1; i <= params.Cheb_order ; i++){
                t_cmplx temp_in;
                temp_in=Cheb_polinoms(z,(i-1)*2);
                for (int j = 0; j < num_amplitudes ; j++){
                    int momentum_inx = threadloc_momentum_inx[omp_get_thread_num()];
                    U_amp[j]+=temp_in*(AMP(momentum_inx,i+j*params.Cheb_order));
                }
            }
        return U_amp;
    }

    t_cmplxMatrix Integrand_vector(t_cmplxArray1D args){
        t_cmplxMatrix result(num_amplitudes,1), pre_result(num_amplitudes,1);
        t_cmplx Kinematic_factor;
        t_cmplx x,y,z;
        x=sqrt(args[0]);
        z=args[1];
        y=args[2];
        C_Kinematics_1loop momenta = SetIntMomenta(x,y,z);
        t_cmplxArray1D U_amp = SetDressing_normal(z);
        Kinematic_factor=-1.0/(8.0*pi*pi*pi*pi)*(args[0]);
        std::vector<t_cmplxDirac> WaveFunctions = SetWaveFunctions(momenta);
        t_cmplxDirac FullWaveFunction = SetFullWaveFunction(WaveFunctions, U_amp);
        pre_result = traceWithKernel_Vector(FullWaveFunction, momenta);
        result = Kinematic_factor*pre_result;
        return result;
    }

    t_cmplxMatrix traceWithKernel_Vector(t_cmplxDirac & FullWaveFunction,C_Kinematics_1loop & momenta){
        t_cmplxMatrix result(num_amplitudes,1), pre_result(num_amplitudes,1);
        t_cmplxVector k_p_P;
        k_p_P = (momenta.k + momenta.p - momenta.P)/2.0;
        bool flag_reset_kernel=true;
        for (int i = 0; i < num_amplitudes; i++) {
            pre_result(i,0)=Kernel->TraceKernelWithoutStoring(threadloc_Projectors[omp_get_thread_num()][i],
                                                          FullWaveFunction,
                                                          momenta.q, momenta.k, k_p_P,
                                                          flag_reset_kernel);
            flag_reset_kernel=false;
        }
        result = DisentangleAmps(&pre_result);
        return result;
    }

    t_cmplxMatrix CalcBSA(t_cmplx _p, t_cmplx _P, int proj, std::function<t_cmplxMatrix(t_cmplxArray1D)> bound_member_fn){
        t_cmplxMatrix BSAs(num_amplitudes,1);
        t_cmplxMatrix BUFFER_dataAmp_ex(num_amplitudes,params.NumCheb_nod2+1);
        threadloc_Momenta[omp_get_thread_num()].SetVector_P(_P);
        for (int j=1;j<=params.NumCheb_nod2;j++){
            t_cmplx z_v = zz_cheb[j];
            threadloc_Momenta[omp_get_thread_num()].SetVectors_p(z_v,_p);
            SetDiracStructures(threadloc_Momenta[omp_get_thread_num()].p,
                               threadloc_Momenta[omp_get_thread_num()].P,
                               threadloc_Projectors[omp_get_thread_num()]);
            SetWeightCoeff();
            BSAs=calcIntegral3D(&bound_member_fn, num_amplitudes, 1);
            for (int i = 0; i < num_amplitudes ; i++)
            {
                t_cmplx Bare_vertex=0.0;
                if (flag_off_shell && i==0) {Bare_vertex=2.0/pi* Parton_P->DressingFactor();}
                BUFFER_dataAmp_ex(i,j)=BSAs(i,0)+Bare_vertex;
            }
        }
        return BUFFER_dataAmp_ex;
    }

    t_cmplxMatrix DeProj(int proj, t_cmplxMatrix & BUFFER_dataAmp_ex){
        t_cmplxMatrix result(num_amplitudes,1);
        for (int i = 0; i < num_amplitudes ; i++) {
            result(i,0)=0.0;
            for (int j=1;j<=params.NumCheb_nod2;j++) {
                t_cmplx z_v=zz_cheb[j];
                result(i,0) += w_cheb[j]*BUFFER_dataAmp_ex(i,j)*Cheb_polinoms(z_v,(proj-1)*2);
            }
        }
        return result;
    }

    void setBufferIn(int i, t_cmplxMatrix & BUFFER_dataAmp_ex){
        for (int ampl = 0; ampl < num_amplitudes ; ampl++) {
            for (int w = 1; w <= params.NumCheb_nod2 ; w++) {
                BUFFER_F_ex(i,w + params.NumCheb_nod2*(ampl))=BUFFER_dataAmp_ex(ampl,w);
            }
        }
    }

    void setBufferOut(int i, t_cmplxMatrix & BUFFER_dataAmp_ex){
        for (int ampl = 0; ampl < num_amplitudes ; ampl++) {
            for (int w = 1; w <=params.NumCheb_nod2; w++) {
                BUFFER_dataAmp_ex(ampl,w)=BUFFER_F_ex(i,w + params.NumCheb_nod2*(ampl))*1.0;
            }
        }
    }

//
//------------------------------------------------------------------
    void BSA_step(t_cmplx P)
    {
        std::vector<t_cmplxMatrix> threadloc_BUFFER_dataAmp_ex(omp_get_max_threads());
        for (int proj_cheb = 1; proj_cheb <=params.Cheb_order ; proj_cheb++)
        {
#pragma omp parallel
 {//start of pragma
                std::function<t_cmplxMatrix(t_cmplxArray1D)> bound_member_fn =
                        std::bind(&C_BSE_TwoBody::Integrand_vector, this, std::placeholders::_1);
#pragma omp for
                for (int i = 1; i <= params.NumRadial; i++) {
                    t_cmplx p2=zz_rad[i];
                    BUFFER_AMP(i,0)=p2;
                    t_cmplxMatrix BUFFER_dataAmp_ex;
                    if(proj_cheb==1) {
                        threadloc_BUFFER_dataAmp_ex[omp_get_thread_num()] = CalcBSA(sqrt(p2),P,proj_cheb,bound_member_fn);
                        setBufferIn(i,threadloc_BUFFER_dataAmp_ex[omp_get_thread_num()]);
                    } else {
                        setBufferOut(i,threadloc_BUFFER_dataAmp_ex[omp_get_thread_num()]);
                    }
                    t_cmplxMatrix Temp_matrix=DeProj(proj_cheb,threadloc_BUFFER_dataAmp_ex[omp_get_thread_num()]);
                    for (int ampl = 0; ampl < num_amplitudes ; ampl++) {
                        BUFFER_AMP(i,proj_cheb+params.Cheb_order*(ampl))=Temp_matrix(ampl,0);
                    }
                }
            }//end of pragma
        }
    }

    t_cmplx DressBSA(t_cmplx P, int steps) {
        t_cmplx Lambda_EV;
        threadloc_Momenta[omp_get_thread_num()].SetVector_P(P);
        preparePropagators();
        PreCalculation();
        std::cout << "Dressing BSA... " << std::endl << std::endl;
        std::cout << setprecision(SHOW_PRECISION);
        setInitialAMP();
        for (int kk = 1; kk <= steps; kk++) {
            std::cout << "Step number - " << kk << std::endl;
            BSA_step(P);
            Lambda_EV=calcLambdaEV();
            PreNormBSA();
            PrintLine('-');
            std::cout  << "Eigenvalue.v1 - " << Lambda_EV << "  at P = " << "  " << P << std::endl;
            PrintLine('-');
            if(flag_amp_desciption) WriteAmplitudeDescyption();
            PrintLine('#');
        }
        flag_precalculation=false;
        return Lambda_EV;
    }

    t_cmplx calcLambdaEV(){
        t_cmplx ff1,ff2;
        ff1=0;ff2=0;
        for (int k = 1; k <=1; k++){
            for (int i = 1; i <= params.NumRadial; i++) {
                ff1+=(BUFFER_AMP(i,k)*AMP(i,k));
                ff2+=(AMP(i,k)*AMP(i,k));
            }
        }
        return ff1/ff2;
    }

    void PreNormBSA() {
        for (int i = 1; i <= params.NumRadial; i++) {
            for (int j = 1; j <= params.Cheb_order*(num_amplitudes); j++) {
                AMP(i,0)=BUFFER_AMP(i,0);
                if(flag_off_shell)AMP(i,j)=BUFFER_AMP(i,j);
                else AMP(i,j)=BUFFER_AMP(i,j)/(real(BUFFER_AMP(1,1)));
            }
        }
    }

    void WriteAmplitudeDescyption() {
        std::cout << fixed;
        std::cout << std::setprecision(SHOW_PRECISION);
        std::cout << "Dirac Amplitudes" << std::endl;
        PrintSpaces((SHOW_PRECISION+2)*2+3+2);
        for (int i = 0; i < num_amplitudes; i++)
        {
            std::cout << i;
            for (int j = 1; j <= params.Cheb_order; j++)
            {
                PrintSpaces((SHOW_PRECISION+2)*2+3);
            }
            std::cout << "";
        }
        std::cout << std::endl;
        std::cout << "Chebyshev Projections" << std::endl << "    p^2";
        PrintSpaces((SHOW_PRECISION+2)*2+3-5);
        for (int j = 0; j <num_amplitudes; j++)
        {
            for (int i = 1; i <=params.Cheb_order ; i++)
            {
                std::cout << (i-1)*2;
                PrintSpaces((SHOW_PRECISION+2)*2+3);
            }
        }
        std::cout << std::endl << AMP << std::endl;
    }

public:
    void dressBSE(t_cmplxVector P){
        DressBSA(P[3],10); // 4th component (rest frame)
    }

    double checkSum_PowerMethod(){
        flag_amp_desciption = false;
        DressBSA(t_cmplx(0,0.1),10);
        t_cmplx checkSum = 0;
        for (int i = 1; i <= params.NumRadial; i++){
            for (int j = 1; j <=params.Cheb_order*num_amplitudes ; j++){
                checkSum += AMP(i,j);
            }
        }
        return norm(checkSum);
    }

    virtual double checkSum_EVMatrixNorm(){  return 0;}

    void setBSAonPath(t_cmplxArray2D & AmplitudePath,t_cmplxArray1D & Path ,t_cmplx P) {
        std::cout << std::endl;
        std::cout << "On Path calculation..." << std::endl;
        flag_amp_desciption=true;
        DressBSA(P,10);
        AmplitudePath.resize(Path.size(), t_cmplxArray1D(num_amplitudes));
        std::cout << "Initialized" << std::endl;
        std::cout << "points on the path - " << Path.size() << std::endl;
#pragma omp parallel
        {//start of pragma
            std::function<t_cmplxMatrix(t_cmplxArray1D)> bound_member_fn =
                    std::bind(&C_BSE_TwoBody::Integrand_vector, this, std::placeholders::_1);
            t_cmplxMatrix Temp_matrix(num_amplitudes,1);
#pragma omp for
            for (int i = 0; i < Path.size(); i++)
            {
                Temp_matrix = CalcBSA(sqrt(Path[i]), P,1,bound_member_fn);
                for (int j = 0; j < num_amplitudes; ++j) {
                    AmplitudePath[i][j]=Temp_matrix(j,0);
                }
            }
        }//end of pragma
        std::cout << "On Path calculation finished." << std::endl;
    }

    void linkToKernel(C_AbstractKernel * _Kernel){
        Kernel=_Kernel;
    }

    void linkToPartons(std::vector<C_Propagator*> _partons){
        Parton_P=_partons[0];
        Parton_M=_partons[1];
    }
};


#endif //DSEPP_BSE_TWOBODY_H