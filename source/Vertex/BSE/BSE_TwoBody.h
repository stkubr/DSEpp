//
// Created by stkubr on 21.04.15.
//

#define NUM_PRECISION 6

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
public:
    // Partons Dirac structure
    t_cmplxDirac S_p, S_m, Y_T;

    /// External Objects
    C_AbstractKernel * Kernel;
    C_Propagator * Parton_P;
    C_Propagator * Parton_M;

    C_BoundState_parameters params;

    std::vector<t_cmplxDirac> Amplitudes,Projectors,WaveFunctions;
    t_cmplxDirac FullWaveFunction,FullAmplitude;
    t_cmplxArray1D U_amp, WeightCoeff;

    t_dArray1D zz_rad, zz_cheb, zz_angleY , zz_cauchy, w_rad, w_cheb, w_angleY, w_cauchy, proj_amp;

    int num_amplitudes,index_p;
    bool flag_off_shell,flag_amp_desciption,flag_precalculation;
    double NormalizationConst;

    void (C_BSE_TwoBody::*SetDressing_ref) (t_cmplx);
    C_Integrator_Path<t_cmplxArray1D,t_cmplxArray3D,t_cmplx> * Integ_cauchy_long;
    t_cmplxMatrix (C_BSE_TwoBody::*integrand) (t_cmplxArray1D);
    t_cmplxMatrix (C_BSE_TwoBody::*GetBSA_ref) ();

    t_cmplxMatrix Zero,dataAmp,AMP,BUFFER_AMP,BUFFER_F_ex,BUFFER_dataAmp_ex;

    // Constructor
    C_BSE_TwoBody() {
        Memory= static_cast<C_DedicMem_BSE *>(DedicMemFactory_BSE->CreateMemory());
        flag_amp_desciption=true;
        flag_precalculation=false;
        NormalizationConst=1.0;
        numIntegDimentions=3;

        params.setDefault();
        params.print();
    }

    void SetWeightCoeff(){
        for (int i = 0; i < num_amplitudes; i++) { WeightCoeff[i]=(Projectors[i]|Projectors[i]).Tr(); }
    }

    void OrthogonalityCheck(){
        for (int i = 0; i < num_amplitudes; i++) {
            for (int j = 0; j < num_amplitudes; j++) {
                t_cmplx product;
                product = (Projectors[i]|Projectors[j]).Tr();
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
        //Integ_cauchy_long=C_Integrator_Path<t_cmplxArray1D,t_cmplxArray3D,t_cmplx>::createIntegrator(params.NumRadial_Contour, sqrt(params.LimDk), sqrt(params.LimUk), num_amplitudes, qcauchyleg_lin_ID);

        integrand_args.resize(3);
        ResizeALL();
        setInitialAMP();
    }

    void ResizeALL(){
        Integrator_momentum->getNodes(zz_rad,w_rad);
        Integrator_angle_Z->getNodes(zz_cheb,w_cheb);
        Integrator_angle_Y->getNodes(zz_angleY,w_angleY);
        //Integ_cauchy_long->getNodes(&zz_cauchy,&w_cauchy);

        Zero.Resize(num_amplitudes,1);
        dataAmp.Resize(num_amplitudes,1);
        BUFFER_AMP.Resize(params.NumRadial+1,params.Cheb_order*num_amplitudes+1);
        BUFFER_F_ex.Resize(params.NumRadial+1,(params.NumCheb_nod2)*num_amplitudes+1);
        BUFFER_dataAmp_ex.Resize(num_amplitudes,params.NumCheb_nod2+1);
        AMP.Resize(params.NumRadial+1,params.Cheb_order*num_amplitudes+1);
        U_amp.resize(num_amplitudes);
        proj_amp.resize(params.Cheb_order+1);
    }

    void setInitialAMP(){
        for (int i = 1; i <= params.NumRadial; i++){
            for (int j = 1; j <=params.Cheb_order*num_amplitudes ; j++){
                AMP(i,j)=1.0;
            }
        }
        std::cout << AMP << std::endl;
    }

    void SetIntMomenta(t_cmplx x, t_cmplx y, t_cmplx z){
        Momenta.SetVectors_k(x,y,z);
        Momenta.SetVectors_q();
        Momenta.SetVestors_k_for_S(params.zetta_part,Momenta.k);
    }

    virtual C_BSE_TwoBody * MakeCopy()
    {
        return new C_BSE_TwoBody(*this);
    }

 /*   t_cmplxMatrix IntegAngleY(){
        return Integ_angle_Y->getResult(&C_BSE_TwoBody::f3);
    }*/

    void PreCalculation(){
        Memory->resizeAmpStorage(num_amplitudes,params.NumRadial*params.NumCheb_nod1*params.NumAngleY);
        int int_ctr=0;
        for (int i = 1; i < zz_rad.size(); i++){
            for (int j = 1; j < zz_cheb.size(); j++){
                for (int k = 1; k < zz_angleY.size(); k++){
                    SetIntMomenta(sqrt(zz_rad[i]),zz_angleY[k],zz_cheb[j]);
                    SetWaveFunctions();
                    for (int amp_ctr = 0; amp_ctr < num_amplitudes; amp_ctr++){ Memory->setAmpStorage(amp_ctr, int_ctr, &WaveFunctions[amp_ctr]); }
                    int_ctr++;
                }
            }
        }
        flag_precalculation=true;
        //cin.get();
    }

 /*   void CalcEVMatrix(Eigen::ComplexEigenSolver<Eigen::MatrixXcf> * ces){
        Memory->ResizeEVMatrix(params.NumRadial,params.NumCheb_nod1,num_amplitudes,1);
        //std::cout << "Memory->EVMatrix call" << std::endl;
#pragma omp parallel
        {//start of pragma
            C_BSE_Hadron_Base * bsa_copy_omp;
            bsa_copy_omp=MakeCopy();
            t_cmplxMatrix Temp_Matrix;
            double _p2,_k2,_z,_y,_zp;
            double _w_zp,_w_k2,_w_z;
            bsa_copy_omp->k_col=0;
            bsa_copy_omp->Int_counter=0;
#pragma omp for
            for (int p_ctr = 1; p_ctr < zz_rad.size() ; p_ctr++){
                bsa_copy_omp->index_zp=0;
                bsa_copy_omp->index_p=p_ctr;
                _p2=zz_rad[p_ctr];
                for (int zp_ctr = 1; zp_ctr < zz_cheb.size() ; zp_ctr++){
                    bsa_copy_omp->index_zp=zp_ctr;
                    _zp=zz_cheb[zp_ctr];
                    bsa_copy_omp->Momenta.SetVectors_p(_zp,sqrt(_p2));
                    bsa_copy_omp->SetDiracStructures(bsa_copy_omp->Momenta.p,bsa_copy_omp->Momenta.P,&bsa_copy_omp->Projectors);
                    bsa_copy_omp->SetWeightCoeff();
                    //FullWaveFunction=Projectors[0];
                    //_w_zp=w10_ch_ex[zp_ctr];
                    for (int k_ctr = 1; k_ctr < zz_rad.size(); k_ctr++){
                        //k_ctr = 1;
                        _k2=zz_rad[k_ctr];
                        _w_k2=w_rad[k_ctr];
                        bsa_copy_omp->integrand_args[0]=_k2;
                        bsa_copy_omp->flag_sigma=false;
                        for (int z_ctr = 1; z_ctr < zz_cheb.size() ; z_ctr++){
                            _z=zz_cheb[z_ctr];
                            _w_z=w_cheb[z_ctr];
                            bsa_copy_omp->integrand_args[1]=_z;
                            //Temp_Matrix=bsa_copy_omp->SetMatrixAtPoint(_P,_p2,_zp,_k2,_z);
                            Temp_Matrix=bsa_copy_omp->IntegAngleY();
                            for (int p_amp_ctr = 0; p_amp_ctr < num_amplitudes ; p_amp_ctr++){
                                for (int k_amp_ctr = 0; k_amp_ctr < num_amplitudes ; k_amp_ctr++){

                                    //std::cout << k_ctr << "  " << z_ctr << "  " << z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr << std::endl;
                                    //std::cout << k_ctr << "  " << z_ctr << "  " << zp_ctr + (p_ctr-1)*(NumRadialPoints-z_ctr*(NumRadialPoints-1))-1+ NumRadialPoints*(num_cheb_nod1)*p_amp_ctr << "  " << z_ctr + (k_ctr-1)*(NumRadialPoints-z_ctr*(NumRadialPoints-1)) + NumRadialPoints*(num_cheb_nod1)*k_amp_ctr-1 << std::endl;
                                    Memory->EVMatrix(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr, z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1) + params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr)=pi/2.0*_w_k2*_w_z*Temp_Matrix(p_amp_ctr,k_amp_ctr);
                                }
                            }
                            //std::cout << Temp_Matrix(0,0) << std::endl;
                            //cin.get();
                        }
                        bsa_copy_omp->k_col++;
                    }
                    bsa_copy_omp->k_col=0;
                    bsa_copy_omp->Int_counter=0;
                }
            }
            delete bsa_copy_omp;
        }// end of pragma
        //std::cout << "EVMatrix is full. EigenValues computation engaged..." << std::endl;
        ces->compute(Memory->EVMatrix);
    }

    t_cmplxArray2D SetEVMatrix(t_cmplx _P){
        Momenta.SetVector_P(_P);
        PreCalculation();
        t_dArray1D zz_rad_temp,zz_cheb_temp;
        SetDressing_ref=&C_BSE_Hadron_Base::SetDressing_normal;
        GetBSA_ref=&C_BSE_Hadron_Base::GetBSA_matrix;

        setInitialAMP();

        Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;

        Memory->ResizeEVMatrix(params.NumRadial,params.NumCheb_nod1,num_amplitudes,1);
        //std::cout << "Memory->EVMatrix call" << std::endl;
#pragma omp parallel
        {//start of pragma
            C_BSE_Hadron_Base * bsa_copy_omp;
            bsa_copy_omp=MakeCopy();
            t_cmplxMatrix Temp_Matrix;
            double _p2,_k2,_z,_y,_zp;
            double _w_zp,_w_k2,_w_z;
            bsa_copy_omp->k_col=0;
            bsa_copy_omp->Int_counter=0;
#pragma omp for
            for (int p_ctr = 1; p_ctr < zz_rad.size() ; p_ctr++){
                bsa_copy_omp->index_zp=0;
                bsa_copy_omp->index_p=p_ctr;
                _p2=zz_rad[p_ctr];
                for (int zp_ctr = 1; zp_ctr < zz_cheb.size() ; zp_ctr++){
                    bsa_copy_omp->index_zp=zp_ctr;
                    _zp=zz_cheb[zp_ctr];
                    bsa_copy_omp->Momenta.SetVectors_p(_zp,sqrt(_p2));
                    bsa_copy_omp->SetDiracStructures(bsa_copy_omp->Momenta.p,bsa_copy_omp->Momenta.P,&bsa_copy_omp->Projectors);
                    bsa_copy_omp->SetWeightCoeff();
                    //FullWaveFunction=Projectors[0];
                    //_w_zp=w10_ch_ex[zp_ctr];
                    for (int k_ctr = 1; k_ctr < zz_rad.size(); k_ctr++){
                        //k_ctr = 1;
                        _k2=zz_rad[k_ctr];
                        _w_k2=w_rad[k_ctr];
                        bsa_copy_omp->integrand_args[0]=_k2;
                        bsa_copy_omp->flag_sigma=false;
                        for (int z_ctr = 1; z_ctr < zz_cheb.size() ; z_ctr++){
                            _z=zz_cheb[z_ctr];
                            _w_z=w_cheb[z_ctr];
                            bsa_copy_omp->integrand_args[1]=_z;
                            //Temp_Matrix=bsa_copy_omp->SetMatrixAtPoint(_P,_p2,_zp,_k2,_z);
                            Temp_Matrix=bsa_copy_omp->IntegAngleY();
                            for (int p_amp_ctr = 0; p_amp_ctr < num_amplitudes ; p_amp_ctr++){
                                for (int k_amp_ctr = 0; k_amp_ctr < num_amplitudes ; k_amp_ctr++){

                                    //std::cout << k_ctr << "  " << z_ctr << "  " << z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr << std::endl;
                                    //std::cout << k_ctr << "  " << z_ctr << "  " << zp_ctr + (p_ctr-1)*(NumRadialPoints-z_ctr*(NumRadialPoints-1))-1+ NumRadialPoints*(num_cheb_nod1)*p_amp_ctr << "  " << z_ctr + (k_ctr-1)*(NumRadialPoints-z_ctr*(NumRadialPoints-1)) + NumRadialPoints*(num_cheb_nod1)*k_amp_ctr-1 << std::endl;
                                    Memory->EVMatrix(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr, z_ctr-1 + (k_ctr-1)*(params.NumCheb_nod1) + params.NumRadial*(params.NumCheb_nod1)*k_amp_ctr)=pi/2.0*_w_k2*_w_z*Temp_Matrix(p_amp_ctr,k_amp_ctr);
                                }
                            }
                            //std::cout << Temp_Matrix(0,0) << std::endl;
                            //cin.get();
                        }
                        bsa_copy_omp->k_col++;
                    }
                    bsa_copy_omp->k_col=0;
                    bsa_copy_omp->Int_counter=0;
                }
            }
            delete bsa_copy_omp;
        }// end of pragma
        //std::cout << "EVMatrix is full. EigenValues computation engaged..." << std::endl;
        ces.compute(Memory->EVMatrix);
        flag_precalculation=false;
        //CalcEVMatrix(&ces);
        Eigen::VectorXcf EV=ces.eigenvalues();


        auto Parity = [&](int num_state) -> t_cmplx {
            t_cmplx parity=0.0;
            for (int p_ctr = 1; p_ctr < 2 ; p_ctr++){
                for (int p_amp_ctr = 0; p_amp_ctr < 1 ; p_amp_ctr++){
                    for (int zp_ctr = 1; zp_ctr < zz_cheb.size() ; zp_ctr++){
                        t_cmplx summ,diff;

                        summ = ces.eigenvectors().col(num_state)(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr) +
                               ces.eigenvectors().col(num_state)(zz_cheb.size() - zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr);

                        diff = ces.eigenvectors().col(num_state)(zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr) -
                               ces.eigenvectors().col(num_state)(zz_cheb.size() - zp_ctr-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*p_amp_ctr);

                        if(fabs(real(summ))<=fabs(real(diff))) {parity = -1.0;}
                        else{ parity = 1.0; }
                    }
                }
            }
            return parity;
        };

        //std::cout << "EigenValues computation is done. The eigenvalues of EVMatrix are obtained." << std::endl;
        int i= EV.size()-1;
        t_cmplxArray2D Dominant_EV_and_parity(2);
        while( i > EV.size()-10)
        {
            //std::cout << EV[i] << "  " << i << std::endl;
            Dominant_EV_and_parity[0].push_back(EV[i]);
            Dominant_EV_and_parity[1].push_back(Parity(i));
            //std::cout << i << "  " << EV[i] << "  " << Parity(i)<< std::endl;
            i--;
        }
        //std::cout << ces.eigenvectors().col(EV.size()-1) << std::endl;
        return Dominant_EV_and_parity;
    }

    void DrawBSA_matrix(t_cmplx _P, int _state, int amp_num){
        int num_state;
        Momenta.SetVector_P(_P);
        PreCalculation();
        t_dArray1D zz_rad_temp,zz_cheb_temp;
        SetDressing_ref=&C_BSE_Hadron_Base::SetDressing_normal;
        GetBSA_ref=&C_BSE_Hadron_Base::GetBSA_matrix;

        flag_precalculation=false;
        Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
        //ces.compute(Memory->EVMatrix);
        CalcEVMatrix(&ces);
        Eigen::VectorXcf EV=ces.eigenvalues();
        num_state = EV.size()- _state;

        ofstream DrawBSA_matrix;
        DrawBSA_matrix.open ("Data_files/BSEs_matrix.dat");

        int i=EV.size()-1;

        while(i > EV.size()-6)
        {
            //std::cout << EV[i] << "  " << i << std::endl;
            i--;
        }

        t_cmplxArray1D result(zz_rad.size());
        for (int p_ctr = 1; p_ctr <= params.NumRadial ; p_ctr++){
            result[p_ctr-1]=0;
            for (int j=1;j<=params.NumCheb_nod1;j++)
            {
                t_cmplx z_v,element;
                element=ces.eigenvectors().col(num_state)(j-1 + (p_ctr-1)*(params.NumCheb_nod1)+ params.NumRadial*(params.NumCheb_nod1)*amp_num);
                z_v=zz_cheb[j];
                result[p_ctr-1] += w_cheb[j]*Cheb_polinoms(z_v,0)*real(element);
                //std::cout << w_cheb[j] << "  " << U_ex(z_v) << "  " << real(element) << std::endl;
            }
            DrawBSA_matrix << zz_rad[p_ctr] << "  " << real(result[p_ctr-1]) << std::endl;
            //std::cout << zz_rad[p_ctr] << "  " << result[p_ctr-1] << std::endl;
            //cin.get();
        }
        DrawBSA_matrix.close();
    }
*/

    void SetWaveFunctions(){
     if(!flag_precalculation){
         SetDiracStructures(Momenta.k,Momenta.P,&Amplitudes);
         setPropagators(&Momenta.k_p, &Momenta.k_m);
         for (int i = 0; i < num_amplitudes; i++){ WaveFunctions[i]=S_p*Amplitudes[i]*S_m;}
     }
     else{
         for (int i = 0; i < num_amplitudes; i++){ WaveFunctions[i]=Memory->AmpStorage[i][threadloc_Integ_ctr[omp_get_thread_num()]];}
     }
 }

    void SetFullWaveFunction(){
        FullWaveFunction.Zero();
        for (int i = 0; i < num_amplitudes; i++){ FullWaveFunction+=WaveFunctions[i]*U_amp[i]; }
    }

    t_cmplxMatrix GetBSA(){
        t_cmplxMatrix result,pre_result;
        t_cmplxVector k_p_P;
        k_p_P=(Momenta.k + Momenta.p - Momenta.P)/2.0;
        result.Resize(num_amplitudes,1);
        pre_result.Resize(num_amplitudes,1);
        SetFullWaveFunction();
        bool flag_reset_kernel=true;
        for (int i = 0; i < num_amplitudes; i++) {
            pre_result(i,0)=Kernel->TraceKernelWithoutStoring((Projectors[i]),
                                                              FullWaveFunction,
                                                              Momenta.q,Momenta.k,
                                                              k_p_P,flag_reset_kernel);
            flag_reset_kernel=false;
        }
        result=DisentangleAmps(&pre_result);
        return result;
    }

    virtual t_cmplxMatrix DisentangleAmps(t_cmplxMatrix * pre_result){t_cmplxMatrix dummy; std::cout << "Error - virtual call" << std::endl; assert(false); return dummy;};
    virtual void SetDiracStructures(t_cmplxVector _k, t_cmplxVector _P, std::vector<t_cmplxDirac> * DiracStructure){std::cout << "Error - virtual call" << std::endl; assert(false);};
    virtual void setPropagators(t_cmplxVector *K_plus, t_cmplxVector *K_minus){std::cout << "Error - virtual call" << std::endl;}
    //virtual t_cmplxMatrix GetBSA(){t_cmplxMatrix dummy; std::cout << "Error - virtual call" << std::endl; assert(false); return dummy;};

    void SetDressing_normal(t_cmplx z){
        /*if (flag_sigma == true){
            flag_sigma = false;*/
            for (int j = 0; j < num_amplitudes ; j++){
                U_amp[j]=0.0;
            }
            for (int i = 1; i <= params.Cheb_order ; i++){
                t_cmplx temp_in;
                temp_in=Cheb_polinoms(z,(i-1)*2);
                for (int j = 0; j < num_amplitudes ; j++){
                    int momentum_inx = threadloc_momentum_inx[omp_get_thread_num()];
                    U_amp[j]+=temp_in*(AMP(momentum_inx,i+j*params.Cheb_order));
                }
            }
        //}
    }
/*
    void SetDressing_shifted(t_cmplx z){
        Momenta.ShiftMomenta(params.zetta_part);
        for (int j = 0; j < num_amplitudes ; j++){
            U_amp[j]=0.0;
        }
        for (int i = 1; i <= params.Cheb_order ; i++){
            t_cmplx temp_in;
            temp_in=Cheb_polinoms(z,(i-1)*2);
            for (int j = 0; j < num_amplitudes ; j++){
                U_amp[j]+=temp_in*(Memory->CauchyGrid[j+1][index_p][Complex_Int_counter]);
            }
        }
    }
*/
    t_cmplxMatrix Integrand(t_cmplxArray1D args){
        t_cmplxMatrix result;
        t_cmplx Kinematic_factor;
        t_cmplx x,y,z;
        x=sqrt(args[0]);
        z=args[1];
        y=args[2];
        SetIntMomenta(x, y, z);

        (this->*SetDressing_ref)(z);

        Kinematic_factor=-1.0/(8.0*pi*pi*pi*pi)*(args[0]);//(1.0 + args[0]/1000.0);//(1.0 + real(Momenta.P2)/1000.0)/(1.0 + real(Momenta.p2)/1000.0);
        SetWaveFunctions();
        result=Kinematic_factor*GetBSA();


        //Complex_Int_counter++;
        return result;
    }

    void CalcBSA(t_cmplx _p, t_cmplx _P, int proj, std::function<t_cmplxMatrix(t_cmplxArray1D)> bound_member_fn){
        t_cmplxMatrix result(num_amplitudes,1);
        Momenta.SetVector_P(_P);
        for (int j=1;j<=params.NumCheb_nod2;j++){
            t_cmplx z_v = zz_cheb[j];
            Momenta.SetVectors_p(z_v,_p);

            SetDiracStructures(Momenta.p,Momenta.P,&Projectors);
            SetWeightCoeff();
            FullWaveFunction=Projectors[0];

            result=calcIntegral3D(&bound_member_fn, num_amplitudes, 1);
            for (int i = 0; i < num_amplitudes ; i++)
            {
                t_cmplx Bare_vertex=0.0;
                if (flag_off_shell && i==0) {Bare_vertex=2.0/pi* Parton_P->DressingFactor();}
                BUFFER_dataAmp_ex(i,j)=result(i,0)+Bare_vertex;
            }
        }
    }

    t_cmplxMatrix DeProj(int proj){
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

    void setBufferIn(int i){
        for (int ampl = 0; ampl < num_amplitudes ; ampl++) {
            for (int w = 1; w <= params.NumCheb_nod2 ; w++) {
                BUFFER_F_ex(i,w + params.NumCheb_nod2*(ampl))=BUFFER_dataAmp_ex(ampl,w);
            }
        }
    }

    void setBufferOut(int i){
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
        for (int proj_cheb = 1; proj_cheb <=params.Cheb_order ; proj_cheb++)
        {
//#pragma omp parallel
 {//start of pragma
                t_cmplx p2;
                std::function<t_cmplxMatrix(t_cmplxArray1D)> bound_member_fn =
                        std::bind(&C_BSE_TwoBody::Integrand, this, std::placeholders::_1);
                t_cmplxMatrix Temp_matrix(num_amplitudes,1);
//#pragma omp for
                for (int i = 1; i <= params.NumRadial; i++) {
                    index_p=i;
                    p2=zz_rad[i];
                    BUFFER_AMP(i,0)=p2;
                    if(proj_cheb==1) {
                        CalcBSA(sqrt(p2),P,proj_cheb,bound_member_fn);
                        setBufferIn(i);
                    } else {
                        setBufferOut(i);
                    }

                    Temp_matrix=DeProj(proj_cheb);

                    for (int ampl = 0; ampl < num_amplitudes ; ampl++) {
                        BUFFER_AMP(i,proj_cheb+params.Cheb_order*(ampl))=Temp_matrix(ampl,0);
                    }
                }
            }//end of pragma
        }
    }


    void dressBSE(t_cmplxVector P){
        DressBSA(t_cmplx(0,0.1),10);
    }

//
//------------------------------------------------------------------
    t_cmplx DressBSA(t_cmplx P, int steps)
    {
        t_cmplx ff1,ff2,ff3,Lambda_EV;
        ff1=0;ff2=0;ff3=0;
        Momenta.SetVector_P(P);
        PreCalculation();
        SetDressing_ref=&C_BSE_TwoBody::SetDressing_normal;
        GetBSA_ref=&C_BSE_TwoBody::GetBSA;
        std::cout << "Dressing BSA... " << std::endl << std::endl;
        std::cout << fixed;
        //std::cout << setprecision (NUM_PRECISION);

        //if(flag_load_AMP)
        setInitialAMP();

        for (int kk = 1; kk <= steps; kk++)
        {
            std::cout << "Step number - " << kk << std::endl;

            BSA_step(P);

            //std::cout << "Time now - " << (Get_Time() - Init_time)/NumThreadsUse << std::endl;

            for (int k = 1; k <=1; k++)
            {
                for (int i = 1; i <= params.NumRadial; i++)
                {
                    ff1+=(BUFFER_AMP(i,k)*AMP(i,k));
                    ff2+=(AMP(i,k)*AMP(i,k));
                }
            }

            Lambda_EV=(ff1/ff2);
            PrintLine('-');
            std::cout  << "Eigenvalue.v1 - " << ff1/ff2 << "  at P = " << "  " << P << std::endl;
            PrintLine('-');
            ff1=0;ff2=0;ff3=0;

            PreNormBSA();
            if(flag_amp_desciption) WriteAmplitudeDescyption();

            PrintLine('#');
        }
        //std::cout << "Time now - " << (Get_Time() - Init_time) << std::endl;
        //DrawBSA(P);

        //writeBSA();

        //flag_precalculation=false;
        return Lambda_EV;
    }

    void PreNormBSA()
    {
        for (int i = 1; i <= params.NumRadial; i++)
        {
            for (int j = 1; j <= params.Cheb_order*(num_amplitudes); j++)
            {
                AMP(i,0)=BUFFER_AMP(i,0);
                if(flag_off_shell)AMP(i,j)=BUFFER_AMP(i,j);//BUFFER_AMP(1,1);
                else AMP(i,j)=BUFFER_AMP(i,j)/(real(BUFFER_AMP(1,1)));

            }
        }
    }

    void WriteAmplitudeDescyption()
    {
        std::cout << fixed;
        std::cout << std::setprecision(NUM_PRECISION);
        std::cout << "Dirac Amplitudes" << std::endl;
        PrintSpaces((NUM_PRECISION+2)*2+3+2);
        for (int i = 0; i < num_amplitudes; i++)
        {
            std::cout << i;
            for (int j = 1; j <= params.Cheb_order; j++)
            {
                PrintSpaces((NUM_PRECISION+2)*2+3);
            }
            std::cout << "";
        }
        std::cout << std::endl;
        std::cout << "Chebyshev Projections" << std::endl << "    p^2";
        PrintSpaces((NUM_PRECISION+2)*2+3-5);
        for (int j = 0; j <num_amplitudes; j++)
        {
            for (int i = 1; i <=params.Cheb_order ; i++)
            {
                std::cout << (i-1)*2;
                PrintSpaces((NUM_PRECISION+2)*2+3);
            }
        }
        std::cout << std::endl << AMP << std::endl;
    }

 /*   void DrawBSA(t_cmplx _P){
        t_cmplxMatrix TempArray(num_amplitudes,1);
        ofstream temp_continuation;
        temp_continuation.open ("Data_files/BSEs.dat");
        for (int i = 1; i < zz_rad.size(); i++)
        {
            TempArray=CalcBSA(sqrt(zz_rad[i]),_P,1);
            temp_continuation << (zz_rad[i]);
            for (int amp = 0; amp < num_amplitudes; amp++) temp_continuation  << '	' <<  real(TempArray(amp,0));
            temp_continuation << std::endl;
        }
        temp_continuation.close();
    }
*/
/*
    double NormalizeBSA(t_cmplx _K){
        std::cout << std::endl;
        PrintLine('=');
        std::cout << std::endl;
        std::cout << "Calculation Norm Factor..." << std::endl;
        std::cout << std::endl;
        double BSE_Norm_factor;
        t_cmplx h=_K/5.0;
        t_cmplx res_int,N,Lambda_EV;
        flag_amp_desciption=false;
        std::vector <t_cmplx> deriv_step_norm(4);
        for (int i = 1; i <= 2; i++)
        {
            double i_count;
            i_count=i;
            std::cout << "Derivative at point " << _K+(i_count)*h << std::endl;
            deriv_step_norm[(i-1)]=DressBSA(_K+(i_count)*h,10);

            std::cout << "Derivative at point " << _K-(i_count)*h << std::endl;
            deriv_step_norm[(i-1)+2]=DressBSA(_K-(i_count)*h,10);
        }
        //flag_normalization=true;

        Lambda_EV=DressBSA(_K,10);

        SetDressing_ref=&C_BSE_Hadron_Base::SetDressing_normal;
        GetBSA_ref=&C_BSE_Hadron_Base::GetBSA_norm;

        res_int=quad3d()(0,0);
        //flag_normalization=false;
        N=(-deriv_step_norm[1] + 8.0*deriv_step_norm[0] - 8.0*deriv_step_norm[2] + deriv_step_norm[3])/(12.0*imag(h));
        std::cout << deriv_step_norm[0] << "  " << deriv_step_norm[1] << "  " << deriv_step_norm[2] << "  " << deriv_step_norm[3] << std::endl;
        std::cout << N << "  " << res_int << std::endl;
        //NORM=sqrt(1.0/fabs(real(N))/3.0/fabs(real(res_int))/Lambda_EV)/sqrt(2.0*imag(_K))*0.5;
        //NORM=sqrt(fabs(real(N))/fabs(real(res_int)))/sqrt(2.0*imag(K_v));

        BSE_Norm_factor=real(sqrt(Lambda_EV/fabs(real(N))/fabs(real(res_int)))*sqrt(2.0*imag(_K)));

        std::cout << "NormFactor = "<< BSE_Norm_factor << std::endl;

        return BSE_Norm_factor;
    }
*/
    void SetBSAonPath(t_cmplxArray1D (*AmplitudePath),t_cmplxArray1D (*Path) ,t_cmplx Q)
    {
        std::cout << std::endl;
        std::cout << "On Path calculation..." << std::endl;
        flag_amp_desciption=true;
        DressBSA(Q,1);
        int num_points;
        double Time;
        num_points=(*Path).size();
        (*AmplitudePath).resize(num_points);
        t_cmplxMatrix Temp(num_amplitudes,1);
        std::cout << "Initialized" << std::endl;
        std::cout << "points on the path - " << num_points << std::endl;
        Time=Get_Time();
#pragma omp parallel
        {//start of pragma
            C_BSE_TwoBody * bsa_copy_omp;
            bsa_copy_omp=MakeCopy();
            t_cmplxMatrix Temp_matrix(num_amplitudes,1);
#pragma omp for
            for (int i = 0; i < num_points; i++)
            {
               // (*AmplitudePath)[i]=NORM*bsa_copy_omp->CalcBSA(sqrt((*Path)[i]),Q,1)(0,0);
            }
            delete bsa_copy_omp;
        }//end of pragma
        std::cout << " Time for continuation spent - " << (Get_Time() - Time)/8.0 << std::endl;
        std::cout << "On Path calculation finished." << std::endl;
    }

public:
    C_DedicMem_BSE * Memory;

    void linkToKernel(C_AbstractKernel * _Kernel){
        Kernel=_Kernel;
    }

    void linkToPartons(std::vector<C_Propagator*> _partons){
        Parton_P=_partons[0];
        Parton_M=_partons[1];
    }
};


#endif //DSEPP_BSE_TWOBODY_H
