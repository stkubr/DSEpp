#pragma once

template <typename T_out, typename F, typename T_in> class C_MultiDimInt{
	public:
	double y_cx;
	T_out (F::*cxfunc) (T_in);
	T_out (C_MultiDimInt::*MultiInt) (T_in);
	vector<C_Integrator_Line<T_out,C_MultiDimInt,double> *>  IntArray;
	F * obj_ref;
	int D,current_dim;
	T_in temp_vector;
	
	//vector<C_Integrator_Line<T_out,F,double> *>  __IntArray
	C_MultiDimInt(int _D){
		IntArray.resize(2);
		IntArray[0]=C_Integrator_Line<t_cmplxMatrix,C_MultiDimInt,double>::createIntegrator(40, 0.0001, 1.0, 1, qgausleg_lin_ID);
		IntArray[1]=C_Integrator_Line<t_cmplxMatrix,C_MultiDimInt,double>::createIntegrator(40, 0.0001, 1.0, 1, qgausleg_lin_ID);
		temp_vector.resize(_D);
	}
	
	void setIntegrators(vector<C_Integrator_Line<T_out,F,double> *>  __IntArray, int _D){
		IntArray=__IntArray;
		temp_vector.resize(_D);
	}	
	
	T_out getResult(T_out (F::*func2) (T_in), F * obj )
	{
		cout << "GetResult in...1" << endl;
		cxfunc=func2;
		cout << "GetResult in...2" << endl;
		obj_ref=obj;
		cout << "GetResult in...3" << endl;
		T_out temp(1,0);
		//cout << cxfunc << "  " << func2 << "  " << obj_ref <<endl;
		
		temp = Multi_INT_cx(cxfunc);

		return temp;
	}
	
	T_out Multi_INT_cx(T_out (F::*func2) (T_in))
	{
		return IntArray[0]->getResult(&C_MultiDimInt::f1,this);
	}
	T_out f1 (double y)
	{
		y_cx=y;
		return IntArray[1]->getResult(&C_MultiDimInt::f2,this);;
	}
	T_out f2 (double z)
	{
		temp_vector[0]=y_cx;
		temp_vector[1]=z;
		return (obj_ref->*cxfunc)(temp_vector);
	}
	
	
	T_out BottomFunc (T_in temp)
	{
		return (obj_ref->*cxfunc)(temp);
	}
	
};




class C_Temp{
	public:
	t_cmplxArray2D Contour;
	
	//C_Integrator_Line<dcxMatrix,C_Temp,double> * integ_line_1;
	//C_Integrator_Line<dcxMatrix,C_Temp,double> * integ_line_2;
	//vector<C_Integrator_Line<dcxMatrix,C_Temp,double> *>  IntArray;
	C_MultiDimInt<t_cmplxMatrix,C_Temp,t_dArray1D> * multi_int;
	//C_Ololo<dcxMatrix,C_Temp,dArray1D> * ololo;
	//C_Integrator_Cauchy<dcxArray1D,dcxArray2D,dcx> * integ_cauchy;
	
	t_cmplxMatrix func (t_dArray1D a){
		t_cmplxMatrix temp(1,1);
		temp(0,0)=1.0;
		return temp;
	}
	
	void getResult(){
		cout << "GetResult" << endl;
		//integ=C_Integrator<dcxMatrix,C_Temp,double>::createIntegrator(40, 0.0001, 1.0, 1, qgauscheb_ID);
		//cout << integ_line->getResult(&C_Temp::func, this) << endl;
		//cout << integ_cauchy->getResult(&Contour,&point)[0] << endl;
		cout << multi_int->getResult(&C_Temp::func,this) << endl;
		//cout << ololo->getResult(&C_Temp::func,this) << endl;
	}
	
	void setContour(){
		t_cmplx ii(0.0,1.0);
		t_cmplx M2(0.0,1.0);
		t_dArray1D x,w;
		//integ_cauchy->getNodes(&x,&w);
		cout << x.size() << "  " << w.size() << endl;
		Contour.resize(3,t_cmplxArray1D());
		for (int i = 0; i < x.size(); i++){
			double t = x[i];
			Contour[0].push_back(exp(ii*t));
			Contour[1].push_back(ii*exp(ii*t));
			Contour[2].push_back(exp(ii*t));
		}
		
	
	}
	
	void SetIntegrators(){
		//integ_line_1=C_Integrator_Line<dcxMatrix,C_Temp,double>::createIntegrator(40, 0.0001, 1.0, 1, qgausleg_lin_ID);
		//integ_line_2=C_Integrator_Line<dcxMatrix,C_Temp,double>::createIntegrator(40, 0.0001, 1.0, 1, qgausleg_lin_ID);
		multi_int= new C_MultiDimInt<t_cmplxMatrix,C_Temp,t_dArray1D>(2);
		//ololo= new C_Ololo<dcxMatrix,C_Temp,dArray1D>();
		//multi_int->setIntegrators(IntArray,2);
		//integ_cauchy=C_Integrator_Cauchy<dcxArray1D,dcxArray2D,dcx>::createIntegrator(20000, 0.0, 6.2831, 1, qcauchyleg_lin_ID);
	}
	
	/*
	C_Temp * Copy(){
		return new C_Temp(*this);
	}*/
	

};

