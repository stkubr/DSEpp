

/// Prints a line accross the screen
void PrintLine(char __c){ 
  for(int i = 0; i < 80; i++) cout << __c; 
  cout << endl;
}

/// Prints "h" spaces 
void PrintSpaces(int h){ 
  for(int i = 0; i < h; i++) cout << " "; 
}


t_cmplx Cheb_polinoms(t_cmplx x, int order)
{
	t_cmplx U[order+10];
	
	U[0]=1.0;
	U[1]=2.0*x;
	U[2]=4.0*x*x - 1.0;
	
	if(order==0)
	{
		return U[0];
	}
	
	if(order==1)
	{
		return U[1];
	}
		
	if(order==2)
	{
		return U[2];
	}
		
	if(order>2)
	{
		for (int i = 3; i <= order; i++)
		{
			U[i]=2.0*x*U[i-1] - U[i-2];
		}
		return U[order];
	}
	
	cout << "Chebys Error!" << endl;
	return 0;
}
/*
double SignFlap(int k)
{
	dcx I,res;
	res=dcx(0.0,0.0);
	I=dcx(0.0,1.0);
	
	if(k==0) return 1;
	if(k>0)
	{
		for (int i = 1; i < k; i++)
		{
			res
		}
		
	}
	
	
}*/


double Get_Time()
{
	return (double)clock()/CLOCKS_PER_SEC;
}

void DebugLine(string _word)
{
//#if DEBUG_MODE==1
	cout <<" Debug Line at place - "<< _word << endl;
//#endif
}

void StopLine()
{
	exit(1);
}

