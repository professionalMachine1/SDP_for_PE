#include "TopRelaxMethod.h"


void TopRelaxMethod::SetUserParametr(double Parametr)
{
	w = Parametr;
}


double TopRelaxMethod::GetParametr()
{
	return w;
}


void TopRelaxMethod::SetParametrs()
{
	double temp1, temp2;
	temp1 = pow(sin(pi / (2.0*n)), 2);
	temp2 = pow(sin(pi / (2.0*m)), 2);

	temp2 = 2 * pow(k, 2)*temp1 + 2 * pow(h, 2)*temp2;
	temp1 = pow(h, 2) + pow(k, 2);
	temp2 = temp2 / temp1;
	w = 2 / (1 + sqrt(temp2*(2 - temp2))); //Оптимальный параметр
	
}


double TopRelaxMethod::Runner()
{
	double temp, prev;
	double errorM = 0; //Точность метода
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			prev = V[j][i];
			//
			temp = A * prev + hE * (V[j][i - 1] + V[j][i + 1]) + kE * (V[j - 1][i] + V[j + 1][i]);
			V[j][i] = prev - w * (temp + F[j][i]) / A;
			//
			temp = fabs(V[j][i] - prev);
			if (temp > errorM)
				errorM = temp;
		}
	return errorM;
}


void TopRelaxMethod::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name)
{
	AllMethods* Solution;
	Solution = new TopRelaxMethod(2 * n, 2 * m, xBorder, yBorder, TaskNumber);

	MainPartOne(Solution, Name);

	if (optimal)
	{
		SetParametrs();
		Solution->SetParametrs();
	}
	else {
		Solution->Optimal(false);
		Solution->SetUserParametr(w);
	}

	TVector<double> accurancy(2); //Точность метода
	//Чтобы зайти в цикл
	for (int i = 0; i < 2; i++)
		accurancy[i] = 1 + eps[i];

	TVector<int> IterationsCount(2);
	while ((accurancy[0] > eps[0]) && (IterationsCount[0] < MaxIterations[0]))
	{
		accurancy[0] = Runner();
		IterationsCount[0]++;
	}
	while ((accurancy[1] > eps[1]) && (IterationsCount[1] < MaxIterations[1]))
	{
		accurancy[1] = Solution->Runner();
		IterationsCount[1]++;
	}

	MainPartTwo(Solution, Solution->GetV());
	
	resultMain[0] = accurancy[0]; //Точность метода на сетке (n+1,m+1)
	resultMain[1] = IterationsCount[0]; //Количество итераций на сетке (n+1,m+1)
	resultMain[2] = accurancy[1]; //Точность метода на сетке (2n+1,2m+1)
	resultMain[3] = IterationsCount[1]; //Количество итераций на сетке (2n+1,2m+1)
	
}

