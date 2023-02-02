#include "ChebyshevMethod.h"


void ChebyshevMethod::SetUserParametr(double Parametr)
{
	K = (int)Parametr;
}


double ChebyshevMethod::GetParametr()
{
	return K;
}


void ChebyshevMethod::SetParametrs()
{
	double Max, Min;
	Min = -4 * hE*pow(sin(pi / (2.0*n)), 2) - 4 * kE*pow(sin(pi / (2.0*m)), 2);
	Max = -4 * hE*pow(cos(pi / (2.0*n)), 2) - 4 * kE*pow(cos(pi / (2.0*m)), 2);

	Tau = TVector<double>(K);

	int step = K / 2;
	TVector<int> sup(K);
	//
	if (optimal)
	{
		for (int i = 2; i <= K; i *= 2)
		{
			for (int j = step; j < K; j += 2 * step)
			{
				sup[j] = (i - 1) - sup[j - step];
			}
			step /= 2;
		}
	}
	else {
		for (int i = 0; i < K; i++)
			sup[i] = i;
	}

	//
	for (int i = 0; i < K; i++)
	{
		Tau[i] = 2 / (Max + Min + (Max - Min)*cos(pi*(2 * sup[i] + 1) / (2.0*K)));
	}
}


void ChebyshevMethod::Runner(int index)
{
	double temp;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
			R[j][i] = A * V[j][i] + hE * (V[j][i + 1] + V[j][i - 1]) + kE * (V[j + 1][i] + V[j - 1][i]) + F[j][i];


	double accurancy = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
			V[j][i] = V[j][i] - Tau[index] * R[j][i];

}


void ChebyshevMethod::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name)
{
	AllMethods* Solution;
	Solution = new ChebyshevMethod(2 * n, 2 * m, xBorder, yBorder, TaskNumber);

	MainPartOne(Solution, Name);

	SetParametrs();
	if (optimal)
		Solution->Optimal(true);
	else
		Solution->Optimal(false);
	Solution->SetUserParametr(K);
	Solution->SetParametrs();

	TVector<double> accurancy(2); //Точность метода
	//Чтобы зайти в цикл
	for (int i = 0; i < 2; i++)
		accurancy[i] = 1 + eps[i];

	TVector<int> IterationsCount(2);
	//
	int TauIndex = 0;
	TMatrix<double> Matr = V;
	while ((accurancy[0] > eps[0]) && (IterationsCount[0] < MaxIterations[0]))
	{
		if (TauIndex < K - 1)
		{
			Runner(TauIndex);
			TauIndex++;
		}
		else
		{
			Runner(TauIndex);
			TauIndex = 0;
			IterationsCount[0]++;
			accurancy[0] = CalculateAccurancy(Matr);
			Matr = V;
		}
	}
	TauIndex = 0;
	Matr = Solution->GetV();
	while ((accurancy[1] > eps[1]) && (IterationsCount[1] < MaxIterations[1]))
	{
		if (TauIndex < K - 1)
		{
			Solution->Runner(TauIndex);
			TauIndex++;
		}
		else
		{
			Solution->Runner(TauIndex);
			TauIndex = 0;
			IterationsCount[1]++;
			accurancy[1] = Solution->CalculateAccurancy(Matr);
			Matr = Solution->GetV();
		}
	}
	//

	MainPartTwo(Solution, Solution->GetV());

	resultMain[0] = accurancy[0]; //Точность метода на сетке (n+1,m+1)
	resultMain[1] = IterationsCount[0]; //Количество итераций на сетке (n+1,m+1)
	resultMain[2] = accurancy[1]; //Точность метода на сетке (2n+1,2m+1)
	resultMain[3] = IterationsCount[1]; //Количество итераций на сетке (2n+1,2m+1)

}


void ChebyshevMethod::MethodError(double eps, int MaxIterations)
{
	//Считаем точное решение
	TMatrix<double> u(m + 1);
	
	TestPartOne(u);

	double accurancy = eps + 1; //Погрешность решения и точость метода соответственно
	int IterationsCount = 0; //Количество итераций

	SetParametrs();

	int TauIndex = 0;
	TMatrix<double> Matr = V;
	while ((accurancy > eps) && (IterationsCount < MaxIterations))
	{
		if (TauIndex < K - 1)
		{
			Runner(TauIndex);
			TauIndex++;
		}
		else
		{
			Runner(TauIndex);
			TauIndex = 0;
			IterationsCount++;
			accurancy = CalculateAccurancy(Matr);
			Matr = V;
		}
	}

	TestPartTwo(u);

	resultTest[0] = accurancy; //Точность метода
	resultTest[1] = IterationsCount; //Количество итераций

}
