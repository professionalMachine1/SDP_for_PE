#include "MPI.h"


void MPI::SetUserParametr(double Parametr)
{
	Tau = Parametr;
}


double MPI::GetParametr()
{
	return Tau;
}


void MPI::SetParametrs()
{
	double Max, Min;
	Min = -4 * hE*pow(sin(pi / (2.0*n)), 2) - 4 * kE*pow(sin(pi / (2.0*m)), 2);
	Max = -4 * hE*pow(cos(pi / (2.0*n)), 2) - 4 * kE*pow(cos(pi / (2.0*m)), 2);

	Tau = 2 / (Max + Min);
}


double MPI::Runner()
{
	double temp;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
			R[j][i] = A * V[j][i] + hE * (V[j][i + 1] + V[j][i - 1]) + kE * (V[j + 1][i] + V[j - 1][i]) + F[j][i];

	double accurancy = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			temp = V[j][i];
			//
			V[j][i] = V[j][i] - Tau * R[j][i];
			//
			temp = fabs(V[j][i] - temp);
			if (temp > accurancy)
				accurancy = temp;
		}
	return accurancy;
}


void MPI::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name)
{
	AllMethods* Solution;
	Solution = new MPI(2 * n, 2 * m, xBorder, yBorder, TaskNumber);

	MainPartOne(Solution, Name);

	if (optimal)
	{
		SetParametrs();
		Solution->SetParametrs();
	}
	else {
		Solution->Optimal(false);
		Solution->SetUserParametr(Tau);
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


