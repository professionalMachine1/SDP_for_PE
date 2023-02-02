#include "MSG.h"


void MSG::SetParametrs()
{
	double b = 0;
	double temp = 0, res = 0;

	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			res = A * H[j][i] + hE * (H[j][i + 1] + H[j][i - 1]) + kE * (H[j + 1][i] + H[j - 1][i]);
			temp += res * R[j][i];
		}

	b = temp / Ahh;

	H = R * (-1) + H * b;

	Ahh = 0; temp = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			res = A * H[j][i] + hE * (H[j][i + 1] + H[j][i - 1]) + kE * (H[j + 1][i] + H[j - 1][i]);
			temp += R[j][i] * H[j][i];
			Ahh += res * H[j][i];
		}
	a = -temp / Ahh;
	
}


double MSG::Runner()
{
	double temp;
	bet = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			R[j][i] = A * V[j][i] + hE * (V[j][i + 1] + V[j][i - 1]) + kE * (V[j + 1][i] + V[j - 1][i]) + F[j][i];
			bet += fabs(R[j][i]);
		}

	double accurancy = 0;
	SetParametrs();
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			temp = V[j][i];
			//
			V[j][i] = V[j][i] + a * H[j][i];
			//
			temp = fabs(V[j][i] - temp);
			if (temp > accurancy)
				accurancy = temp;
		}
	return accurancy;
}


void MSG::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name)
{
	AllMethods* Solution;
	Solution = new MSG(2 * n, 2 * m, xBorder, yBorder, TaskNumber);
	
	MainPartOne(Solution, Name);

	double sup; //Вспомогательная переменная
	TVector<double> accurancy(2); //Точность метода
	//Чтобы зайти в цикл
	for (int i = 0; i < 2; i++)
		accurancy[i] = 1 + eps[i];

	TVector<int> IterationsCount(2);

	//
	while ((accurancy[0] > eps[0]) && (IterationsCount[0] < MaxIterations[0]))
	{
		sup = Runner();
		if (bet <= pow(10, -12))
			break;
		accurancy[0] = sup;
		IterationsCount[0]++;
	}

	while ((accurancy[1] > eps[1]) && (IterationsCount[1] < MaxIterations[1]))
	{
		sup = Solution->Runner();
		if (Solution->GetBet() <= pow(10, -12))
			break;
		accurancy[1] = sup;
		IterationsCount[1]++;
	}
	
	MainPartTwo(Solution, Solution->GetV());
		
	resultMain[0] = accurancy[0]; //Точность метода на сетке (n+1,m+1)
	resultMain[1] = IterationsCount[0]; //Количество итераций на сетке (n+1,m+1)
	resultMain[2] = accurancy[1]; //Точность метода на сетке (2n+1,2m+1)
	resultMain[3] = IterationsCount[1]; //Количество итераций на сетке (2n+1,2m+1)

}


void MSG::MethodError(double eps, int MaxIterations)
{
	//Считаем точное решение
	TMatrix<double> u(m + 1);
	
	TestPartOne(u);

	double accurancy = eps + 1; //Погрешность решения и точость метода соответственно
	double sup; //Переменная помощник
	int IterationsCount = 0; //Количество итераций

	while ((accurancy > eps) && (IterationsCount < MaxIterations))
	{
		sup = Runner();
		if (bet <= pow(10, -12))
			break;
		accurancy = sup;
		IterationsCount++;
	}

	TestPartTwo(u);

	resultTest[0] = accurancy; //Точность метода
	resultTest[1] = IterationsCount; //Количество итераций

}
