#include "MethodOutside.h"


void MethodOutside::MethodError(double eps, int MaxIterations)
{
	//Считаем точное решение
	TMatrix<double> u(m + 1);
	
	TestPartOne(u);

	double accurancy = eps + 1; //Погрешность решения и точость метода соответственно
	int IterationsCount = 0; //Количество итераций

	while ((accurancy > eps) && (IterationsCount < MaxIterations))
	{
		accurancy = Runner();
		IterationsCount++;
	}

	TestPartTwo(u);

	resultTest[0] = accurancy; //Точность метода
	resultTest[1] = IterationsCount; //Количество итераций

}
