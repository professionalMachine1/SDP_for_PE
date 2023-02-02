#include "MethodOutside.h"


void MethodOutside::MethodError(double eps, int MaxIterations)
{
	//������� ������ �������
	TMatrix<double> u(m + 1);
	
	TestPartOne(u);

	double accurancy = eps + 1; //����������� ������� � ������� ������ ��������������
	int IterationsCount = 0; //���������� ��������

	while ((accurancy > eps) && (IterationsCount < MaxIterations))
	{
		accurancy = Runner();
		IterationsCount++;
	}

	TestPartTwo(u);

	resultTest[0] = accurancy; //�������� ������
	resultTest[1] = IterationsCount; //���������� ��������

}
