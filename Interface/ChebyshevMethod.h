#pragma once
#include "AllMethods.h"
#include "MethodInside.h"

#ifndef CHEBYSHEVMETHOD_H
#define CHEBYSHEVMETHOD_H

using namespace std;

class ChebyshevMethod :
	public MethodInside
{
private:
	int K; //���������� ���������� 
	TVector<double> Tau;
	// optimal ����� �������� ��������� ���������� ����������
public:
	ChebyshevMethod() : MethodInside() {}
	ChebyshevMethod(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber) : MethodInside(
		N, M, XBorder, YBorder, TASKNumber) {
		Tau = TVector<double>(1);
		optimal = false;
		K = 8;
	}

	virtual double GetParametr(); //��������� �������� ���������
	virtual void SetUserParametr(double Parametr); //�������� ���������� ����������

	virtual void SetParametrs(); //�������� ������

	virtual void Runner(int index); //������������ �����

	virtual void MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name); //�������� �������, ��� �������� �����
	virtual void MethodError(double eps, int MaxIterations); //����������� �������, ��� �������� �����

};

#endif