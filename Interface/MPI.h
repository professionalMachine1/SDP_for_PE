#pragma once
#include "AllMethods.h"
#include "MethodInside.h"

#ifndef MPI_H
#define MPI_H

using namespace std;

class MPI :
	public MethodInside
{
private:
	double Tau;
public:
	MPI() : MethodInside() {}
	MPI(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber) : MethodInside(
		N, M, XBorder, YBorder, TASKNumber) {
		Tau = 0;
		optimal = true;
	}

	virtual void SetUserParametr(double Parametr); //�������� ���������������� ��������
	virtual double GetParametr(); //��������� �������� ���������

	virtual void SetParametrs(); //�������� ������
	virtual double Runner(); //������������ �����

	virtual void MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name);

};

#endif