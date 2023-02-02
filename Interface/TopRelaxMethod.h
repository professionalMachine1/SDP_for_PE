#pragma once
#include "AllMethods.h"
#include "MethodInside.h"

#ifndef TOPRELAXMETHOD_H
#define TOPRELAXMETHOD_H

using namespace std;

class TopRelaxMethod :
	public MethodInside
{
private:
	double w;
public:
	TopRelaxMethod() : MethodInside() {}
	TopRelaxMethod(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber) : MethodInside(
		N, M, XBorder, YBorder, TASKNumber) {
		w = 1;
		optimal = true;
	}

	virtual void SetUserParametr(double Parametr); //�������� ���������������� ��������
	virtual double GetParametr(); //��������� �������� ���������

	virtual void SetParametrs(); //�������� ������

	virtual double Runner(); //������������ �����

	virtual void MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name);

};

#endif