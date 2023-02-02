#pragma once
#include "AllMethods.h"
#include "MethodOutside.h"

#ifndef MMN_H
#define MMN_H

using namespace std;


class MMN :
	public MethodOutside
{
private:
	double Tau; //�������� ������
public:
	MMN() : MethodOutside() {}
	MMN(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber) : MethodOutside(
		N, M, XBorder, YBorder, TASKNumber) {
		Tau = 0;
	}

	virtual void SetParametrs(); //�������� ������
	virtual double Runner(); //������������ �����

	virtual void MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name);
};

#endif