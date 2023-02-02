#pragma once
#include "AllMethods.h"

#ifndef METHODINSIDE_H
#define METHODINSIDE_H

class MethodInside :
	public AllMethods
{
protected:
	bool optimal;
public:
	MethodInside() : AllMethods() {}
	MethodInside(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber) : AllMethods(
		N, M, XBorder, YBorder, TASKNumber) {
	}

	virtual void Optimal(bool wish); //��������� ��� ��������� ������������ ������������ ���������

	virtual void MethodError(double eps, int MaxIterations); //����������� �������, ��� �������� �����
};

#endif
