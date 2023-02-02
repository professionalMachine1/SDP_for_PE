#pragma once
#include "AllMethods.h"

#ifndef METHODOUTSIDE_H
#define METHODOUTSIDE_H

class MethodOutside :
	public AllMethods
{
public:
	MethodOutside() : AllMethods() {}
	MethodOutside(int N, int M, TVector<double> XBorder, TVector<double> YBorder, int TASKNumber) : AllMethods(
		N, M, XBorder, YBorder, TASKNumber) {}

	virtual void MethodError(double eps, int MaxIterations); //ѕогрешность решени€, дл€ тестовых задач
};

#endif
