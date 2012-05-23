// vectors.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <tchar.h>
#include <float.h>
#include <assert.h>

#include "Win32specific.h"
#include "Cry_Vector3.h"
#include "VecMath.h"
#include "Recast.h"

int _tmain(int argc, _TCHAR* argv[])
{
	Vec3 A(1,1,1);
	Vec3 A2(1,2,3);

	Vec3 A3 = A + A2;

	Vec3f B = {1,1,1};
	Vec3f B2 = {1,2,3};
	Vec3f B3 = B + B2;

	float C[3] = {1,1,1};
	float C2[3] = {1,2,3};
	float C3[3];

	rcVadd(C3, C, C2);

	return 0;
}

