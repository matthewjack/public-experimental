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


class Vec3Adaptor
{
	Vec3Adaptor() {};
};

/*
template <class C> add( C &result, const C &a, const C &b)
{
	result.x = a.x + b.x;
	result.y = a.y + b.y;
	result.z = a.z + b.z;
}

template <class C, class D> add( C &result, const C &a, const D &b)
{
	result.x = a.x + b.x;
	result.y = a.y + b.z;
	result.z = a.z + b.x;
}
*/

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

	Vec3fRecast D = B;
	rcVadd(C, C2, B.Transfrom());



	return 0;
}
