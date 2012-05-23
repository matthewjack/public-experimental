//////////////////////////////////////////////////////////////////////
//
//	Crytek Common Source code
//	
//	File: Cry_Vector3.h
//	Description: Common vector class
//
//	History:
//	-Feb 27,2003: Created by Ivo Herzeg
//
//////////////////////////////////////////////////////////////////////

#ifndef VECTOR_H
#define VECTOR_H

#if _MSC_VER > 1000
# pragma once
#endif


#define ILINE _inline

#define IF(a, b) if((a))
#define WHILE(a, b) while((a))
#define IF_UNLIKELY(a) if((a))
#define IF_LIKELY(a) if((a))


#ifndef CRY_MATH_ARG_REF
	#ifndef PS3
		#define CRY_MATH_ARG_REF &


	#endif
#endif

// some constants
template<typename T> struct Vec3_tpl;

template<typename T> struct Vec3Constants
{
	static const Vec3_tpl<T> fVec3_Zero;
	static const Vec3_tpl<T> fVec3_OneX;
	static const Vec3_tpl<T> fVec3_OneY;
	static const Vec3_tpl<T> fVec3_OneZ;
	static const Vec3_tpl<T> fVec3_One;
};


template<typename T>
struct VecPrecisionValues
{
  ILINE static bool CheckGreater(const T value)
  {
    return value > 0;
  }
};

template<>
struct VecPrecisionValues<float>
{
  ILINE static bool CheckGreater(const float value)
  {
    return value > FLT_EPSILON;
  }
};

template <typename F> struct Vec3s_tpl
{
	F x,y,z;

	ILINE Vec3s_tpl( F vx, F vy, F vz ) : x(vx),y(vy),z(vz){}
	ILINE F &operator [] (int32 index)		  { assert(index>=0 && index<=2); return ((F*)this)[index]; }
	ILINE F operator [] (int32 index) const { assert(index>=0 && index<=2); return ((F*)this)[index]; }
};


//#ifdef XENON_INTRINSICS
//#include "Console\Cry_Xenon_Vector3.h"
//#else

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// class Vec3_tpl
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template <typename F> struct Vec3_tpl
{
	typedef F value_type;

	F x,y,z;

	ILINE Vec3_tpl()	{};

	/*!
	* template specailization to initialize a vector 
	* 
	* Example:
	*  Vec3 v0=Vec3(ZERO);
	*  Vec3 v1=Vec3(MIN);
	*  Vec3 v2=Vec3(MAX);
	*/
//	Vec3_tpl(type_zero) : x(0),y(0),z(0) {}
//	Vec3_tpl(type_min);
//	Vec3_tpl(type_max);

	/*!
	* constrctors and bracket-operator to initialize a vector 
	* 
	* Example:
	*  Vec3 v0=Vec3(1,2,3);
	*  Vec3 v1(1,2,3);
	*  v2.Set(1,2,3);
	*/
	ILINE Vec3_tpl( F vx, F vy, F vz ) : x(vx),y(vy),z(vz){ assert(this->IsValid()); }
	ILINE void operator () ( F vx, F vy, F vz ) { x=vx; y=vy; z=vz; assert(this->IsValid()); }
	ILINE Vec3_tpl<F>& Set(const F xval,const F yval, const F zval) { x=xval; y=yval; z=zval; assert(this->IsValid()); return *this; }

	explicit ILINE Vec3_tpl( F f ) : x(f),y(f),z(f) { assert(this->IsValid()); }

	/*!
	* the copy/casting/assignement constructor 
	* 
	* Example:
	*  Vec3 v0=v1;
	*/
	template <class T> ILINE  Vec3_tpl( const Vec3_tpl<T>& v ) : x((F)v.x), y((F)v.y), z((F)v.z) { assert(this->IsValid()); }
	//explicit ILINE Vec3_tpl(const Ang3_tpl<F>& v) : x((F)v.x), y((F)v.y), z((F)v.z) { assert(this->IsValid()); }
//	ILINE Vec3_tpl(const Vec2_tpl<F>& v) : x((F)v.x), y((F)v.y), z(0) { assert(this->IsValid()); }

	/*!
	* overloaded arithmetic operator  
	* 
	* Example:
	*  Vec3 v0=v1*4;
	*/
	ILINE Vec3_tpl<F> operator * (F k) const
	{
		const Vec3_tpl<F> v = *this;
		return Vec3_tpl<F>(v.x*k,v.y*k,v.z*k);
	}
	ILINE Vec3_tpl<F> operator / (F k) const
	{
		k=(F)1.0/k; return Vec3_tpl<F>(x*k,y*k,z*k);
	}
	ILINE friend Vec3_tpl<F> operator * (f32 f, const Vec3_tpl &vec)
	{
		return Vec3_tpl((F)(f*vec.x), (F)(f*vec.y), (F)(f*vec.z));
	}

	ILINE Vec3_tpl<F>& operator *= (F k)
	{
		//#ifndef XENON_INTRINSICS
		x*=k;y*=k;z*=k; return *this;
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMStoreFloat3(&m, X*k);
		//    return *this;
		//#endif
	}
	ILINE Vec3_tpl<F>& operator /= (F k)
	{
		//#ifndef XENON_INTRINSICS
		k=(F)1.0/k; x*=k;y*=k;z*=k; return *this;
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMStoreFloat3(&m, X/k);
		//    return *this;
		//#endif
	}

	ILINE Vec3_tpl<F> operator - ( void ) const
	{
		//#ifndef XENON_INTRINSICS
		return Vec3_tpl<F>(-x,-y,-z);
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    Vec3_tpl<F> res;
		//    XMStoreFloat3(&res.m, -X);
		//    return res;
		//#endif
	}
	ILINE Vec3_tpl<F>& Flip()
	{
		//#ifndef XENON_INTRINSICS
		x=-x; y=-y; z=-z; return *this;
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMStoreFloat3(&m, -X);
		//    return *this;
		//#endif
	}


	/*!
	* bracked-operator   
	* 
	* Example:
	*  uint32 t=v[0];
	*/
	ILINE F &operator [] (int32 index)		  { assert(index>=0 && index<=2); return ((F*)this)[index]; }
	ILINE F operator [] (int32 index) const { assert(index>=0 && index<=2); return ((F*)this)[index]; }


	/*!
	* functions and operators to compare vector   
	* 
	* Example:
	*  if (v0==v1) dosomething;
	*/
	ILINE bool operator==(const Vec3_tpl<F> &vec)
	{
		//#ifndef XENON_INTRINSICS
		return x == vec.x && y == vec.y && z == vec.z;
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR Y = XMLoadFloat3(&vec.m);
		//    return static_cast<bool>(XMVector3Equal(X, Y));
		//#endif
	}
	ILINE bool operator!=(const Vec3_tpl<F> &vec) { return !(*this == vec); }

	ILINE friend bool operator ==(const Vec3_tpl<F> &v0, const Vec3_tpl<F> &v1)
	{
		//#ifndef XENON_INTRINSICS
		return ((v0.x==v1.x) && (v0.y==v1.y) && (v0.z==v1.z));
		//#else
		//    XMVECTOR X = XMLoadFloat3(&v0.m);
		//    XMVECTOR Y = XMLoadFloat3(&v1.m);
		//    return static_cast<bool>(XMVector3Equal(X, Y));
		//#endif
	}
	ILINE friend bool operator !=(const Vec3_tpl<F> &v0, const Vec3_tpl<F> &v1)	{	return !(v0==v1);	}

	ILINE bool IsZero(F e = (F)0.0) const
	{
		//#ifndef XENON_INTRINSICS
		return  (fabs_tpl(x) <= e) && (fabs_tpl(y) <= e) && (fabs_tpl(z) <= e);
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    return static_cast<bool>(XMVector3Equal(X, XMVectorZero()));
		//#endif
	}

	// IsZeroSlow would be a better name [5/20/2010 evgeny]
	ILINE bool IsZeroFast(F e = (F)0.0003) const
	{
		//#ifndef XENON_INTRINSICS
		return (fabs_tpl(x) + fabs_tpl(y) + fabs_tpl(z)) <= e;
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    return static_cast<bool>(XMVector3Equal(X, XMVectorZero()));
		//#endif
	}


	ILINE bool IsEquivalent(const Vec3_tpl<F> CRY_MATH_ARG_REF v1, F epsilon=VEC_EPSILON) const 
	{
		assert(v1.IsValid()); 
		assert(this->IsValid()); 
		//#ifndef XENON_INTRINSICS
		return  ((fabs_tpl(x-v1.x) <= epsilon) &&	(fabs_tpl(y-v1.y) <= epsilon)&&	(fabs_tpl(z-v1.z) <= epsilon));	
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR Y = XMLoadFloat3(&v1.m);
		//    return static_cast<bool>(XMVector3NearEqual(X, Y, XMVectorReplicate(epsilon)));
		//#endif
	}

	ILINE bool IsUnit(F epsilon=VEC_EPSILON) const 
	{
		return (fabs_tpl(1 - GetLengthSquared()) <= epsilon);
	}

	bool IsValid() const
	{
		/*
		if (!NumberValid(x)) return false;
		if (!NumberValid(y)) return false;
		if (!NumberValid(z)) return false;
		*/
		return true;
	}

	//! force vector length by normalizing it
	ILINE void SetLength(F fLen)
	{ 
		//#ifndef XENON_INTRINSICS
		F fLenMe = GetLengthSquared();
		if(fLenMe<0.00001f*0.00001f)
			return;
		fLenMe = fLen * isqrt_tpl(fLenMe);
		x*=fLenMe; y*=fLenMe; z*=fLenMe;
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR fLenMe = XMVector3LengthSq(X);
		//    if (XMVector4Less(fLenMe, XMVectorReplicate(0.00001f*0.00001f)))
		//      return;
		//    //FIXME: do we need high precision here? If yes we have to use XMVectorReciprocalSqrt (slower)
		//    fLen *= XMVectorReciprocalSqrt(fLenMe).v[0];
		//    XMStoreFloat3(&m, X*fLen);
		//#endif
	}

	ILINE void ClampLength(F maxLength)
	{
		F sqrLength = GetLengthSquared();
		if (sqrLength > (maxLength * maxLength))
		{
			F scale = maxLength * isqrt_tpl(sqrLength);
			x *= scale; y *= scale; z *= scale;
		}
	}

	//! calculate the length of the vector
	ILINE F	GetLength() const
	{
		//#ifndef XENON_INTRINSICS
		return sqrt_tpl(x*x+y*y+z*z);
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    return XMVector3Length(X).v[0];
		//#endif
	}		

	ILINE F	GetLengthFloat() const
	{
		return GetLength();
	}

	ILINE F	GetLengthFast() const
	{
		//#ifndef XENON_INTRINSICS
		return sqrt_fast_tpl(x*x+y*y+z*z);
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    return XMVector3LengthEst(X).v[0];
		//#endif
	}		

	//! calculate the squared length of the vector
	ILINE F GetLengthSquared() const
	{
		//#ifndef XENON_INTRINSICS
		return x*x+y*y+z*z;
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    return XMVector3LengthSq(X).v[0];
		//#endif
	}

	ILINE F GetLengthSquaredFloat() const
	{
		return GetLengthSquared();
	}

	//! calculate the length of the vector ignoring the z component
	ILINE F	GetLength2D() const
	{
		//#ifndef XENON_INTRINSICS
		return sqrt_tpl(x*x+y*y);
		//#else
		//    XMVECTOR X = XMLoadVector2(&x);
		//    return XMVector2Length(X).v[0];
		//#endif
	}		

	//! calculate the squared length of the vector ignoring the z component
	ILINE F	GetLengthSquared2D() const
	{
		//#ifndef XENON_INTRINSICS
		return x*x+y*y;
		//#else
		//    XMVECTOR X = XMLoadVector2(&x);
		//    return XMVector2LengthSq(X).v[0];
		//#endif
	}		

	ILINE F GetDistance(const Vec3_tpl<F> vec1) const
	{
		//#ifndef XENON_INTRINSICS
		return  sqrt_tpl((x-vec1.x)*(x-vec1.x)+(y-vec1.y)*(y-vec1.y)+(z-vec1.z)*(z-vec1.z)); 
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR Y = XMLoadFloat3(&vec1.m);
		//    return XMVector3Length(X-Y).v[0];
		//#endif
	}
	ILINE F GetSquaredDistance ( const Vec3_tpl<F> &v) const
	{		
		//#ifndef XENON_INTRINSICS
		return  (x-v.x)*(x-v.x) + (y-v.y)*(y-v.y) + (z-v.z)*(z-v.z); 
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR Y = XMLoadFloat3(&v.m);
		//    return XMVector3LengthSq(X-Y).v[0];
		//#endif
	}
	ILINE F GetSquaredDistance2D ( const Vec3_tpl<F> &v) const
	{		
		//#ifndef XENON_INTRINSICS
		return  (x-v.x)*(x-v.x) + (y-v.y)*(y-v.y);
		//#else
		//    XMVECTOR X = XMLoadVector2(&x);
		//    XMVECTOR Y = XMLoadVector2(&v.x);
		//    return XMVector2LengthSq(X-Y).v[0];
		//#endif
	}

	//! normalize the vector
	// The default Normalize function is in fact "safe". 0 vectors remain unchanged.
	ILINE void	Normalize() 
	{
		assert(this->IsValid()); 
		//#ifndef XENON_INTRINSICS
		F fInvLen = isqrt_safe_tpl( x*x+y*y+z*z );
		x*=fInvLen; y*=fInvLen; z*=fInvLen; 
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    X = XMVector3Normalize(X);
		//    XMStoreFloat3(&m, X);
		//#endif
	}

	//! may be faster and less accurate
	ILINE void NormalizeFast() 
	{
		assert(this->IsValid()); 
		//#ifndef XENON_INTRINSICS
		F fInvLen = isqrt_fast_tpl( x*x+y*y+z*z );
		x*=fInvLen; y*=fInvLen; z*=fInvLen; 
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVector3NormalizeEst(X);
		//    XMStoreFloat3(&m, X);
		//#endif
	}

	//! normalize the vector
	// check for null vector - set to the passed in vector (which should be normalised!) if it is null vector
	// returns the original length of the vector
	ILINE F NormalizeSafe(const struct Vec3_tpl<F>& safe = Vec3Constants<F>::fVec3_Zero) 
	{ 
		assert(this->IsValid()); 
		//#ifndef XENON_INTRINSICS
		F fLen2 = x*x+y*y+z*z;
    IF (VecPrecisionValues<F>::CheckGreater(fLen2), 1)
		{
			F fInvLen = isqrt_tpl(fLen2);
			x*=fInvLen; y*=fInvLen; z*=fInvLen; 
			return F(1) / fInvLen;
		}
		else
		{
			*this = safe;
			return F(0);
		}
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR fLen2 = XMVector3LengthSq(X);
		//    if (fLen2.v[0]>0.0f)
		//    {
		//      XMVECTOR fInvLen = XMVectorReciprocalSqrt(fLen2);
		//      XMStoreFloat3(&m, X*fInvLen);
		//      return F(1) / fInvLen.v[0];
		//    }
		//    else
		//    {
		//      *this = safe;
		//      return F(0);
		//    }
		//#endif
	}

	ILINE Vec3_tpl GetNormalizedFloat() const
	{
		return GetNormalized();
	}

	//! return a normalized vector
	ILINE Vec3_tpl GetNormalized() const 
	{ 
		//#ifndef XENON_INTRINSICS
		F fInvLen = isqrt_safe_tpl( x*x+y*y+z*z );
		return *this * fInvLen;
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR Norm = XMVector3Normalize(X);
		//    Vec3_tpl <F> v;
		//    XMStoreFloat3(&v.m, Norm);
		//    return v;
		//#endif
	}

	//! return a normalized vector
	ILINE Vec3_tpl GetNormalizedFast() const 
	{ 
		F fInvLen = isqrt_fast_tpl( x*x+y*y+z*z );
		return *this * fInvLen;
	}

	/*
	//! return a safely normalized vector - returns safe vector (should be normalised) if original is zero length
	ILINE Vec3_tpl GetNormalizedSafe(const struct Vec3_tpl<F>& safe = Vec3Constants<F>::fVec3_OneX) const 
	{ 
		//#ifndef XENON_INTRINSICS
		F fLen2 = x*x+y*y+z*z;	
		IF (VecPrecisionValues<F>::CheckGreater(fLen2), 1)
		{
			F fInvLen = isqrt_tpl(fLen2);
			return *this * fInvLen;
		}
		else
		{
			return safe;
		}
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR fLen2 = XMVector3LengthSq(X);
		//    if (fLen2.v[0]>0.0f)
		//    {
		//      XMVECTOR fInvLen = XMVectorReciprocalSqrt(fLen2);
		//      Vec3_tpl <F> v;
		//      XMStoreFloat3(&v.m, X*fInvLen);
		//      return v;
		//    }
		//    else
		//    {
		//      return safe;
		//    }
		//#endif
	}
	*/

	//! return a safely normalized vector - returns safe vector (should be normalised) if original is zero length
	ILINE Vec3_tpl GetNormalizedSafeFloat(const struct Vec3_tpl<F>& safe = Vec3Constants<F>::fVec3_OneX) const 
	{
		return GetNormalizedSafe(safe);
	}

	// permutate coordinates so that z goes to new_z slot
	ILINE Vec3_tpl GetPermutated(int new_z) const { return Vec3_tpl(*(&x+inc_mod3[new_z]), *(&x+dec_mod3[new_z]), *(&x+new_z)); }

	// returns volume of a box with this vector as diagonal 
	ILINE F GetVolume() const { return x*y*z; }

	// returns a vector that consists of absolute values of this one's coordinates
	ILINE Vec3_tpl<F> abs() const
	{
		//#ifndef XENON_INTRINSICS
		return Vec3_tpl(fabs_tpl(x),fabs_tpl(y),fabs_tpl(z));
		//#else
		//    XMVECTOR X = XMLoadVector3(&x);
		//    XMVECTOR AbsV = XMVectorAbs(X);
		//    Vec3_tpl<F> v;
		//    XMStoreVector3(&v.x, AbsV);
		//    return v;
		//#endif
	}

	//! check for min bounds
	ILINE void CheckMin(const Vec3_tpl<F> CRY_MATH_ARG_REF other)
	{ 
		//#ifndef XENON_INTRINSICS
		x = min(other.x,x);
		y = min(other.y,y);
		z = min(other.z,z);
		/*
		if (other.x<x) x=other.x;
		if (other.y<y) y=other.y;
		if (other.z<z) z=other.z;
		*/
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR Y = XMLoadFloat3(&other.m);
		//    XMStoreFloat3(&m, XMVectorMin(X, Y));
		//#endif
	}			
	//! check for max bounds
	ILINE void CheckMax(const Vec3_tpl<F> &other)
	{
		//#ifndef XENON_INTRINSICS
		x = max(other.x,x);
		y = max(other.y,y);
		z = max(other.z,z);
		/*
		if (other.x>x) x=other.x;
		if (other.y>y) y=other.y;
		if (other.z>z) z=other.z;
		*/
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR Y = XMLoadFloat3(&other.m);
		//    XMStoreFloat3(&m, XMVectorMax(X, Y));
		//#endif
	}


	/*!
	* sets a vector orthogonal to the input vector 
	* 
	* Example:
	*  Vec3 v;
	*  v.SetOrthogonal( i );
	*/
	ILINE void SetOrthogonal( const Vec3_tpl<F>& v )
	{
		//#ifndef XENON_INTRINSICS
		int i = isneg(square((F)0.9)*v.GetLengthSquared()-v.x*v.x);
		(*this)[i]=0; (*this)[incm3(i)]= v[decm3(i)];	(*this)[decm3(i)]=-v[incm3(i)];
		/*#else
		XMVECTOR X = XMLoadVector3(&v.x);
		XMVECTOR Orth = XMVector3Orthogonal(X);
		XMStoreVector3(&x, Orth);
		#endif*/
	}
	// returns a vector orthogonal to this one
	ILINE Vec3_tpl GetOrthogonal() const
	{
		//#ifndef XENON_INTRINSICS
		int i = isneg(square((F)0.9)*GetLengthSquared()-x*x);
		Vec3_tpl<F> res;
		res[i]=0; res[incm3(i)]=(*this)[decm3(i)]; res[decm3(i)]=-(*this)[incm3(i)];
		return res;
		/*#else
		XMVECTOR X = XMLoadVector3(&x);
		XMVECTOR Orth = XMVector3Orthogonal(X);
		Vec3_tpl<F> res;
		XMStoreVector3(&res.x, Orth);
		return res;
		#endif*/
	}

	/*!
	* Project a point/vector on a (virtual) plane 
	* Consider we have a plane going through the origin. 
	* Because d=0 we need just the normal. The vector n is assumed to be a unit-vector.
	* 
	* Example:
	*  Vec3 result=Vec3::CreateProjection( incident, normal );
	*/
	ILINE void SetProjection( const Vec3_tpl& i, const Vec3_tpl& n ) { *this = i-n*(n|i); }
	ILINE static Vec3_tpl<F> CreateProjection( const Vec3_tpl& i, const Vec3_tpl& n ) { return i-n*(n|i); }

	/*!
	* Calculate a reflection vector. Vec3 n is assumed to be a unit-vector.
	* 
	* Example:
	*  Vec3 result=Vec3::CreateReflection( incident, normal );
	*/
	ILINE void SetReflection( const Vec3_tpl& i, const Vec3_tpl &n )
	{
		//#ifndef XENON_INTRINSICS
		*this=(n*(i|n)*2)-i;
		/*#else
		XMVECTOR I = XMLoadVector3(&i.x);
		XMVECTOR N = XMLoadVector3(&n.x);
		XMVECTOR ReflV = XMVector3Reflect(I, N);
		XMStoreVector3(&x, ReflV);
		#endif*/
	}
	ILINE static Vec3_tpl<F> CreateReflection( const Vec3_tpl& i, const Vec3_tpl &n )
	{
		//#ifndef XENON_INTRINSICS
		return (n*(i|n)*2)-i;
		/*#else
		XMVECTOR I = XMLoadVector3(&i.x);
		XMVECTOR N = XMLoadVector3(&n.x);
		XMVECTOR ReflV = XMVector3Reflect(I, N);
		Vec3_tpl<F> res;
		XMStoreVector3(&res.x, ReflV);
		return res;
		#endif*/
	}

	/*!
	* Linear-Interpolation between Vec3 (lerp)
	* 
	* Example:
	*  Vec3 r=Vec3::CreateLerp( p, q, 0.345f );
	*/
	ILINE void SetLerp( const Vec3_tpl<F> &p, const Vec3_tpl<F> &q, F t )
	{
		//#ifndef XENON_INTRINSICS
		const Vec3_tpl<F> diff = q-p;
		*this = p + (diff*t);
		//#else
		//    XMVECTOR P = XMLoadVector3(&p.x);
		//    XMVECTOR Q = XMLoadVector3(&q.x);
		//    XMVECTOR LerpV = XMVectorLerp(P, Q, t);
		//    XMStoreVector3(&x, LerpV);
		//#endif
	}
	ILINE static Vec3_tpl<F> CreateLerp( const Vec3_tpl<F> &p, const Vec3_tpl<F> &q, F t )
	{
		//#ifndef XENON_INTRINSICS
		const Vec3_tpl<F> diff = q-p;
		return p+(diff*t);
		//#else
		//    XMVECTOR P = XMLoadVector3(&p.x);
		//    XMVECTOR Q = XMLoadVector3(&q.x);
		//    XMVECTOR LerpV = XMVectorLerp(P, Q, t);
		//    Vec3_tpl<F> res;
		//    XMStoreVector3(&res.x, LerpV);
		//    return res;
		//#endif
	}


	/*!
	* Spherical-Interpolation between 3d-vectors (geometrical slerp)
	* both vectors are assumed to be normalized.  
	*
	* Example:
	*  Vec3 r=Vec3::CreateSlerp( p, q, 0.5674f );
	*/
	void SetSlerp( const Vec3_tpl<F>& p, const Vec3_tpl<F>& q, F t ) {
		assert(p.IsUnit(0.005f));
		assert(q.IsUnit(0.005f));
		// calculate cosine using the "inner product" between two vectors: p*q=cos(radiant)
		F cosine = clamp_tpl((p|q), F(-1), F(1));
		//we explore the special cases where the both vectors are very close together, 
		//in which case we approximate using the more economical LERP and avoid "divisions by zero" since sin(Angle) = 0  as   Angle=0
		if(cosine>=(F)0.99) {
			SetLerp(p,q,t); //perform LERP:
			Normalize();
		}	else {
			//perform SLERP: because of the LERP-check above, a "division by zero" is not possible
			F rad				= acos_tpl(cosine);
			F scale_0   = sin_tpl((1-t)*rad);
			F scale_1   = sin_tpl(t*rad);
			*this=(p*scale_0 + q*scale_1) / sin_tpl(rad);
			Normalize();
		}
	}
	ILINE static Vec3_tpl<F> CreateSlerp( const Vec3_tpl<F>& p, const Vec3_tpl<F>& q, F t ) {
		Vec3_tpl<F> v;	v.SetSlerp(p,q,t); return v;	
	}

	/*!
	* set random normalized vector (=random position on unit sphere) 
	* 
	* Example:
	*  Vec3 v;
	*  v.SetRandomDirection(); 
	*/
	void SetRandomDirection( void )
	{ 
		int nMax = 5; // To prevent hanging with bad randoms.
		F Length2;
		do {
			x = 1.0f - 2.0f*cry_frand();
			y = 1.0f - 2.0f*cry_frand();
			z = 1.0f - 2.0f*cry_frand();
			Length2 = len2();
			nMax--;
		} while((Length2>1.0f || Length2<0.0001f) && nMax > 0);
		F InvScale = isqrt_tpl(Length2);				// dividion by 0 is prevented
		x *= InvScale; y *= InvScale; z *= InvScale;
	}	


	/*!
	* rotate a vector using angle&axis 
	* 
	* Example:
	*  Vec3 r=v.GetRotated(axis,angle);
	*/
	ILINE Vec3_tpl GetRotated(const Vec3_tpl &axis, F angle) const
	{ 
		return GetRotated(axis,cos_tpl(angle),sin_tpl(angle)); 
	}
	ILINE Vec3_tpl GetRotated(const Vec3_tpl &axis, F cosa,F sina) const {
		Vec3_tpl zax = axis*(*this|axis); 
		Vec3_tpl xax = *this-zax; 
		Vec3_tpl yax = axis%xax;
		return xax*cosa + yax*sina + zax;
	}

	/*!
	* rotate a vector around a center using angle&axis 
	* 
	* Example:
	*  Vec3 r=v.GetRotated(axis,angle);
	*/
	ILINE Vec3_tpl GetRotated(const Vec3_tpl &center,const Vec3_tpl &axis, F angle) const { 
		return center+(*this-center).GetRotated(axis,angle); 
	}
	ILINE Vec3_tpl GetRotated(const Vec3_tpl &center,const Vec3_tpl &axis, F cosa,F sina) const { 
		return center+(*this-center).GetRotated(axis,cosa,sina); 
	}

	/*!
	* component wise multiplication of two vectors
	*/
	ILINE Vec3_tpl CompMul( const Vec3_tpl& rhs ) const { 
		return( Vec3_tpl( x * rhs.x, y * rhs.y, z * rhs.z ) ); 
	}

	//DEPRICATED ILINE friend F GetDistance(const Vec3_tpl<F> &vec1, const Vec3_tpl<F> &vec2) { 
	//return  sqrt_tpl((vec2.x-vec1.x)*(vec2.x-vec1.x)+(vec2.y-vec1.y)*(vec2.y-vec1.y)+(vec2.z-vec1.z)*(vec2.z-vec1.z)); 
	//}	
	//DEPRICATED ILINE friend F	GetSquaredDistance(const Vec3_tpl<F> &vec1, const Vec3_tpl<F> &vec2)	{		
	//return (vec2.x-vec1.x)*(vec2.x-vec1.x)+(vec2.y-vec1.y)*(vec2.y-vec1.y)+(vec2.z-vec1.z)*(vec2.z-vec1.z);
	//}
	//three methods for a "dot-product" operation
	ILINE F Dot (const Vec3_tpl<F> & v)	const
	{
		//#ifndef XENON_INTRINSICS
		const Vec3_tpl<F> vec2 = v;
		return x*vec2.x + y*vec2.y + z*vec2.z;
		//#else
		//    XMVECTOR X = XMLoadVector3(&x);
		//    XMVECTOR Y = XMLoadVector3(&vec2.x);
		//    return XMVector3Dot(X, Y).v[0];
		//#endif
	}
	//two methods for a "cross-product" operation
	ILINE Vec3_tpl<F> Cross (const Vec3_tpl<F> CRY_MATH_ARG_REF vec2) const
	{
		//#ifndef XENON_INTRINSICS
		return Vec3_tpl<F>( y*vec2.z  -  z*vec2.y,     z*vec2.x -    x*vec2.z,   x*vec2.y  -  y*vec2.x);
		/*#else
		XMVECTOR X = XMLoadVector3(&x);
		XMVECTOR Y = XMLoadVector3(&vec2.x);
		XMVECTOR CrossV = XMVector3Cross(X, Y);
		Vec3_tpl<F> res;
		XMStoreVector3(&res.x, CrossV);
		return res;
		#endif*/
	}	

	//f32* fptr=vec;
	DEPRICATED operator F* ()					{ return (F*)this; }
	template <class T>	explicit DEPRICATED Vec3_tpl(const T *src) { x=src[0]; y=src[1]; z=src[2]; }

	ILINE Vec3_tpl& zero() { x=y=z=0; return *this; }
	ILINE F len() const
	{
		//#ifndef XENON_INTRINSICS
		return sqrt_tpl(x*x+y*y+z*z);
		//#else
		//    XMVECTOR X = XMLoadVector3(&x);
		//    return XMVector3Length(X).v[0];
		//#endif
	}
	ILINE F len2() const
	{
		//#ifndef XENON_INTRINSICS
		return x*x +y*y + z*z;
		//#else
		//    XMVECTOR X = XMLoadVector3(&x);
		//    return XMVector3LengthSq(X).v[0];
		//#endif
	}

	ILINE Vec3_tpl& normalize()
	{ 
		//#ifndef XENON_INTRINSICS
		F len2 = x*x+y*y+z*z;
#if defined(XENON) || (defined(PS3) && !defined(__SPU__))
		const F rlen = isqrt_tpl(len2+1E-20f);
		x *= rlen; y *= rlen; z *= rlen;
		const F cRev = -len2;
		x = (F)__fsel(cRev, (F)0.f, x);
		y = (F)__fsel(cRev, (F)0.f, y);
		z = (F)__fsel(cRev, (F)1.f, z);
#else
		if (len2>(F)1e-20f) 
		{ 
			F rlen = isqrt_tpl(len2);
			x*=rlen; y*=rlen; z*=rlen; 
		} 
		else 
			Set(0,0,1); 
#endif
		return *this; 
		//#else
		//    XMVECTOR X = XMLoadVector3(&x);
		//    XMVECTOR NormV = XMVector3Normalize(X);
		//    XMStoreVector3(&x, NormV);
		//    return *this; 
		//#endif
	}
	ILINE Vec3_tpl normalized() const
	{ 
		//#ifndef XENON_INTRINSICS
		F len2 = x*x+y*y+z*z;
#if defined(XENON) || (defined(PS3) && !defined(__SPU__))
		Vec3_tpl ret;
		const F rlen = isqrt_tpl(len2+1E-20f);
		ret.x = x*rlen; ret.y = y*rlen; ret.z = z*rlen;
		const F cRev = -len2;
		ret.x = (F)__fsel(cRev, (F)0.f, ret.x);
		ret.y = (F)__fsel(cRev, (F)0.f, ret.y);
		ret.z = (F)__fsel(cRev, (F)1.f, ret.z);
		return ret;
#else
		if (len2>(F)1e-20f) 
		{ 
			F rlen = isqrt_tpl(len2); 
			return Vec3_tpl(x*rlen,y*rlen,z*rlen);
		}
		else 
			return Vec3_tpl(0,0,1);
#endif
		//#else
		//    XMVECTOR X = XMLoadVector3(&x);
		//    XMVECTOR NormV = XMVector3Normalize(X);
		//    Vec3_tpl res;
		//    XMStoreVector3(&res.x, NormV);
		//    return res;
		//#endif
	}

	// function overrides which should be used for uncached XMox memory!!!
	//vector subtraction
	template<class F1>
	ILINE Vec3_tpl<F1> sub(const Vec3_tpl<F1> &v) const
	{
		return Vec3_tpl<F1>(x-v.x, y-v.y, z-v.z);
	}

	//vector dot product
	template<class F1>
	ILINE F1 dot(const Vec3_tpl<F1> &v) const
	{
		return (F1)(x*v.x+y*v.y+z*v.z); 
	}
	//vector scale
	template<class F1>
	ILINE Vec3_tpl<F1> scale(const F1 k)
	{
		return Vec3_tpl<F>(x*k,y*k,z*k);
	}

	//vector cross product
	template<class F1>
	ILINE Vec3_tpl<F1> cross(const Vec3_tpl<F1> &v) const
	{
		return Vec3_tpl<F1>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); 
	}

	//AUTO_STRUCT_INFO
};


// dot product (2 versions)
template<class F1,class F2> 
ILINE F1 operator * (const Vec3_tpl<F1> &v0, const Vec3_tpl<F2> &v1)
{ 
	return v0.Dot(v1);
} 

template<class F1,class F2> 
ILINE F1 operator | (const Vec3_tpl<F1> &v0, const Vec3_tpl<F2> &v1)
{ 
	return v0.Dot(v1);
} 

// cross product
template<class F1,class F2> 
ILINE Vec3_tpl<F1> operator ^ (const Vec3_tpl<F1> &v0, const Vec3_tpl<F2> &v1)
{
	return v0.Cross(v1);
}

template<class F1,class F2> 
ILINE Vec3_tpl<F1> operator % (const Vec3_tpl<F1> &v0, const Vec3_tpl<F2> &v1)
{
	return v0.Cross(v1);
} 


template <class F> 
ILINE bool IsEquivalent(const Vec3_tpl<F> CRY_MATH_ARG_REF v0, const Vec3_tpl<F> CRY_MATH_ARG_REF v1, f32 epsilon=VEC_EPSILON )
{
	//#ifndef XENON_INTRINSICS
	return  ((fabs_tpl(v0.x-v1.x) <= epsilon) &&	(fabs_tpl(v0.y-v1.y) <= epsilon)&&	(fabs_tpl(v0.z-v1.z) <= epsilon));	
	//#else
	//  XMVECTOR X = XMLoadVector3(&v0.x);
	//  XMVECTOR Y = XMLoadVector3(&v1.x);
	//  return static_cast<bool>(XMVector3NearEqual(X, Y, XMVectorReplicate(epsilon)));
	//#endif
}
//---------------------------------------------------------------------------

//vector addition
template<class F1,class F2>
ILINE Vec3_tpl<F1> operator + (const Vec3_tpl<F1> CRY_MATH_ARG_REF v0, const Vec3_tpl<F2> CRY_MATH_ARG_REF v1)
{
	//	return Vec3_tpl<F1>(v0 + v1);
	//#ifndef XENON_INTRINSICS
	return Vec3_tpl<F1>(v0.x+v1.x, v0.y+v1.y, v0.z+v1.z);
	//#else
	//  XMVECTOR X = XMLoadVector3(&v0.x);
	//  XMVECTOR Y = XMLoadVector3(&v1.x);
	//  XMVECTOR AddV = XMVectorAdd(X, Y);
	//  Vec3_tpl<F1> res;
	//  XMStoreVector3(&res.x, AddV);
	//  return res;
	//#endif
}
/*
//vector addition
template<class F1,class F2>
ILINE Vec3_tpl<F1> operator + (const Vec2_tpl<F1> &v0, const Vec3_tpl<F2> &v1) {
	return Vec3_tpl<F1>(v0.x+v1.x, v0.y+v1.y, v1.z);
}
//vector addition
template<class F1,class F2>
ILINE Vec3_tpl<F1> operator + (const Vec3_tpl<F1> &v0, const Vec2_tpl<F2> &v1) {
	return Vec3_tpl<F1>(v0.x+v1.x, v0.y+v1.y, v0.z);
}
*/
//vector subtraction
template<class F1,class F2>
ILINE Vec3_tpl<F1> operator - (const Vec3_tpl<F1> & v0, const Vec3_tpl<F2> & r)
{
	//#ifndef XENON_INTRINSICS
	const Vec3_tpl<F2> v1 = r;
	return Vec3_tpl<F1>((F1)(v0.x-v1.x), (F1)(v0.y-v1.y), (F1)(v0.z-v1.z));
	//#else
	//  XMVECTOR X = XMLoadVector3(&v0.x);
	//  XMVECTOR Y = XMLoadVector3(&v1.x);
	//  XMVECTOR SubV = XMVectorSubtract(X, Y);
	//  Vec3_tpl<F1> res;
	//  XMStoreVector3(&res.x, SubV);
	//  return res;
	//#endif
}
/*
template<class F1,class F2>
ILINE Vec3_tpl<F1> operator - (const Vec2_tpl<F1> &v0, const Vec3_tpl<F2> &v1) {
	return Vec3_tpl<F1>(v0.x-v1.x, v0.y-v1.y, 0.0f-v1.z);
}
template<class F1,class F2>
ILINE Vec3_tpl<F1> operator - (const Vec3_tpl<F1> &v0, const Vec2_tpl<F2> &v1) {
	return Vec3_tpl<F1>(v0.x-v1.x, v0.y-v1.y, v0.z);
}
*/

//---------------------------------------------------------------------------


//vector self-addition
template<class F1,class F2>
ILINE Vec3_tpl<F1>& operator += (Vec3_tpl<F1> &v0, const Vec3_tpl<F2> & r)
{
	//#ifndef XENON_INTRINSICS
	const Vec3_tpl<F2> v1 = r;
  v0.x+=v1.x; v0.y+=v1.y; v0.z+=v1.z; return v0;
	//#else
	//  XMVECTOR X = XMLoadVector3(&v0.x);
	//  XMVECTOR Y = XMLoadVector3(&v1.x);
	//  XMVECTOR AddV = XMVectorAdd(X, Y);
	//  XMStoreVector3(&v0.x, AddV);
	//  return v0;
	//#endif
}

//vector self-subtraction
template<class F1,class F2>
ILINE Vec3_tpl<F1>& operator -= (Vec3_tpl<F1> &v0, const Vec3_tpl<F2> CRY_MATH_ARG_REF v1)
{
	//#ifndef XENON_INTRINSICS
  v0.x-=v1.x; v0.y-=v1.y; v0.z-=v1.z; return v0;
	//#else
	//  XMVECTOR X = XMLoadVector3(&v0.x);
	//  XMVECTOR Y = XMLoadVector3(&v1.x);
	//  XMVECTOR SubV = XMVectorSubtract(X, Y);
	//  XMStoreVector3(&v0.x, SubV);
	//  return v0;
	//#endif
}






///////////////////////////////////////////////////////////////////////////////
// Typedefs                                                                  //
///////////////////////////////////////////////////////////////////////////////
typedef Vec3_tpl<f32>		Vec3;			//we will use only this throughout the project
#ifndef REAL_IS_FLOAT
typedef Vec3_tpl<real>	Vec3r;		//for systems with float precision higher then 64bit
#else
typedef Vec3_tpl<f32>	Vec3r;		//for systems with float precision higher then 64bit
#endif
typedef Vec3_tpl<f64>	  Vec3_f64; //for double-precision
typedef Vec3_tpl<int>	  Vec3i;		//for integers
/*
template<> inline Vec3_tpl<f64>::Vec3_tpl(type_min) { x=y=z=-1.7E308; }
template<> inline Vec3_tpl<f64>::Vec3_tpl(type_max) { x=y=z=1.7E308; }
template<> inline Vec3_tpl<f32>::Vec3_tpl(type_min) { x=y=z=-3.3E38f; }
template<> inline Vec3_tpl<f32>::Vec3_tpl(type_max) { x=y=z=3.3E38f; }
*/

//////////////////////////////////////////////////////////////////////////
// Random vector functions.
/*
ILINE Vec3 Random(const Vec3& v)
{
	return Vec3( Random(v.x), Random(v.y), Random(v.z) );
}
ILINE Vec3 Random(const Vec3& a, const Vec3& b)
{
	return Vec3( Random(a.x,b.x), Random(a.y,b.y), Random(a.z,b.z) );
}

// Return a normalized copy of the vector flattened to the XY plane.
static inline Vec3 Vec3FlattenXY(const Vec3 & vIn)
{
	return Vec3(vIn.x,vIn.y,0.0f).GetNormalizedSafe(Vec3Constants<float>::fVec3_OneY);
}

*/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// class Vec4_tpl
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

































template <typename F> struct Vec4_tpl
{
	typedef F value_type;

	F x,y,z,w;

#if defined(_DEBUG) && !defined(__SPU__)
	ILINE Vec4_tpl() 
	{
		if (sizeof(F)==4)
		{
			uint32* p=alias_cast<uint32*>(&x);		p[0]=F32NAN;	p[1]=F32NAN; p[2]=F32NAN; p[3]=F32NAN;
		}
		if (sizeof(F)==8)
		{
			uint64* p=alias_cast<uint64*>(&x);		p[0]=F64NAN;	p[1]=F64NAN; p[2]=F64NAN; p[3]=F64NAN;
		}

	}
#else
	ILINE Vec4_tpl()	{CHECK_SIMD_ALIGNMENT};
#endif

	template<typename F2>
	ILINE Vec4_tpl<F>& operator = (const Vec4_tpl<F2> &v1) 
	{
		(*this).x=F(v1.x); (*this).y=F(v1.y); (*this).z=F(v1.z); (*this).w=F(v1.w);
		return (*this);
	}

	ILINE Vec4_tpl( F vx, F vy, F vz, F vw ) { x=vx; y=vy; z=vz; w=vw; };
	ILINE Vec4_tpl( const Vec3_tpl<F> &v, F vw ) {  x=v.x; y=v.y; z=v.z; w=vw; };

	ILINE void operator () ( F vx, F vy, F vz, F vw ) { x=vx; y=vy; z=vz; w=vw; };
	ILINE void operator () ( const Vec3_tpl<F> &v, F vw ) {  x=v.x; y=v.y; z=v.z; vw=vw; };

	ILINE F &operator [] (int index)		  { assert(index>=0 && index<=3); CHECK_SIMD_ALIGNMENT; return ((F*)this)[index]; }
	ILINE F operator [] (int index) const { assert(index>=0 && index<=3); CHECK_SIMD_ALIGNMENT; return ((F*)this)[index]; }
	template <class T> ILINE  Vec4_tpl( const Vec4_tpl<T>& v ) : x((F)v.x), y((F)v.y), z((F)v.z), w((F)v.w) { assert(this->IsValid()); }

	ILINE bool IsEquivalent(const Vec4_tpl<F>& v1, F epsilon=VEC_EPSILON) const 
	{
		assert(v1.IsValid()); 
		assert(this->IsValid()); 
		//#ifndef XENON_INTRINSICS
		return  ((fabs_tpl(x-v1.x) <= epsilon) &&	(fabs_tpl(y-v1.y) <= epsilon)&&	(fabs_tpl(z-v1.z) <= epsilon) && (fabs_tpl(w-v1.w) <= epsilon));	
		//#else
		//    XMVECTOR X = XMLoadFloat3(&m);
		//    XMVECTOR Y = XMLoadFloat3(&v1.m);
		//    return static_cast<bool>(XMVector3NearEqual(X, Y, XMVectorReplicate(epsilon)));
		//#endif
	}






























	ILINE Vec4_tpl<F> operator * (F k) const 
	{ 





		return Vec4_tpl<F>(x*k,y*k,z*k,w*k); 

	}
	ILINE Vec4_tpl<F> operator / (F k) const 
	{ 
		k=(F)1.0/k; 





		return Vec4_tpl<F>(x*k,y*k,z*k,w*k); 

	}

	ILINE Vec4_tpl<F>& operator *= (F k) 
	{ 




		x*=k;
		y*=k;
		z*=k;
		w*=k; 

		return *this; 
	}
	ILINE Vec4_tpl<F>& operator /= (F k) 
	{ 
		k=(F)1.0/k; 




		x*=k;		y*=k;		z*=k;		w*=k; 

		return *this; 
	}

	ILINE F Dot (const Vec4_tpl<F> &vec2)	const	
	{ 















		return x*vec2.x + y*vec2.y + z*vec2.z + w*vec2.w; 

	}
	ILINE F GetLength() const 
	{ 
		CHECK_SIMD_ALIGNMENT
			return sqrt_tpl(Dot(*this)); 
	}

	bool IsValid() const
	{
		CHECK_SIMD_ALIGNMENT
			if (!NumberValid(x)) return false;
		if (!NumberValid(y)) return false;
		if (!NumberValid(z)) return false;
		if (!NumberValid(w)) return false;
		return true;
	}

	/*!
	* functions and operators to compare vector   
	* 
	* Example:
	*  if (v0==v1) dosomething;
	*/
	ILINE bool operator==(const Vec4_tpl<F> &vec)
	{
		return x == vec.x && y == vec.y && z == vec.z && w == vec.w;
	}
	ILINE bool operator!=(const Vec4_tpl<F> &vec) { return !(*this == vec); }

	ILINE friend bool operator ==(const Vec4_tpl<F> &v0, const Vec4_tpl<F> &v1)
	{
		return ((v0.x==v1.x) && (v0.y==v1.y) && (v0.z==v1.z) && (v0.w==v1.w));
	}
	ILINE friend bool operator !=(const Vec4_tpl<F> &v0, const Vec4_tpl<F> &v1)	{	return !(v0==v1);	}

	//! normalize the vector
	// The default Normalize function is in fact "safe". 0 vectors remain unchanged.
	ILINE void	Normalize() 
	{
		assert(this->IsValid()); 
		F fInvLen = isqrt_safe_tpl( x*x + y*y + z*z + w*w );
		x*=fInvLen; y*=fInvLen; z*=fInvLen; w*= fInvLen;
	}

	//! may be faster and less accurate
	ILINE void NormalizeFast() 
	{
		assert(this->IsValid()); 
		F fInvLen = isqrt_fast_tpl( x*x + y*y + z*z + w*w );
		x*=fInvLen; y*=fInvLen; z*=fInvLen; w*= fInvLen;
	}

	ILINE void SetLerp( const Vec4_tpl<F> &p, const Vec4_tpl<F> &q, F t ) {	*this = p*(1.0f-t) + q*t;}
	ILINE static Vec4_tpl<F> CreateLerp( const Vec4_tpl<F> &p, const Vec4_tpl<F> &q, F t ) {	return p*(1.0f-t) + q*t;}


	//AUTO_STRUCT_INFO
} 



;


typedef Vec4_tpl<f32>		Vec4;			//we will use only this throughout the project
typedef DEFINE_ALIGNED_DATA(Vec4, Vec4A, 16); // typedef __declspec(align(16)) Vec4_tpl<f32>		Vec4A;			//we will use only this throughout the project
//before enabling it: add size check into Vec4 impl. for PS3
//typedef Vec4_tpl<f64>	  Vec4_f64; //for double-precision
typedef Vec4_tpl<real>	Vec4r;		//for systems with float precision higher then 64bit
//typedef Vec4_tpl<int>	  Vec4i;		//for integers


//vector self-addition
template<class F1,class F2>
ILINE Vec4_tpl<F1>& operator += (Vec4_tpl<F1> &v0, const Vec4_tpl<F2> CRY_MATH_ARG_REF v1) 
{





	v0.x+=v1.x; v0.y+=v1.y; v0.z+=v1.z; v0.w+=v1.w;

	return v0;
}

//vector addition
template<class F1,class F2>
ILINE Vec4_tpl<F1> operator + (const Vec4_tpl<F1> v0, const Vec4_tpl<F2> CRY_MATH_ARG_REF v1) 
{






	return Vec4_tpl<F1>(v0.x+v1.x, v0.y+v1.y, v0.z+v1.z, v0.w+v1.w);

}
//vector subtraction
template<class F1,class F2>
ILINE Vec4_tpl<F1> operator - (const Vec4_tpl<F1> CRY_MATH_ARG_REF v0, const Vec4_tpl<F2> CRY_MATH_ARG_REF v1) 
{






	return Vec4_tpl<F1>(v0.x-v1.x, v0.y-v1.y, v0.z-v1.z, v0.w-v1.w);

}

//vector multiplication
template<class F1,class F2>
ILINE Vec4_tpl<F1> operator * (const Vec4_tpl<F1> CRY_MATH_ARG_REF v0, const Vec4_tpl<F2> CRY_MATH_ARG_REF v1) 
{






	return Vec4_tpl<F1>(v0.x*v1.x, v0.y*v1.y, v0.z*v1.z, v0.w*v1.w);

}

//vector division
template<class F1,class F2>
ILINE Vec4_tpl<F1> operator / (const Vec4_tpl<F1> CRY_MATH_ARG_REF v0, const Vec4_tpl<F2> CRY_MATH_ARG_REF v1) 
{






	return Vec4_tpl<F1>(v0.x/v1.x, v0.y/v1.y, v0.z/v1.z, v0.w/v1.w);

}

//-----------------------------------------------------------------
// define the constants

template <typename T> const Vec3_tpl<T> Vec3Constants<T>::fVec3_Zero(0, 0, 0);
template <typename T> const Vec3_tpl<T> Vec3Constants<T>::fVec3_OneX(1, 0, 0);
template <typename T> const Vec3_tpl<T> Vec3Constants<T>::fVec3_OneY(0, 1, 0);
template <typename T> const Vec3_tpl<T> Vec3Constants<T>::fVec3_OneZ(0, 0, 1);
template <typename T> const Vec3_tpl<T> Vec3Constants<T>::fVec3_One(1, 1, 1);


//#endif



//reenable intrinsics




#endif //vector
