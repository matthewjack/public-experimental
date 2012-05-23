#include <math.h>


struct Vec2f
{
	//POD type, do not add constructor, copy assignment, destructor, virtual functions
    union
    {
        struct
        {
            float x, y;
        };
        float v[2];
    };

	//casting
	operator float*() { return v; };
	operator const float*() { return v; };

	float Dot( const Vec2f& rhs_ ) const;
	float Length() const;
	float FastLength() const;
	Vec2f GetNormalized() const;
	Vec2f GetFastNormalized() const;


	//mutators
	Vec2f& Normalize();
	Vec2f& FastNormalize();

};

struct Vec3fRecast;

struct Vec3f
{
	//POD type, do not add constructor, copy assignment, destructor, virtual functions
    union
    {
        struct
        {
            float x, y, z;
        };
        float v[3];
    };

	//casting
	operator float*() { return v; };
	operator const float*() { return v; };

	float Dot( const Vec3f& rhs_ ) const;
	float Length() const;
	float FastLength() const;
	Vec3f GetNormalized() const;
	Vec3f GetFastNormalized() const;


	//mutators
	Vec3f& Normalize();
	Vec3f& FastNormalize();

	Vec3fRecast Transfrom() const;

};

struct Vec3fRecast
{
	float v[3];
	Vec3fRecast( const Vec3f& in_ )
	{
		v[0] = in_.x;
		v[1] = in_.z;
		v[2] = in_.y;
	}
	Vec3fRecast& operator=( const Vec3f& rhs_ )
	{

		v[0] = rhs_.x;
		v[1] = rhs_.z;
		v[2] = rhs_.y;
	}

	//casting
	operator float*() { return v; };
	operator const float*() { return v; };
	

};

Vec3fRecast Vec3f::Transfrom() const
{
	Vec3fRecast retval( *this );
	return retval;
}

//Fast Inverse Square Root
float FastInvSqrtf( float in_ )
{
	long i;
	float x2, y;
	const float threehalfs = 1.5F;
	
	x2 = in_ * 0.5f;
	y  = in_;
	i  = * ( long * ) &y;
	i  = 0x5f3759df - ( i >> 1 );
	y  = * ( float * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) );
	
	return y;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Vec3f Functions

//assignment operators
Vec2f& operator+=( Vec2f& lhs_, const Vec2f& rhs_ )
{
	lhs_.x += rhs_.x;
	lhs_.y += rhs_.y;
	return lhs_;
};

Vec2f& operator-=( Vec2f& lhs_, const Vec2f& rhs_ )
{
	lhs_.x -= rhs_.x;
	lhs_.y -= rhs_.y;
	return lhs_;
};

Vec2f& operator*=( Vec2f& lhs_, float rhs_ )
{
	lhs_.x *= rhs_;
	lhs_.y *= rhs_;
	return lhs_;
};

Vec2f& operator/=( Vec2f& lhs_, float rhs_ )
{
	lhs_.x /= rhs_;
	lhs_.y /= rhs_;
	return lhs_;
};

//binary operators
Vec2f operator+( const Vec2f& lhs_, const Vec2f& rhs_ )
{
	Vec2f result( lhs_ );
	result += rhs_;
	return result;
}

Vec2f operator-( const Vec2f& lhs_, const Vec2f& rhs_ )
{
	Vec2f result( lhs_ );
	result -= rhs_;
	return result;
}

//Dot Product
float Vec2f::Dot( const Vec2f& rhs_ ) const
{
	float result( x * rhs_.x + y * rhs_.y );
	return result;
};

Vec2f Vec2f::GetNormalized() const
{
	Vec2f result( *this );
	return result.Normalize();
}

Vec2f& Vec2f::Normalize()
{
	*this /= sqrtf( this->Dot( *this ) );
	return *this;
}

Vec2f& Vec2f::FastNormalize()
{
	*this *= FastInvSqrtf( this->Dot( *this ) );
	return *this;
}

Vec2f Vec2f::GetFastNormalized() const
{
	Vec2f result( *this );
	return result.FastNormalize();
}

float Vec2f::Length() const
{
	return sqrtf( this->Dot( *this ) );
}

float Vec2f::FastLength() const
{
	return 1.0f / FastInvSqrtf( this->Dot( *this ) );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Vec3f Functions

//assignment operators
Vec3f& operator+=( Vec3f& lhs_, const Vec3f& rhs_ )
{
	lhs_.x += rhs_.x;
	lhs_.y += rhs_.y;
	lhs_.z += rhs_.z;
	return lhs_;
};

Vec3f& operator-=( Vec3f& lhs_, const Vec3f& rhs_ )
{
	lhs_.x -= rhs_.x;
	lhs_.y -= rhs_.y;
	lhs_.z -= rhs_.z;
	return lhs_;
};

Vec3f& operator *=( Vec3f& lhs_, const Vec3f& rhs_ )
{
	lhs_.x = lhs_.y * rhs_.z - lhs_.z * rhs_.y;
	lhs_.y = lhs_.z * rhs_.x - lhs_.x * rhs_.z;
	lhs_.z = lhs_.x * rhs_.y - lhs_.y * rhs_.x;
	return lhs_;
}

Vec3f& operator/=( Vec3f& lhs_, float rhs_ )
{
	lhs_.x /= rhs_;
	lhs_.y /= rhs_;
	lhs_.z /= rhs_;
	return lhs_;
};

//cross product assignment
Vec3f& operator*=( Vec3f& lhs_, float rhs_ )
{
	lhs_.x *= rhs_;
	lhs_.y *= rhs_;
	lhs_.z *= rhs_;
	return lhs_;
};


//binary operators
Vec3f operator+( const Vec3f& lhs_, const Vec3f& rhs_ )
{
	Vec3f result( lhs_ );
	result += rhs_;
	return result;
}

Vec3f operator-( const Vec3f& lhs_, const Vec3f& rhs_ )
{
	Vec3f result( lhs_ );
	result -= rhs_;
	return result;
}

//Cross product
Vec3f operator*( const Vec3f& lhs_, const Vec3f& rhs_ )
{
	Vec3f result( lhs_ );
	result *= rhs_;
	return result;
}

//Dot Product
float Vec3f::Dot( const Vec3f& rhs_ ) const
{
	float result( x * rhs_.x + y * rhs_.y + z * rhs_.z );
	return result;
};

float Vec3f::Length() const
{
	return sqrtf( Dot( *this ) );
}

float Vec3f::FastLength() const
{
	return 1.0f / FastInvSqrtf( Dot( *this ) );
}

Vec3f Vec3f::GetNormalized() const
{
	Vec3f result( *this );
	return result.Normalize();
}

Vec3f Vec3f::GetFastNormalized() const
{
	Vec3f result( *this );
	return result.FastNormalize();
}
Vec3f& Vec3f::Normalize()
{
	*this /= sqrtf( Dot( *this ) );
	return *this;
}

Vec3f& Vec3f::FastNormalize()
{
	*this *= FastInvSqrtf( Dot( *this ) );
	return *this;
}