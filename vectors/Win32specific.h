////////////////////////////////////////////////////////////////////////////
//
//  Crytek Engine Source File.
//  Copyright (C), Crytek Studios, 2002.
// -------------------------------------------------------------------------
//  File name:   Win32specific.h
//  Version:     v1.00
//  Created:     31/03/2003 by Sergiy.
//  Compilers:   Visual Studio.NET
//  Description: Specific to Win32 declarations, inline functions etc.
// -------------------------------------------------------------------------
//  History:
//
////////////////////////////////////////////////////////////////////////////

#ifndef _CRY_COMMON_WIN32_SPECIFIC_HDR_
#define _CRY_COMMON_WIN32_SPECIFIC_HDR_

#define _CPU_X86
#define _CPU_SSE
#ifdef _DEBUG
#define ILINE _inline
#else
#define ILINE __forceinline
#endif

#define DEBUG_BREAK _asm { int 3 }
#define RC_EXECUTABLE "rc.exe"
#define DEPRICATED __declspec(deprecated)
#define TYPENAME(x) typeid(x).name()
#define SIZEOF_PTR 4

#ifndef _WIN32_WINNT
# define _WIN32_WINNT 0x501
#endif

//////////////////////////////////////////////////////////////////////////
// Standart includes.
//////////////////////////////////////////////////////////////////////////
#include <malloc.h>
#include <io.h>
#include <new.h>
#include <direct.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
//////////////////////////////////////////////////////////////////////////

// Special intrinsics
#include <intrin.h> // moved here to prevent assert from being redefined when included elsewhere


//////////////////////////////////////////////////////////////////////////
// Define platform independent types.
//////////////////////////////////////////////////////////////////////////
typedef signed char         int8;
typedef signed short        int16;
typedef signed int					int32;
typedef signed __int64			int64;
typedef unsigned char				uint8;
typedef unsigned char				BYTE;
typedef unsigned short			uint16;
typedef unsigned int				uint32;
typedef unsigned __int64		uint64;
typedef float               f32;
typedef double              f64;
typedef double              real;  //biggest float-type on this machine

typedef unsigned long       DWORD;  //biggest float-type on this machine

#if defined(_WIN64)
    typedef __int64 INT_PTR, *PINT_PTR;
    typedef unsigned __int64 UINT_PTR, *PUINT_PTR;

    typedef __int64 LONG_PTR, *PLONG_PTR;
    typedef unsigned __int64 ULONG_PTR, *PULONG_PTR;

		typedef ULONG_PTR DWORD_PTR, *PDWORD_PTR;
#else
    typedef __w64 int INT_PTR, *PINT_PTR;
    typedef __w64 unsigned int UINT_PTR, *PUINT_PTR;

    typedef __w64 long LONG_PTR, *PLONG_PTR;
    typedef __w64 unsigned long ULONG_PTR, *PULONG_PTR;

		typedef ULONG_PTR DWORD_PTR, *PDWORD_PTR;
#endif

typedef void *THREAD_HANDLE;
typedef void *EVENT_HANDLE;

//////////////////////////////////////////////////////////////////////////
// Multi platform Hi resolution ticks function, should only be used for profiling.
//////////////////////////////////////////////////////////////////////////

#if defined(_CPU_SSE) && !defined(NO_XENON_INTRINSICS)
//# define XENON_INTRINSICS
//# include "Cry_XMath.h"
#endif //defined(_CPU_SSE) && !defined(NO_XENON_INTRINSICS)

int64 CryGetTicks();

#ifndef SAFE_DELETE
#define SAFE_DELETE(p)			{ if(p) { delete (p);		(p)=NULL; } }
#endif

#ifndef SAFE_DELETE_ARRAY
#define SAFE_DELETE_ARRAY(p)	{ if(p) { delete [] (p);		(p)=NULL; } }
#endif

#ifndef SAFE_RELEASE
#define SAFE_RELEASE(p)			{ if(p) { (p)->Release();	(p)=NULL; } }
#endif

#define _MS_ALIGN(num) __declspec(align(num))

#define DEFINE_ALIGNED_DATA( type, name, alignment ) _declspec(align(alignment)) type name;
#define DEFINE_ALIGNED_DATA_STATIC( type, name, alignment ) static _declspec(align(alignment)) type name;
#define DEFINE_ALIGNED_DATA_CONST( type, name, alignment ) const _declspec(align(alignment)) type name;

#ifndef FILE_ATTRIBUTE_NORMAL 
	#define FILE_ATTRIBUTE_NORMAL 0x00000080
#endif

#if !defined(CRYTEK_VISUALIZATION)
#define FP16_MESH 1
#endif 

#define TARGET_DEFAULT_ALIGN (0x4U)

//#define ALLOW_CAMERA_SPACE 1  // Allow camera-space objects (for Kiew project)

#endif //_CRY_COMMON_WIN32_SPECIFIC_HDR_
