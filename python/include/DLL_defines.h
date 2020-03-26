#ifndef _ctmm_DLLDEFINES_H_
#define _ctmm_DLLDEFINES_H_

#if defined (_WIN32)
    #if defined (ctmm_EXPORTS)
        #define CTMM_EXPORT __declspec(dllexport)
    #else
        #define CTMM_EXPORT __declspec(dllimport)
    #endif
#else
    #define CTMM_EXPORT
#endif

#endif