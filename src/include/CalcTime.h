#pragma once
#ifdef _OSX
#include <mach/clock.h>
#include <mach/mach.h>
#endif

void InitTimer();
inline void StartTimer(int n);
inline void StopTimer(int n);
void OutputTimer(struct BindStruct *X);
