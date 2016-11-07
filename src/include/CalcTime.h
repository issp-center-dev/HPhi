#pragma once
#ifdef _OSX
#include <mach/clock.h>
#include <mach/mach.h>
#endif

void InitTimer();
void StartTimer(int n);
void StopTimer(int n);
void OutputTimer(struct BindStruct *X);
