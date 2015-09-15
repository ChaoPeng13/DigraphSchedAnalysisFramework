#include "Timer.h"
//#include <Windows.h>

Timer::Timer() {
	QueryPerformanceFrequency(&m_i64CPUFreq);
}

Timer::~Timer() {}

void Timer::start() {
	reset();
	QueryPerformanceCounter(&m_i64Begin);
}

void Timer::end() {
	QueryPerformanceCounter(&m_i64End);
}

double Timer::getTime() const {
	return (double)(m_i64End.QuadPart-m_i64Begin.QuadPart)/(double)m_i64CPUFreq.QuadPart;
}

void Timer::reset() {
	m_i64Begin = m_i64End;
}