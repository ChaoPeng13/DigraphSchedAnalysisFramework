/* \file Timer.h
*  this file implements a timer class to count the runtime. 
*  \author Chao Peng
*  
*  Changes
*  ------
*  14-Sept-2015 : initial revision (CP)
*
*/

#ifndef TIMER_H_
#define TIMER_H_

#include <windows.h>

/* \brief A Timer class implements microsecond-level timer
 * From http://blog.csdn.net/hoya5121/article/details/3778487
 */
class Timer {
public:
	Timer();
	~Timer();

	void start();
	void end();
	// return time in second
	double getTime() const;

private:
	LARGE_INTEGER m_i64CPUFreq;
	LARGE_INTEGER m_i64Begin;
	LARGE_INTEGER m_i64End;
	void reset();
};

#endif