#include <sys/time.h>

class SimpleTimer
{
	struct timeval starttime;
public:
	void start();
	double end(); 
};

void SimpleTimer::start()
{
	gettimeofday(&starttime,0);
}

double SimpleTimer::end()
{
	struct timeval endtime;
	gettimeofday(&endtime,0);

	return (endtime.tv_sec - starttime.tv_sec)*1000.0 + (endtime.tv_usec - starttime.tv_usec)/1000.0;
}
