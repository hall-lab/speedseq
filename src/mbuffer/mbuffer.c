/*
 *  Copyright (C) 2000-2014, Thomas Maier-Komor
 *
 *  This is the source code of mbuffer.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#ifdef S_SPLINT_S
typedef int caddr_t;
#include <sys/_types.h>
#include <cygwin/types.h>
#include <cygwin/in.h>
#endif

#define _GNU_SOURCE 1	/* needed for O_DIRECT */
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <float.h>
#include <libgen.h>
#include <limits.h>
#include <math.h>
#include <netdb.h>
#include <pthread.h>
#include <semaphore.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <termios.h>
#include <unistd.h>


#ifdef __FreeBSD__
#include <sys/sysctl.h>
#endif

#ifdef HAVE_SENDFILE
#ifdef HAVE_SENDFILE_H
#include <sys/sendfile.h>
#endif
/* if this sendfile implementation does not support sending from buffers,
   disable sendfile support */
	#ifndef SFV_FD_SELF
	#ifdef __GNUC__
	#warning sendfile is unable to send from buffers
	#endif
	#undef HAVE_SENDFILE
	#endif
#endif

#ifndef EBADRQC
#define EBADRQC EINVAL
#endif

#ifdef HAVE_LIBMHASH
#include <mhash.h>
#define HAVE_MD5 1
#elif defined HAVE_LIBMD5
#include <md5.h>
static MD5_CTX MD5ctxt;
#define MD5_INIT(ctxt)		MD5Init(&ctxt);
#define MD5_UPDATE(ctxt,at,num) MD5Update(&ctxt,(unsigned char *)(at),(unsigned int)(num))
#define MD5_END(hash,ctxt)	MD5Final(hash,&(ctxt))
#define HAVE_MD5 1
#elif defined HAVE_LIBCRYPTO
#include <openssl/md5.h>
static MD5_CTX MD5ctxt;
#define MD5_INIT(ctxt)		MD5_Init(&ctxt);
#define MD5_UPDATE(ctxt,at,num)	MD5_Update(&ctxt,at,num)
#define MD5_END(hash,ctxt)	MD5_Final(hash,&(ctxt))
#define HAVE_MD5 1
#endif

/*
 * _POSIX_THREADS is only defined if full thread support is available.
 * We don't need full thread support, so we skip this test...
 * #ifndef _POSIX_THREADS
 * #error posix threads are required
 * #endif
 */

#ifndef S_SPLINT_S
#ifndef _POSIX_SEMAPHORES
#error posix sempahores are required
#endif
#endif

#ifdef O_LARGEFILE
#define LARGEFILE O_LARGEFILE
#else
#define LARGEFILE 0
#endif

#ifdef O_DIRECT
#define DIRECT O_DIRECT
#else
#define DIRECT 0
#endif


#include "dest.h"
#include "network.h"
#include "log.h"

char
	*Prefix;
int 
	In = -1;
size_t
	PrefixLen = 0;

static pthread_t
	Reader, Watchdog;
static long
	Tmp = -1,
	OptSync = 0;
static unsigned long
	Outsize = 10240, Pause = 0, Timeout = 0;
static volatile int
	Terminate = 0,	/* abort execution, because of error or signal */
	EmptyCount = 0,	/* counter incremented when buffer runs empty */
	FullCount = 0,	/* counter incremented when buffer gets full */
	NumSenders = -1,/* number of sender threads */
	Done = 0,
	MainOutOK = 1;	/* is the main outputThread still writing or just coordinating senders */
static unsigned long long
	Totalmem = 0, PgSz = 0, NumP = 0, Blocksize = 10240,
	MaxReadSpeed = 0, MaxWriteSpeed = 0, OutVolsize = 0;
static volatile unsigned long long
	Rest = 0, Numin = 0, Numout = 0;
static double
	StartWrite = 0, StartRead = 1;
static char
	*Tmpfile = 0, **Buffer;
static const char
	*Infile = 0, *AutoloadCmd = 0;
static unsigned int
	AutoloadTime = 0;
static int
	Memlock = 0, TermQ[2],
	Memmap = 0, Quiet = 0, Status = 1, StatusLog = 1,
	Hashers = 0, Direct = 0, SetOutsize = 0;
static long
	NumVolumes = 1,		/* number of input volumes, 0 for interactive prompting */
	Finish = -1,		/* this is for graceful termination */
	Numblocks = 512;	/* number of buffer blocks */
static long long
	TickTime = 0;

static clockid_t
	ClockSrc = CLOCK_REALTIME;

#ifdef __sun
#include <synch.h>
#define sem_t sema_t
#define sem_init(a,b,c) sema_init(a,c,USYNC_THREAD,0)
#define sem_post sema_post
#define sem_getvalue(a,b) ((*(b) = (a)->count), 0)
#if defined(__SunOS_5_8) || defined(__SunOS_5_9)
#define sem_wait SemWait
int SemWait(sema_t *s)
{
	int err;
	do {
		err = sema_wait(s);
	} while (err == EINTR);
	return err;
}
#else
#define sem_wait sema_wait
#endif
#endif

static sem_t Dev2Buf, Buf2Dev;
static pthread_cond_t
	PercLow = PTHREAD_COND_INITIALIZER,	/* low watermark */
	PercHigh = PTHREAD_COND_INITIALIZER,	/* high watermark */
	SendCond = PTHREAD_COND_INITIALIZER;
static pthread_mutex_t
	TermMut = PTHREAD_MUTEX_INITIALIZER,	/* prevents statusThread from interfering with request*Volume */
	LowMut = PTHREAD_MUTEX_INITIALIZER,
	HighMut = PTHREAD_MUTEX_INITIALIZER,
	SendMut = PTHREAD_MUTEX_INITIALIZER;
static int Terminal = 1, Autoloader = 0;
static struct timeval Starttime;
static dest_t *Dest = 0;
static char *volatile SendAt = 0;
static volatile int SendSize = 0, ActSenders = 0;

#ifdef __CYGWIN__
#include <malloc.h>
#undef assert
#define assert(x) ((x) || (*(char *) 0 = 1))
#endif



static int kb2str(char *s, double v)
{
	const char *dim = "KMGT", *f;

	while (v > 10000.0) {
		v /= 1024.0;
		++dim;
		if (*dim == 0) {
			v *= 1024.0*1024.0*1024.0*1024.0;
			break;
		}
	}
	if (v < 0)
		f = " ??? ";
	else if (v < 100)
		f = "%4.1f %ci";
	else if (v < 10000) {
		v = rint(v);
		f = "%4.0f %ci";
	} else
		f = "%5.lg ";
	return sprintf(s,f,v,*dim);
}



static void summary(unsigned long long numb, int numthreads)
{
	int h,m;
	double secs,av;
	char buf[256], *msg = buf;
	struct timeval now;
	
	(void) gettimeofday(&now,0);
	if (Status)
		*msg++ = '\n';
	if ((Terminate == 1) && (numthreads == 0))
		numthreads = 1;
	secs = now.tv_sec - Starttime.tv_sec + (double) now.tv_usec / 1000000 - (double) Starttime.tv_usec / 1000000;
	assert(secs > 0);
	numb >>= 10;
	av = (double)(numb)/secs*numthreads;
	h = (int) secs/3600;
	m = (int) (secs - h * 3600)/60;
	secs -= m * 60 + h * 3600;
	if (numthreads > 1)
		msg += sprintf(msg,"summary: %dx ",numthreads);
	else
		msg += sprintf(msg,"summary: ");
	msg += kb2str(msg,numb);
	msg += sprintf(msg,"Byte in ");
	if (h > 0)
		msg += sprintf(msg,"%dh %02dmin %04.1fsec - average of ",h,m,secs);
	else if (m > 0)
		msg += sprintf(msg,"%2dmin %04.1fsec - average of ",m,secs);
	else
		msg += sprintf(msg,"%4.1fsec - average of ",secs);
	msg += kb2str(msg,av);
	msg += sprintf(msg,"B/s");
	if (EmptyCount != 0)
		msg += sprintf(msg,", %dx empty",EmptyCount);
	if (FullCount != 0)
		msg += sprintf(msg,", %dx full",FullCount);
	*msg++ = '\n';
	*msg = '\0';
	if (Log != STDERR_FILENO)
		(void) write(Log,buf,msg-buf);
	if (Status)
		(void) write(STDERR_FILENO,buf,msg-buf);
}



static void cancelAll(void)
{
	dest_t *d = Dest;
	do {
		(void) pthread_cancel(d->thread);
		if (d->result == 0)
			d->result = "canceled";
		d = d->next;
	} while (d);
	if (Status)
		(void) pthread_cancel(Reader);
}



static RETSIGTYPE sigHandler(int signr)
{
	switch (signr) {
	case SIGHUP:
	case SIGINT:
		ErrorOccurred = 1;
		Terminate = 1;
		(void) close(In);
		if (TermQ[1] != -1) {
			(void) write(TermQ[1],"0",1);
		}
		if (StartWrite > 0)
			(void) pthread_cond_signal(&PercHigh);
		if (StartRead < 1)
			(void) pthread_cond_signal(&PercLow);
		break;
	default:
		(void) raise(SIGABRT);
	}
}



/* Thread-safe replacement for usleep. Argument must be a whole
 * number of microseconds to sleep.
 */
static int mt_usleep(unsigned long sleep_usecs)
{
	struct timespec tv;
	tv.tv_sec = sleep_usecs / 1000000;
	tv.tv_nsec = (sleep_usecs % 1000000) * 1000;

	do {
		/* Sleep for the time specified in tv. If interrupted by a
		 * signal, place the remaining time left to sleep back into tv.
		 */
		if (0 == nanosleep(&tv, &tv)) 
			return 0;
	} while (errno == EINTR);
	return -1;
}



static void *watchdogThread(void *ignored)
{
	unsigned long ni = Numin, no = Numout;
	for (;;) {
		sleep(Timeout);
		if ((ni == Numin) && (Finish == -1)) {
			errormsg("watchdog timeout: input stalled; sending SIGINT\n");
			kill(getpid(),SIGINT);
		}
		if (no == Numout) {
			errormsg("watchdog timeout: output stalled; sending SIGINT\n");
			kill(getpid(),SIGINT);
		}
		ni = Numin;
		no = Numout;
	}
#ifdef __GNUC__
	return 0;	// suppresses a gcc warning
#endif
}


static void statusThread(void) 
{
	struct timeval last, now;
	double in = 0, out = 0, total, diff, fill;
	unsigned long long lin = 0, lout = 0;
	int unwritten = 1;	/* assumption: initially there is at least one unwritten block */
 	fd_set readfds;
 	struct timeval timeout = {0,200000};
	int maxfd = 0;
  
	last = Starttime;
#ifdef __alpha
	(void) mt_usleep(1000);	/* needed on alpha (stderr fails with fpe on nan) */
#endif
	if (TermQ[0] != -1)
		maxfd = TermQ[0]+1;
 	while ((Numin == 0) && (Terminate == 0) && (Finish == -1)) {
 		timeout.tv_sec = 0;
 		timeout.tv_usec = 200000;
		FD_ZERO(&readfds);
		if (TermQ[0] != -1)
			FD_SET(TermQ[0],&readfds);
 		switch (select(maxfd,&readfds,0,0,&timeout)) {
 		case 0: continue;
 		case 1: break;
 		case -1:
 			if (errno == EINTR)
 				break;
 		default: abort();
 		}
 	}
	while (!Done) {
		int err,numsender;
		ssize_t nw = 0;
		char buf[256], *b = buf;

 		timeout.tv_sec = 0;
 		timeout.tv_usec = 500000;
		FD_ZERO(&readfds);
		if (TermQ[0] != -1)
			FD_SET(TermQ[0],&readfds);
 		err = select(maxfd,&readfds,0,0,&timeout);
 		switch (err) {
 		case 0: break;
 		case 1: return;
 		case -1:
 			if (errno == EINTR)
 				break;
 		default: abort();
 		}
		(void) gettimeofday(&now,0);
		diff = now.tv_sec - last.tv_sec + (double) (now.tv_usec - last.tv_usec) / 1000000;
		err = pthread_mutex_lock(&TermMut);
		assert(0 == err);
		err = sem_getvalue(&Buf2Dev,&unwritten);
		assert(0 == err);
		fill = (double)unwritten / (double)Numblocks * 100.0;
		in = (double)(((Numin - lin) * Blocksize) >> 10);
		in /= diff;
		out = (double)(((Numout - lout) * Blocksize) >> 10);
		out /= diff;
		lin = Numin;
		lout = Numout;
		last = now;
		total = (double)((Numout * Blocksize) >> 10);
		fill = (fill < 0.0) ? 0.0 : fill;
		b += sprintf(b,"\rin @ ");
		b += kb2str(b,in);
		numsender = NumSenders + MainOutOK - Hashers;
		b += sprintf(b,"B/s, out @ ");
		b += kb2str(b, out * numsender);
		if (numsender != 1)
			b += sprintf(b,"B/s, %d x ",numsender);
		else
			b += sprintf(b,"B/s, ");
		b += kb2str(b,total);
		b += sprintf(b,"B total, buffer %3.0f%% full",fill);
		if (Quiet == 0) {
#ifdef NEED_IO_INTERLOCK
			if (Log == STDERR_FILENO) {
				int e;
				e = pthread_mutex_lock(&LogMut);
				assert(e == 0);
				nw = write(STDERR_FILENO,buf,strlen(buf));
				e = pthread_mutex_unlock(&LogMut);
				assert(e == 0);
			} else
#endif
				nw = write(STDERR_FILENO,buf,strlen(buf));
		}
		if ((StatusLog != 0) && (Log != STDERR_FILENO))
			infomsg("%s\n",buf+1);
		err = pthread_mutex_unlock(&TermMut);
		assert(0 == err);
		if (nw == -1)	/* stop trying to print status messages after a write error */
			break;
	}
}



static inline long long timediff(struct timespec *restrict t1, struct timespec *restrict t2)
{
	long long tdiff;
	tdiff = (t1->tv_sec - t2->tv_sec) * 1000000;
	tdiff += (t1->tv_nsec - t2->tv_nsec) / 1000;
	if (tdiff < 0)
		tdiff = 0;
	return tdiff;
}



static long long enforceSpeedLimit(unsigned long long limit, long long num, struct timespec *last)
{
	struct timespec now;
	long long tdiff;
	double dt;
	long self = (long) pthread_self();
	
	num += Blocksize;
	if (num < 0) {
		debugmsg("enforceSpeedLimit(%lld,%lld): thread %ld\n",limit,num,self);
		return num;
	}
	(void) clock_gettime(ClockSrc,&now);
	tdiff = timediff(&now,last);
	dt = (double)tdiff * 1E-6;
	if (((double)num/dt) > (double)limit) {
		double req = (double)num/limit - dt;
		long long w = (long long) (req * 1E6);
		if (w >= TickTime) {
			long long slept, ret;
			(void) mt_usleep(w);
			(void) clock_gettime(ClockSrc,last);
			slept = timediff(last,&now);
			ret = -(long long)((double)limit * (double)(slept-w) * 1E-6);
			debugmsg("thread %ld: slept for %lld usec (planned for %lld), ret = %lld\n",self,slept,w,ret);
			return ret;
		} else {
			debugmsg("thread %ld: request for sleeping %lld usec delayed\n",self,w);
			/* 
			 * Sleeping now would cause too much of a slowdown. So
			 * we defer this sleep until the sleeping time is
			 * longer than the tick time. Like this we can stay as
			 * close to the speed limit as possible.
			 */
			return num;
		}
	}
	debugmsg("thread %ld: %lld/%g (%g) <= %g\n",self,num,dt,num/dt,(double)limit);
	return num;
}



static int promptInteractive(unsigned at, unsigned num)
{
	static const char prompt[] = "\nContinue with next volume? Press 'y' to continue or 'n' to finish...";
	static const char contmsg[] = "\nyes - continuing with next volume...\n";
	static const char donemsg[] = "\nno - input done, waiting for output to finish...\n";
	int err;

	err = pthread_mutex_lock(&TermMut);
	assert(0 == err);
	if (-1 == write(STDERR_FILENO,prompt,sizeof(prompt))) {
		errormsg("error accessing controlling terminal for manual volume change request: %s\nConsider using autoload option, when running mbuffer without terminal.\n",strerror(errno));
		Terminate = 1;
		pthread_exit((void *) -1);
	}
	for (;;) {
		char c = 0;
		if (-1 == read(STDERR_FILENO,&c,1) && (errno != EINTR)) {
			errormsg("error accessing controlling terminal for manual volume change request: %s\nConsider using autoload option, when running mbuffer without terminal.\n",strerror(errno));
			Terminate = 1;
			pthread_exit((void *) -1);
		}
		debugmsg("prompt input %c\n",c);
		switch (c) {
		case 'n':
		case 'N':
			Rest = num;
			(void) write(STDERR_FILENO,donemsg,sizeof(donemsg));
			err = pthread_mutex_lock(&HighMut);
			assert(err == 0);
			err = sem_post(&Buf2Dev);
			assert(err == 0);
			err = pthread_cond_signal(&PercHigh);
			assert(err == 0);
			err = pthread_mutex_unlock(&HighMut);
			assert(err == 0);
			err = pthread_mutex_unlock(&TermMut);
			assert(0 == err);
			Finish = at;
			if (Status)
				pthread_exit(0);
			return 0;
		case 'y':
		case 'Y':
			(void) write(STDERR_FILENO,contmsg,sizeof(contmsg));
			err = pthread_mutex_unlock(&TermMut);
			assert(0 == err);
			return 1;
		default:;
		}
	}
}



static int requestInputVolume(unsigned at, unsigned num)
{
	static struct timeval volstart = {0,0};
	const char *cmd;
	struct timeval now;
	double diff;
	unsigned min,hr;
	char cmd_buf[15+strlen(Infile)];

	debugmsg("requesting new volume for input\n");
	(void) gettimeofday(&now,0);
	if (volstart.tv_sec) 
		diff = now.tv_sec - volstart.tv_sec + (double) (now.tv_usec - volstart.tv_usec) * 1E-6;
	else
		diff = now.tv_sec - Starttime.tv_sec + (double) (now.tv_usec - Starttime.tv_usec) * 1E-6;
	if (diff > 3600) {
		hr = (unsigned) (diff / 3600);
		diff -= hr * 3600;
		min = (unsigned) (diff / 60);
		diff -= min * 60;
		infomsg("time for reading volume: %u:%02u:%02f\n",hr,min,diff);
	} else if (diff > 60) {
		min = (unsigned) (diff / 60);
		diff -= min * 60;
		infomsg("time for reading volume: %02u:%02f\n",min,diff);
	} else
		infomsg("time for reading volume: %02fsec.\n",diff);
	if (-1 == close(In))
		errormsg("error closing input: %s\n",strerror(errno));
	do {
		if ((Autoloader) && (Infile)) {
			int ret;
			if (AutoloadCmd) {
				cmd = AutoloadCmd;
			} else {
				(void) snprintf(cmd_buf, sizeof(cmd_buf), "mt -f %s offline", Infile);
				cmd = cmd_buf;
			}
			infomsg("requesting new input volume with command '%s'\n",cmd);
			ret = system(cmd);
			if (0 < ret) {
				warningmsg("error running \"%s\" to change volume in autoloader: exitcode %d\n",cmd,ret);
				Terminate = 1;
				pthread_exit((void *) 0);
			} else if (0 > ret) {
				errormsg("error starting \"%s\" to change volume in autoloader: %s\n", cmd, strerror(errno));
				Terminate = 1;
				pthread_exit((void *) -1);
			}
			if (AutoloadTime) {
				infomsg("waiting for drive to get ready...\n");
				(void) sleep(AutoloadTime);
			}
		} else {
			if (0 == promptInteractive(at,num))
				return 0;
		}
		In = open(Infile, O_RDONLY | LARGEFILE | Direct);
		if ((-1 == In) && (errno == EINVAL))
			In = open(Infile, O_RDONLY | Direct);
		if (-1 == In)
			errormsg("could not reopen input: %s\n",strerror(errno));
#ifdef __sun
		if (-1 == directio(In,DIRECTIO_ON))
			infomsg("direct I/O hinting failed for input: %s\n",strerror(errno));
#endif
	} while (In == -1);
	(void) gettimeofday(&volstart,0);
	diff = volstart.tv_sec - now.tv_sec + (double) (volstart.tv_usec - now.tv_usec) * 1E-6;
	infomsg("tape-change took %fsec. - continuing with next volume\n",diff);
	NumVolumes--;
	if (Terminal && ! Autoloader) {
		char msg[] = "\nOK - continuing...\n";
		(void) write(STDERR_FILENO,msg,sizeof(msg));
	}
	return 1;
}



static void releaseLock(void *l)
{
	int err = pthread_mutex_unlock((pthread_mutex_t *)l);
	assert(err == 0);
}



static void *inputThread(void *ignored)
{
	int fill = 0;
	unsigned long long num;
	int at = 0;
	long long xfer = 0;
	const double startread = StartRead, startwrite = StartWrite;
	struct timespec last;
#ifndef __sun
	int maxfd = TermQ[0] > In ? TermQ[0] + 1 : In + 1;

	if (Status != 0)
		assert(TermQ[0] != -1);
#endif
	(void) clock_gettime(ClockSrc,&last);
	assert(ignored == 0);
	infomsg("inputThread: starting with threadid %ld...\n",(long)pthread_self());
	for (;;) {
		int err;

		if (startread < 1) {
			err = pthread_mutex_lock(&LowMut);
			assert(err == 0);
			err = sem_getvalue(&Buf2Dev,&fill);
			assert(err == 0);
			if (fill == Numblocks - 1) {
				debugmsg("inputThread: buffer full, waiting for it to drain.\n");
				pthread_cleanup_push(releaseLock,&LowMut);
				err = pthread_cond_wait(&PercLow,&LowMut);
				assert(err == 0);
				pthread_cleanup_pop(0);
				++FullCount;
				debugmsg("inputThread: low watermark reached, continuing...\n");
			}
			err = pthread_mutex_unlock(&LowMut);
			assert(err == 0);
		}
		if (Terminate) {	/* for async termination requests */
			debugmsg("inputThread: terminating early upon request...\n");
			if (-1 == close(In))
				errormsg("error closing input: %s\n",strerror(errno));
			if (Status)
				pthread_exit((void *)1);
			return (void *) 1;
		}
		err = sem_wait(&Dev2Buf); /* Wait for one or more buffer blocks to be free */
		assert(err == 0);
		num = 0;
		do {
			int in;
#ifndef __sun
			if (Status != 0) {
				fd_set readfds;
				FD_ZERO(&readfds);
				FD_SET(TermQ[0],&readfds);
				FD_SET(In,&readfds);
				err = select(maxfd,&readfds,0,0,0);
				debugiomsg("inputThread: select(%d, {%d,%d}, 0, 0, 0) = %d\n", maxfd,In,TermQ[0],err);
				assert((err > 0) || (errno == EBADF));
				if (FD_ISSET(TermQ[0],&readfds))
					return (void *)-1;
				assert(FD_ISSET(In,&readfds));
			}
#endif
			in = read(In,Buffer[at] + num,Blocksize - num);
			debugiomsg("inputThread: read(In, Buffer[%d] + %llu, %llu) = %d\n", at, num, Blocksize - num, in);
			if (in > 0) {
				num += in;
			} else if ((0 == in) && (Terminal||Autoloader) && (NumVolumes != 1)) {
				if (0 == requestInputVolume(at,num))
					return 0;
			} else if ((-1 == in) && (errno == EIO) && (Terminal||Autoloader) && (NumVolumes != 1)) {
				requestInputVolume(at,num);
			} else if (in <= 0) {
				/* error or end-of-file */
				if ((-1 == in) && (Terminate == 0))
					errormsg("inputThread: error reading at offset 0x%llx: %s\n",Numin*Blocksize,strerror(errno));
				Rest = num;
				Finish = at;
				debugmsg("inputThread: last block has %llu bytes\n",num);
				err = pthread_mutex_lock(&HighMut);
				assert(err == 0);
				err = sem_post(&Buf2Dev);
				assert(err == 0);
				err = pthread_cond_signal(&PercHigh);
				assert(err == 0);
				err = pthread_mutex_unlock(&HighMut);
				assert(err == 0);
				infomsg("inputThread: exiting...\n");
				if (Status)
					pthread_exit((void *) in);
				return (void *) in;
			}
		} while (num < Blocksize);
		if (MaxReadSpeed)
			xfer = enforceSpeedLimit(MaxReadSpeed,xfer,&last);
		err = sem_post(&Buf2Dev);
		assert(err == 0);
		if (startwrite > 0) {
			err = pthread_mutex_lock(&HighMut);
			assert(err == 0);
			err = sem_getvalue(&Buf2Dev,&fill);
			assert(err == 0);
			if (((double) fill / (double) Numblocks) + DBL_EPSILON >= startwrite) {
				err = pthread_cond_signal(&PercHigh);
				assert(err == 0);
			}
			err = pthread_mutex_unlock(&HighMut);
			assert(err == 0);
		}
		if (++at == Numblocks)
			at = 0;
		Numin++;
	}
}



static inline int syncSenders(char *b, int s)
{
	static volatile int size = 0, skipped = 0;
	static char *volatile buf = 0;
	int err;

	err = pthread_mutex_lock(&SendMut);
	assert(err == 0);
	if (b) {
		buf = b;
		size = s;
	}
	if (s < 0)
		--NumSenders;
	if (--ActSenders) {
		debugiomsg("syncSenders(%p,%d): ActSenders = %d\n",b,s,ActSenders);
		pthread_cleanup_push(releaseLock,&SendMut);
		err = pthread_cond_wait(&SendCond,&SendMut);
		assert(err == 0);
		pthread_cleanup_pop(1);
		debugiomsg("syncSenders(): continue\n");
		return 0;
	} else {
		ActSenders = NumSenders + 1;
		assert((buf != 0) || Terminate);
		SendAt = buf;
		SendSize = size;
		buf = 0;
		if (skipped) {
			// after the first time, always give a buffer free after sync
			err = sem_post(&Dev2Buf);
			assert(err == 0);
		} else {
			// the first time no buffer has been given free
			skipped = 1;
		}
		err = pthread_mutex_unlock(&SendMut);
		assert(err == 0);
		debugiomsg("syncSenders(): send %d@%p, BROADCAST\n",SendSize,SendAt);
		err = pthread_cond_broadcast(&SendCond);
		assert(err == 0);
		return 1;
	}
}



static inline void terminateSender(int fd, dest_t *d, int ret)
{
	debugmsg("terminating operation on %s\n",d->arg);
	if (-1 != fd) {
		int err;
		infomsg("syncing %s...\n",d->arg);
		do 
			err = fsync(fd);
		while ((err != 0) && (errno == EINTR));
		if (err != 0) {
			if ((errno == EINVAL) || (errno == EBADRQC)) {
				infomsg("syncing unsupported on %s: omitted.\n",d->arg);
			} else {
				warningmsg("unable to sync %s: %s\n",d->arg,strerror(errno));
			}
		}
		if (-1 == close(fd))
			errormsg("error closing file %s: %s\n",d->arg,strerror(errno));
	}
	if (ret != 0) {
		ret = syncSenders(0,-1);
		debugmsg("terminateSender(%s): sendSender(0,-1) = %d\n",d->arg,ret);
	}
	pthread_exit((void *) ret);
}



static void *senderThread(void *arg)
{
	unsigned long long outsize = Blocksize;
	dest_t *dest = (dest_t *)arg;
	int out = dest->fd;
#ifdef HAVE_SENDFILE
	int sendout = 1;
#endif
#ifdef HAVE_STRUCT_STAT_ST_BLKSIZE
	struct stat st;

	debugmsg("sender(%s): checking output device...\n",dest->arg);
	if (-1 == fstat(out,&st))
		warningmsg("could not stat output %s: %s\n",dest->arg,strerror(errno));
	else if (S_ISBLK(st.st_mode) || S_ISCHR(st.st_mode)) {
		infomsg("blocksize is %d bytes on output device\n",st.st_blksize);
		if ((Blocksize < st.st_blksize) || (Blocksize % st.st_blksize != 0)) {
			warningmsg("Blocksize should be a multiple of the blocksize of the output device!\n"
				"This can cause problems with some device/OS combinations...\n"
				"Blocksize on output device %s is %d (transfer block size is %lld)\n", dest->arg, st.st_blksize, Blocksize);
			if (SetOutsize) {
				errormsg("unable to set output blocksize\n");
				dest->result = strerror(errno);
				terminateSender(out,dest,1);
			}
		} else {
			if (SetOutsize) {
				infomsg("setting output blocksize to %d\n",st.st_blksize);
				outsize = st.st_blksize;
			}
		}
	} else
		infomsg("no device on output stream %s\n",dest->arg);
#endif
	debugmsg("sender(%s): starting...\n",dest->arg);
	for (;;) {
		int size, num = 0;
		(void) syncSenders(0,0);
		size = SendSize;
		if (0 == size) {
			debugmsg("senderThread(\"%s\"): done.\n",dest->arg);
			terminateSender(out,dest,0);
			return 0;	/* for lint */
		}
		if (Terminate) {
			infomsg("senderThread(\"%s\"): terminating early upon request...\n",dest->arg);
			dest->result = "canceled";
			terminateSender(out,dest,1);
		}
		do {
			unsigned long long rest = size - num;
			int ret;
			assert(size >= num);
#ifdef HAVE_SENDFILE
			if (sendout) {
				off_t baddr = (off_t) (SendAt+num);
				unsigned long long n = SetOutsize ? (rest > Outsize ? (rest/Outsize)*Outsize : rest) : rest;
				ret = sendfile(out,SFV_FD_SELF,&baddr,n);
				debugiomsg("sender(%s): sendfile(%d, SFV_FD_SELF, &%p, %llu) = %d\n", dest->arg, dest->fd, (void*)baddr, n, ret);
				if ((ret == -1) && ((errno == EINVAL) || (errno == EOPNOTSUPP))) {
					sendout = 0;
					debugmsg("sender(%s): sendfile unsupported - falling back to write\n", dest->arg);
					continue;
				}
			} else
#endif
			{
				char *baddr = SendAt+num;
				ret = write(out,baddr,rest > outsize ? outsize :rest);
				debugiomsg("sender(%s): writing %llu@0x%p: ret = %d\n",dest->arg,rest,(void*)baddr,ret);
			}
			if (-1 == ret) {
				errormsg("error writing to %s: %s\n",dest->arg,strerror(errno));
				dest->result = strerror(errno);
				terminateSender(out,dest,1);
			}
			num += ret;
		} while (num != size);
	}
}



static void *hashThread(void *arg)
{
#ifdef HAVE_MD5
	dest_t *dest = (dest_t *) arg;
#ifdef HAVE_LIBMHASH
	int algo = dest->fd;

	MHASH ctxt = mhash_init(algo);
	assert(ctxt != MHASH_FAILED);
#else
	MD5_INIT(MD5ctxt);
#endif
	debugmsg("hashThread(): starting...\n");
	for (;;) {
		int size;

		(void) syncSenders(0,0);
		size = SendSize;
		if (0 == size) {
			size_t ds;
			unsigned char hashvalue[128];
			char *msg, *m;
			const char *an;
			int i;
			
			msg = malloc(300);
			m = msg;
			debugmsg("hashThread(): done.\n");
#ifdef HAVE_LIBMHASH
			mhash_deinit(ctxt,hashvalue);
			an = (const char *) mhash_get_hash_name_static(algo);
			ds = mhash_get_block_size(algo);
#else
			MD5_END(hashvalue,MD5ctxt);
			an = "md5";
			ds = 16;
#endif
			assert(sizeof(hashvalue) >= ds);
			m += sprintf(m,"%s hash: ",an);
			for (i = 0; i < ds; ++i)
				m += sprintf(m,"%02x",(unsigned int)hashvalue[i]);
			*m++ = '\n';
			*m = 0;
			dest->result = msg;
			pthread_exit((void *) msg);
			return 0;	/* for lint */
		}
		if (Terminate) {
			(void) syncSenders(0,-1);
			infomsg("hashThread(): terminating early upon request...\n");
			pthread_exit((void *) 0);
		}
		debugiomsg("hashThread(): hashing %d@0x%p\n",size,(void*)SendAt);
#ifdef HAVE_LIBMHASH
		mhash(ctxt,SendAt,size);
#else
		MD5_UPDATE(MD5ctxt,SendAt,size);
#endif
	}
#endif
}


static int requestOutputVolume(int out, const char *outfile)
{
	static struct timeval volstart = {0,0};
	struct timeval now;
	double diff;
	unsigned min,hr;

	if (!outfile) {
		errormsg("End of volume, but not end of input:\n"
			"Output file must be given (option -o) for multi volume support!\n");
		return -1;
	}
	infomsg("end of volume - last block on volume: %lld\n",Numout);
	(void) gettimeofday(&now,0);
	if (volstart.tv_sec) 
		diff = now.tv_sec - volstart.tv_sec + (double) (now.tv_usec - volstart.tv_usec) * 1E-6;
	else
		diff = now.tv_sec - Starttime.tv_sec + (double) (now.tv_usec - Starttime.tv_usec) * 1E-6;
	if (diff > 3600) {
		hr = (unsigned) (diff / 3600);
		diff -= hr * 3600;
		min = (unsigned) (diff / 60);
		diff -= min * 60;
		infomsg("time for writing volume: %u:%02u:%02f\n",hr,min,diff);
	} else if (diff > 60) {
		min = (unsigned) (diff / 60);
		diff -= min * 60;
		infomsg("time for writing volume: %02u:%02f\n",min,diff);
	} else
		infomsg("time for writing volume: %02fsec.\n",diff);
	if (-1 == close(out))
		errormsg("error closing output %s: %s\n",outfile,strerror(errno));
	do {
		mode_t mode;
		if (Autoloader) {
			const char default_cmd[] = "mt -f %s offline";
			char cmd_buf[sizeof(default_cmd)+strlen(outfile)];
			const char *cmd = AutoloadCmd;
			int err;

			if (cmd == 0) {
				(void) snprintf(cmd_buf, sizeof(cmd_buf), default_cmd, Infile);
				cmd = cmd_buf;
			}
			infomsg("requesting new output volume with command '%s'\n",cmd);
			err = system(cmd);
			if (0 < err) {
				errormsg("error running \"%s\" to change volume in autoloader - exitcode %d\n", cmd, err);
				Autoloader = 0;
				return -1;
			} else if (0 > err) {
				errormsg("error starting \"%s\" to change volume in autoloader: %s\n", cmd, strerror(errno));
				Autoloader = 0;
				return -1;
			}
			if (AutoloadTime) {
				infomsg("waiting for drive to get ready...\n");
				(void) sleep(AutoloadTime);
			}
		} else {
			int err;
			char c = 0, msg[] = "\nvolume full - insert new media and press return when ready...\n";
			if (Terminal == 0) {
				errormsg("End of volume, but not end of input.\n"
					"Specify an autoload command, if you are working without terminal.\n");
				return -1;
			}
			err = pthread_mutex_lock(&TermMut);
			assert(0 == err);
			if (-1 == write(STDERR_FILENO,msg,sizeof(msg))) {
				errormsg("error accessing controlling terminal for manual volume change request: %s\nConsider using autoload option, when running mbuffer without terminal.\n",strerror(errno));
				return -1;
			}
			do {
				if (-1 == read(STDERR_FILENO,&c,1) && (errno != EINTR)) {
					errormsg("error accessing controlling terminal for manual volume change request: %s\nConsider using autoload option, when running mbuffer without terminal.\n",strerror(errno));
					return -1;
				}
			} while (c != '\n');
			err = pthread_mutex_unlock(&TermMut);
			assert(0 == err);
		}
		mode = O_WRONLY|O_TRUNC|OptSync|LARGEFILE|Direct;
		if (strncmp(outfile,"/dev/",5))
			mode |= O_CREAT;
		out = open(outfile,mode,0666);
		if (-1 == out)
			errormsg("error reopening output file: %s\n",strerror(errno));
#ifdef __sun
		if (-1 == directio(out,DIRECTIO_ON))
			infomsg("direct I/O hinting failed for output: %s\n",strerror(errno));
#endif
	} while (-1 == out);
	(void) gettimeofday(&volstart,0);
	diff = volstart.tv_sec - now.tv_sec + (double) (volstart.tv_usec - now.tv_usec) * 1E-6;
	infomsg("tape-change took %fsec. - continuing with next volume\n",diff);
	if (Terminal && ! Autoloader) {
		char msg[] = "\nOK - continuing...\n";
		(void) write(STDERR_FILENO,msg,sizeof(msg));
	}
	return out;
}



static void terminateOutputThread(dest_t *d, int status)
{
	int err;

	infomsg("outputThread: syncing %s...\n",d->arg);
	do 
		err = fsync(d->fd);
	while ((err != 0) && (errno == EINTR));
	if (err != 0) {
		if ((errno == EINVAL) || (errno == EBADRQC)) {
			infomsg("syncing unsupported on %s: omitted.\n",d->arg);
		} else {
			warningmsg("unable to sync %s: %s\n",d->arg,strerror(errno));
		}
	}
	infomsg("outputThread: finished - exiting...\n");
	if (-1 == close(d->fd))
		errormsg("error closing %s: %s\n",d->arg,strerror(errno));
	if (TermQ[1] != -1) {
		err = write(TermQ[1],"0",1);
		if (err == -1)
			errormsg("error writing to termination queue: %s\n",strerror(errno));
	}
	if (status) {
		(void) sem_post(&Dev2Buf);
		(void) pthread_cond_broadcast(&SendCond);
	}
	Done = 1;
	pthread_exit((void *)status);
}



static void *outputThread(void *arg)
{
	dest_t *dest = (dest_t *) arg;
	unsigned at = 0;
	int fill = 0, haderror = 0, out, multipleSenders;
#ifdef HAVE_SENDFILE
	int sendout = 1;
#endif
	const double startwrite = StartWrite, startread = StartRead;
	unsigned long long blocksize = Blocksize;
	long long xfer = 0;
	struct timespec last;

	assert(NumSenders >= 0);
	if (dest->next) {
		int ret;
		dest_t *d = dest->next;
		debugmsg("NumSenders = %d\n",NumSenders);
		ActSenders = NumSenders + 1;
		ret = pthread_mutex_init(&SendMut,0);
		assert(ret == 0);
		ret = pthread_cond_init(&SendCond,0);
		assert(ret == 0);
		do {
			if (d->arg == 0) {
				debugmsg("creating hash thread with algorithm %s\n",d->name);
				ret = pthread_create(&d->thread,0,hashThread,d);
				assert(ret == 0);
			} else if (d->fd != -1) {
				debugmsg("creating sender for %s\n",d->arg);
				ret = pthread_create(&d->thread,0,senderThread,d);
				assert(ret == 0);
			} else {
				debugmsg("outputThread: ignoring destination %s\n",d->arg);
				d->name = 0;
			}
			d = d->next;
		} while (d);
	}
	multipleSenders = (NumSenders > 0);
	dest->result = 0;
	out = dest->fd;
	if (startwrite > 0) {
		int err;
		err = pthread_mutex_lock(&HighMut);
		assert(err == 0);
		debugmsg("outputThread: delaying start until buffer reaches high watermark\n");
		pthread_cleanup_push(releaseLock,&HighMut);
		err = pthread_cond_wait(&PercHigh,&HighMut);
		assert(err == 0);
		pthread_cleanup_pop(0);
		debugmsg("outputThread: high watermark reached, starting...\n");
		err = pthread_mutex_unlock(&HighMut);
		assert(err == 0);
	} else
		infomsg("outputThread: starting output on %s...\n",dest->arg);
	/* initialize last to 0, because we don't want to wait initially */
	(void) clock_gettime(ClockSrc,&last);
	for (;;) {
		unsigned long long rest = blocksize;
		int err;

		if ((startwrite > 0) && (fill <= 0)) {
			assert(fill == 0);
			err = pthread_mutex_lock(&HighMut);
			assert(err == 0);
			err = sem_getvalue(&Buf2Dev,&fill);
			assert(err == 0);
			if (fill == 0) {
				debugmsg("outputThread: buffer empty, waiting for it to fill\n");
				pthread_cleanup_push(releaseLock,&HighMut);
				err = pthread_cond_wait(&PercHigh,&HighMut);
				assert(err == 0);
				pthread_cleanup_pop(0);
				++EmptyCount;
				debugmsg("outputThread: high watermark reached, continuing...\n");
				(void) clock_gettime(ClockSrc,&last);
			}
			err = pthread_mutex_unlock(&HighMut);
			assert(err == 0);
		} else
			--fill;
		err = sem_wait(&Buf2Dev);
		assert(err == 0);
		if (Terminate) {
			infomsg("outputThread: terminating upon termination request...\n");
			dest->result = "canceled";
			terminateOutputThread(dest,1);
		}
		if (Finish == at) {
			err = sem_getvalue(&Buf2Dev,&fill);
			assert(err == 0);
			if ((fill == 0) && (0 == Rest)) {
				if (multipleSenders)
					(void) syncSenders((char*)0xdeadbeef,0);
				infomsg("outputThread: finished - exiting...\n");
				terminateOutputThread(dest,haderror);
			} else {
				blocksize = rest = Rest;
				debugmsg("outputThread: last block has %llu bytes\n",(unsigned long long)Rest);
			}
		}
		if (multipleSenders)
			(void) syncSenders(Buffer[at],blocksize);
		/* switch output volume if -D <size> has been reached */
		if ( (OutVolsize != 0) && (Numout > 0) && (Numout % (OutVolsize/Blocksize)) == 0 ) {
			/* Sleep to let status thread "catch up" so that the displayed total is a multiple of OutVolsize */
			(void) mt_usleep(500000);
			out = requestOutputVolume(out,dest->name);
			if (out == -1) {
				haderror = 1;
				dest->result = strerror(errno);
			}
		}
		do {
			/* use Outsize which could be the blocksize of the device (option -d) */
			unsigned long long n = rest > Outsize ? Outsize : rest;
			int num;
			if (haderror) {
				if (NumSenders == 0)
					Terminate = 1;
				num = (int)rest;
			} else
#ifdef HAVE_SENDFILE
			if (sendout) {
				off_t baddr = (off_t) (Buffer[at] + blocksize - rest);
				num = sendfile(out,SFV_FD_SELF,&baddr,n);
				debugiomsg("outputThread: sendfile(%d, SFV_FD_SELF, &(Buffer[%d] + %llu), %llu) = %d\n", out, at, blocksize - rest, n, num);
				if ((num == -1) && ((errno == EOPNOTSUPP) || (errno == EINVAL))) {
					infomsg("sendfile not supported - falling back to write...\n");
					sendout = 0;
					continue;
				}
			} else
#endif
			{
				num = write(out,Buffer[at] + blocksize - rest, n);
				debugiomsg("outputThread: writing %lld@0x%p: ret = %d\n", n, Buffer[at] + blocksize - rest, num);
			}
			if (Terminal||Autoloader) {
				if (((-1 == num) && ((errno == ENOMEM) || (errno == ENOSPC)))
					|| (0 == num)) {
					/* request a new volume */
					out = requestOutputVolume(out,dest->name);
					if (out == -1)
						haderror = 1;
					continue;
				}
			} else if (-1 == num) {
				dest->result = strerror(errno);
				errormsg("outputThread: error writing to %s at offset 0x%llx: %s\n",dest->arg,(long long)Blocksize*Numout+blocksize-rest,strerror(errno));
				MainOutOK = 0;
				if (NumSenders == 0) {
					debugmsg("outputThread: terminating...\n");
					Terminate = 1;
					err = sem_post(&Dev2Buf);
					assert(err == 0);
					terminateOutputThread(dest,1);
				}
				debugmsg("outputThread: %d senders remaining - continuing...\n",NumSenders);
				haderror = 1;
			}
			rest -= num;
		} while (rest > 0);
		if (multipleSenders == 0) {
			err = sem_post(&Dev2Buf);
			assert(err == 0);
		}
		if (MaxWriteSpeed)
			xfer = enforceSpeedLimit(MaxWriteSpeed,xfer,&last);
		if (Pause)
			(void) mt_usleep(Pause);
		if (Finish == at) {
			err = sem_getvalue(&Buf2Dev,&fill);
			assert(err == 0);
			if (fill == 0) {
				if (multipleSenders)
					(void) syncSenders((char*)0xdeadbeef,0);
				terminateOutputThread(dest,0);
				return 0;	/* make lint happy */
			}
		}
		if (Numblocks == ++at)
			at = 0;
		if (startread < 1) {
			err = pthread_mutex_lock(&LowMut);
			assert(err == 0);
			err = sem_getvalue(&Buf2Dev,&fill);
			assert(err == 0);
			if (((double)fill / (double)Numblocks) < startread) {
				err = pthread_cond_signal(&PercLow);
				assert(err == 0);
			}
			err = pthread_mutex_unlock(&LowMut);
			assert(err == 0);
		}
		Numout++;
	}
}



static void version(void)
{
	(void) fprintf(stderr,
		"mbuffer version "VERSION"\n"\
		"Copyright 2001-2014 - T. Maier-Komor\n"\
		"License: GPLv3 - see file LICENSE\n"\
		"This program comes with ABSOLUTELY NO WARRANTY!!!\n"
		"Donations via PayPal to thomas@maier-komor.de are welcome and support this work!\n"
		"\n"
		);
	exit(EXIT_SUCCESS);
}



static void usage(void)
{
	const char *dim = "bkMGTP";
	unsigned long long m = Numblocks * Blocksize;
	while (m >= 10000) {
		m >>= 10;
		++dim;
	}
	(void) fprintf(stderr,
		"usage: mbuffer [Options]\n"
		"Options:\n"
		"-b <num>   : use <num> blocks for buffer (default: %ld)\n"
		"-s <size>  : use blocks of <size> bytes for processing (default: %llu)\n"
#if defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGESIZE) && !defined(__CYGWIN__) || defined(__FreeBSD__)
		"-m <size>  : memory <size> of buffer in b,k,M,G,%% (default: 2%% = %llu%c)\n"
#else
		"-m <size>  : memory <size> of buffer in b,k,M,G,%% (default: %llu%c)\n"
#endif
#ifdef _POSIX_MEMLOCK_RANGE
		"-L         : lock buffer in memory (unusable with file based buffers)\n"
#endif
		"-d         : use blocksize of device for output\n"
		"-D <size>  : assumed output device size (default: infinite/auto-detect)\n"
		"-P <num>   : start writing after buffer has been filled more than <num>%%\n"
		"-p <num>   : start reading after buffer has been filled less than <num>%%\n"
		"-i <file>  : use <file> for input\n"
		"-o <file>  : use <file> for output (this option can be passed MULTIPLE times)\n"
		"--append   : append to output file (must be passed before -o)\n"
		"--truncate : truncate next file (must be passed before -o)\n"
		"-I <h>:<p> : use network port <port> as input, allow only host <h> to connect\n"
		"-I <p>     : use network port <port> as input\n"
		"-O <h>:<p> : output data to host <h> and port <p> (MUTLIPLE outputs supported)\n"
		"-n <num>   : <num> volumes for input, '0' to prompt interactively\n"
		"-t         : use memory mapped temporary file (for huge buffer)\n"
		"-T <file>  : as -t but uses <file> as buffer\n"
		"-l <file>  : use <file> for logging messages\n"
		"-u <num>   : pause <num> milliseconds after each write\n"
		"-r <rate>  : limit read rate to <rate> B/s, where <rate> can be given in b,k,M,G\n"
		"-R <rate>  : same as -r for writing; use eiter one, if your tape is too fast\n"
		"-f         : overwrite existing files\n"
		"-a <time>  : autoloader which needs <time> seconds to reload\n"
		"-A <cmd>   : issue command <cmd> to request new volume\n"
		"-v <level> : set verbose level to <level> (valid values are 0..6)\n"
		"-q         : quiet - do not display the status on stderr\n"
		"-Q         : quiet - do not log the status\n"
		"-c         : write with synchronous data integrity support\n"
		"-e         : stop processing on any kind of error\n"
#ifdef O_DIRECT
		"--direct   : open input and output with O_DIRECT\n"
#endif
#if defined HAVE_LIBCRYPTO || defined HAVE_LIBMD5 || defined HAVE_LIBMHASH
		"-H\n"
		"--md5      : generate md5 hash of transfered data\n"
		"--hash <a> : use alogritm <a>, if <a> is 'list' possible algorithms are listed\n"
#endif
		"-4         : force use of IPv4\n"
		"-6         : force use of IPv6\n"
		"-0         : use IPv4 or IPv6\n"
		"--tcpbuffer: size for TCP buffer\n"
		"-V\n"
		"--version  : print version information\n"
		"Unsupported buffer options: -t -Z -B\n"
		,Numblocks
		,Blocksize
		,m
		,*dim
		);
	exit(EXIT_SUCCESS);
}


static unsigned long long calcint(const char **argv, int c, unsigned long long def)
{
	char ch;
	double d = (double)def;
	
	switch (sscanf(argv[c],"%lf%c",&d,&ch)) {
	default:
		assert(0);
		break;
	case 2:
		if (d <= 0)
			fatal("invalid argument - must be > 0\n");
		switch (ch) {
		case 'k':
		case 'K':
			d *= 1024.0;
			return (unsigned long long) d;
		case 'm':
		case 'M':
			d *= 1024.0*1024.0;
			return (unsigned long long) d;
		case 'g':
		case 'G':
			d *= 1024.0*1024.0*1024.0;
			return (unsigned long long) d;
		case 't':
		case 'T':
			d *= 1024.0*1024.0*1024.0*1024.0;
			return (unsigned long long) d;
		case '%':
			if ((d >= 90) || (d <= 0))
				fatal("invalid value for percentage (must be 0..90)\n");
			return (unsigned long long) d;
		case 'b':
		case 'B':
			if (d < 128)
				fatal("invalid value for number of bytes\n");
			return (unsigned long long) d;
		default:
			if (argv[c][-2] == '-')
				fatal("unrecognized size charakter \"%c\" for option \"%s\"\n",ch,&argv[c][-2]);
			else
				fatal("unrecognized size charakter \"%c\" for option \"%s\"\n",ch,argv[c-1]);
			return d;
		}
	case 1:
		if (d <= 0)
			fatal("invalid argument - must be > 0\n");
		if (d <= 100) {
			if (argv[c][-2] == '-')
				fatal("invalid low value for option \"%s\" - missing suffix?\n",&argv[c][-2]);
			else
				fatal("invalid low value for option \"%s\" - missing suffix?\n",argv[c-1]);
		}
		return d;
	case 0:
		break;
	}
	errormsg("unrecognized argument \"%s\" for option \"%s\"\n",argv[c],argv[c-1]);
	return d;
}



static int argcheck(const char *opt, const char **argv, int *c, int argc)
{
	if (strncmp(opt,argv[*c],strlen(opt))) 
		return 1;
	if (strlen(argv[*c]) > 2)
		argv[*c] += 2;
	else {
		(*c)++;
		if (*c == argc)
			fatal("missing argument to option %s\n",opt);
	}
	return 0;
}



static void addHashAlgorithm(const char *name)
{
	const char *algoname = "";
	int algo = 0;
#if HAVE_LIBMHASH
	int numalgo = mhash_count();

	while (algo <= numalgo) {
		algoname = (const char *) mhash_get_hash_name_static(algo);
		if (algoname && (strcasecmp(algoname,name) == 0))
			break;
		++algo;
	}
#else
	algoname = "MD5";
#endif
	if (strcasecmp(algoname,name) == 0) {
		dest_t *dest = malloc(sizeof(dest_t));
		bzero(dest,sizeof(dest_t));
		dest->name = algoname;
		dest->fd = algo;
		if (Dest) {
			dest->next = Dest->next;
			Dest->next = dest;
		} else {
			Dest = dest;
			dest->next = 0;
		}
		debugmsg("enabled hash algorithm %s\n",name);
		++NumSenders;
		++Hashers;
	} else
		fatal("invalid or unsupported hash function %s\n",name);
}


static void openDestinationFiles(dest_t *d)
{
	unsigned errs = ErrorOccurred;
	while (d) {
		if (d->fd == -1) {
			if (0 == strncmp(d->arg,"/dev/",5))
				d->mode &= ~O_EXCL;
			d->fd = open(d->arg,d->mode,0666);
			if ((-1 == d->fd) && (errno == EINVAL)) {
				d->mode &= ~LARGEFILE;
				d->fd = open(d->arg,d->mode,0666);
			}
			if ((-1 == d->fd) && (errno == EINVAL)) {
				d->mode &= ~O_TRUNC;
				d->fd = open(d->arg,d->mode,0666);
			}
			if (-1 == d->fd) {
				d->result = strerror(errno);
				errormsg("unable to open output %s: %s\n",d->arg,strerror(errno));
			} else {
				debugmsg("successfully opened destination file %s with fd %d\n",d->arg,d->fd);
			}
		}
		if (-1 == d->fd) {
			d->name = 0;	/* tag destination as unstartable */
			--NumSenders;
		}
#ifdef __sun
		else if (d->arg) {
			if (0 == directio(d->fd,DIRECTIO_ON))
				infomsg("direct I/O hinting enabled for output to %s\n",d->arg);
			else
				infomsg("direct I/O hinting failed for output to %s: %s\n",d->arg,strerror(errno));
		}
#endif
		d = d->next;
	}
	if (ErrorOccurred != errs)
		fatal("unable to open all outputs\n");
}



static const char *calcval(const char *arg, unsigned long long *res)
{
	char ch;
	double d;
	
	switch (sscanf(arg,"%lf%c",&d,&ch)) {
	default:
		assert(0);
		break;
	case 2:
		if (d <= 0)
			return "negative value out of range";
		switch (ch) {
		case 'k':
		case 'K':
			d *= 1024.0;
			*res = d;
			return 0;
		case 'm':
		case 'M':
			d *= 1024.0*1024.0;
			*res = d;
			return 0;
		case 'g':
		case 'G':
			d *= 1024.0*1024.0*1024.0;
			*res = d;
			return 0;
		case 't':
		case 'T':
			d *= 1024.0*1024.0*1024.0*1024.0;
			*res = d;
			return 0;
		case '%':
			if ((d >= 90) || (d <= 0))
				return "invalid value for percentage (must be 0..90)";
			*res = d;
			return 0;
		case 'b':
		case 'B':
			if (d < 128)
				return "invalid value for number of bytes";
			*res = d;
			return 0;
		default:
			return "invalid dimension";
		}
	case 1:
		if (d <= 0)
			return "value out of range";
		if (d <= 100)
			return "value out of range";
		*res = d;
		return 0;
	case 0:
		break;
	}
	return "unrecognized argument";
}


static void initDefaults()
{
	char dfname[PATH_MAX+1], line[256];
	const char *home = getenv("HOME");
	size_t l;
	int df;
	FILE *dfstr;
	struct stat st;

	if (home == 0) {
		warningmsg("HOME environment variable not set - unable to find defaults file\n");
		return;
	}
	strncpy(dfname,home,PATH_MAX);
	dfname[sizeof(dfname)-1] = 0;
	l = strlen(dfname);
	if (l + 12 > PATH_MAX) {
		warningmsg("path to defaults file breaks PATH_MAX\n");
		return;
	}
	strcat(dfname,"/.mbuffer.rc");
	df = open(dfname,O_RDONLY);
	if (df == -1) {
		if (errno == ENOENT)
			infomsg("no defaults file ~/.mbuffer.rc\n");
		else
			warningmsg("error opening defaults file %s: %s\n",dfname,strerror(errno));
		return;
	}
	if (-1 == fstat(df,&st)) {
		warningmsg("unable to stat defaults file %s: %s\n",dfname,strerror(errno));
		close(df);
		return;
	}
	if (getuid() != st.st_uid) {
		warningmsg("ignoring defaults file from different user\n");
		close(df);
		return;
	}
	infomsg("reading defaults file %s\n",dfname);
	dfstr = fdopen(df,"r");
	assert(dfstr);
	while (!feof(dfstr)) {
		char key[64],valuestr[64];
		fscanf(dfstr,"%255[^\n]\n",line);
		char *pound = strchr(line,'#');
		unsigned long long value;
		int a;

		if (pound)
			*pound = 0;
		a = sscanf(line,"%63[A-Za-z]%*[ \t=:]%63[0-9a-zA-Z]",key,valuestr);
		if (a != 2) {
			warningmsg("unable to parse line '%s' in .mbuffer.rc; %d arguments\n",line,a);
			continue;
		}
		debugmsg("parsing key/value pair %s=%s\n",key,valuestr);
		if (strcasecmp(key,"numblocks") == 0) {
			long nb = strtol(valuestr,0,0);
			if ((nb == 0) && (errno == EINVAL)) {
				warningmsg("invalid argument for %s: \"%s\"\n",key,valuestr);
			} else {
				Numblocks = nb;
				debugmsg("Numblocks = %llu\n",Numblocks);
			}
		} else if (strcasecmp(key,"pause") == 0) {
			long p = strtol(valuestr,0,0);
			if ((p == 0) && (errno == EINVAL)) {
				warningmsg("invalid argument for %s: \"%s\"\n",key,valuestr);
			} else {
				Pause = p;
				debugmsg("Pause = %d\n",Pause);
			}
		} else if (strcasecmp(key,"autoloadtime") == 0) {
			long at = strtol(valuestr,0,0) - 1;
			if ((at == 0) && (errno == EINVAL)) 
				warningmsg("invalid argument for %s: \"%s\"\n",key,valuestr);
			else {
				AutoloadTime = at;
				debugmsg("Autoloader time = %d\n",AutoloadTime);
			}
		} else if (strcasecmp(key,"startread") == 0) {
			double sr = 0;
			if (1 == sscanf(valuestr,"%lf",&sr))
				sr /= 100;
			if ((sr <= 1) && (sr > 0)) {
				StartRead = sr;
				debugmsg("StartRead = %1.2lf\n",StartRead);
			}
		} else if (strcasecmp(key,"startwrite") == 0) {
			double sw = 0;
			if (1 == sscanf(valuestr,"%lf",&sw))
				sw /= 100;
			if ((sw <= 1) && (sw > 0)) {
				StartWrite = sw;
				debugmsg("StartWrite = %1.2lf\n",StartWrite);
			}
		} else if (strcasecmp(key,"timeout") == 0) {
			long t = strtol(valuestr,0,0);
			if (((t == 0) && (errno == EINVAL)) || (t < 0)) 
				warningmsg("invalid argument for %s: \"%s\"\n",key,valuestr);
			else {
				Timeout = t;
				debugmsg("Timeout = %lu\n",Timeout);
			}
		} else if (strcasecmp(key,"showstatus") == 0) {
			if ((strcasecmp(valuestr,"yes") == 0) || (strcasecmp(valuestr,"on") == 0) || (strcmp(valuestr,"1") == 0)) {
				Quiet = 0;
				debugmsg("showstatus = yes\n");
			} else if ((strcasecmp(valuestr,"no") == 0) || (strcasecmp(valuestr,"off") == 0) || (strcmp(valuestr,"0") == 0)) {
				Quiet = 1;
				debugmsg("showstatus = no\n");
			} else 
				warningmsg("invalid argument for %s: \"%s\"\n",key,valuestr);
			continue;
		} else if (strcasecmp(key,"logstatus") == 0) {
			if ((strcasecmp(valuestr,"yes") == 0) || (strcasecmp(valuestr,"on") == 0) || (strcmp(valuestr,"1") == 0)) {
				StatusLog = 1;
				debugmsg("logstatus = yes\n");
			} else if ((strcasecmp(valuestr,"no") == 0) || (strcasecmp(valuestr,"off") == 0) || (strcmp(valuestr,"0") == 0)) {
				StatusLog = 0;
				debugmsg("logstatus = no\n");
			} else 
				warningmsg("invalid argument for %s: \"%s\"\n",key,valuestr);
			continue;
		} else if (strcasecmp(key,"memlock") == 0) {
			if ((strcasecmp(valuestr,"yes") == 0) || (strcasecmp(valuestr,"on") == 0) || (strcmp(valuestr,"1") == 0)) {
				Memlock = 1;
				debugmsg("Memlock = %lu\n",Memlock);
			} else if ((strcasecmp(valuestr,"no") == 0) || (strcasecmp(valuestr,"off") == 0) || (strcmp(valuestr,"0") == 0)) {
				Memlock = 0;
				debugmsg("Memlock = %lu\n",Memlock);
			} else 
				warningmsg("invalid argument for %s: \"%s\"\n",key,valuestr);
			continue;
		}
		const char *argerror = calcval(valuestr,&value);
		if (argerror) {
			warningmsg("ignoring key/value pair from defaults file (%s = %s): %s\n",key,valuestr,argerror);
			continue;
		}
		if (strcasecmp(key,"blocksize") == 0) {
			Blocksize = value;
		} else if (strcasecmp(key,"maxwritespeed") == 0) {
			MaxWriteSpeed = value;
		} else if (strcasecmp(key,"maxreadspeed") == 0) {
			MaxReadSpeed = value;
		} else if (strcasecmp(key,"Totalmem") == 0) {
			if (value < 100) {
#if defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGESIZE) && !defined(__CYGWIN__) || defined(__FreeBSD__)
				Totalmem = ((unsigned long long) NumP * PgSz * value) / 100 ;
				debugmsg("Totalmem = %lluk\n",Totalmem>>10);
#else
				warningmsg("Unable to determine page size or amount of available memory - please specify an absolute amount of memory.\n");
#endif
			}
		} else if (strcasecmp(key,"tcpbuffer") == 0) {
			TCPBufSize = value;
		} else {
			warningmsg("unknown key: %s\n",key);
			continue;
		}
		infomsg("setting %s to %lld\n",key,value);
	}
	fclose(dfstr);
	close(df);
}


int main(int argc, const char **argv)
{
	int optMset = 0, optSset = 0, optBset = 0, optMode = O_EXCL, numOut = 0;
	int  numstdout = 0, numthreads = 0;
	long mxnrsem;
	int c, fl, err;
	sigset_t       signalSet;
#ifdef HAVE_STRUCT_STAT_ST_BLKSIZE
	struct stat st;
#endif
	unsigned short netPortIn = 0;
	unsigned short netPortOut = 0;
	char *argv0 = strdup(argv[0]), *progname, null;
	const char *outfile = 0;
	struct sigaction sig;
	dest_t *dest = 0;

	/* setup logging prefix */
	progname = basename(argv0);
	PrefixLen = strlen(progname) + 2;
	Prefix = malloc(PrefixLen);
	(void) strcpy(Prefix,progname);
	Prefix[PrefixLen - 2] = ':';
	Prefix[PrefixLen - 1] = ' ';

	/* set verbose level before parsing defaults and options */
	for (c = 1; c < argc; c++) {
		const char *arg = argv[c];
		if ((arg[0] == '-') && (arg[1] == 'v')) {
			long verb;
			if (arg[2])
				verb = strtol(arg+2,0,0);
			else
				verb = strtol(argv[++c],0,0);
			if ((verb == 0) && (errno == EINVAL))
				errormsg("invalid argument to option -v: \"%s\"\n",argv[c]);
			else
				Verbose = verb;
			debugmsg("Verbose = %d\n",Verbose);
		}
	}
	
	/* gather system parameters */
	TickTime = 1000000 / sysconf(_SC_CLK_TCK);
#if defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGESIZE) && !defined(__CYGWIN__)
	PgSz = sysconf(_SC_PAGESIZE);
	assert(PgSz > 0);
	NumP = sysconf(_SC_AVPHYS_PAGES);
	assert(NumP > 0);
	Blocksize = PgSz;
	debugmsg("total # of phys pages: %li (pagesize %li)\n",NumP,PgSz);
	Numblocks = NumP/50;
#elif defined(__FreeBSD__)
	size_t nump_size = sizeof(nump_size);
	sysctlbyname("hw.availpages", &NumP, &nump_size, NULL, 0);
	PgSz = sysconf(_SC_PAGESIZE);
	assert(PgSz > 0);
#endif
#if defined(_POSIX_MONOTONIC_CLOCK) && (_POSIX_MONOTONIC_CLOCK >= 0) && defined(CLOCK_MONOTONIC)
	if (sysconf(_SC_MONOTONIC_CLOCK) > 0)
		ClockSrc = CLOCK_MONOTONIC;
#endif

	/* setup parameters */
	initDefaults();
	debugmsg("default buffer set to %d blocks of %lld bytes\n",Numblocks,Blocksize);
	for (c = 1; c < argc; c++) {
		if (!argcheck("-s",argv,&c,argc)) {
			Blocksize = Outsize = calcint(argv,c,Blocksize);
			optSset = 1;
			debugmsg("Blocksize = %llu\n",Blocksize);
			if (Blocksize < 100)
				fatal("cannot set blocksize as percentage of total physical memory\n");
		} else if (!strcmp("--append",argv[c])) {
			optMode |= O_APPEND;
			debugmsg("append to next file\n");
		} else if (!strcmp("--truncate",argv[c])) {
			optMode |= O_TRUNC;
			debugmsg("truncate next file\n");
		} else if (!argcheck("-m",argv,&c,argc)) {
			Totalmem = calcint(argv,c,Totalmem);
			optMset = 1;
			if (Totalmem < 100) {
#if defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGESIZE) && !defined(__CYGWIN__) || defined(__FreeBSD__)
				Totalmem = ((unsigned long long) NumP * PgSz * Totalmem) / 100 ;
#else
				fatal("Unable to determine page size or amount of available memory - please specify an absolute amount of memory.\n");
#endif
			}
			debugmsg("Totalmem = %lluk\n",Totalmem>>10);
		} else if (!argcheck("-b",argv,&c,argc)) {
			long nb = strtol(argv[c],0,0);
			if ((nb == 0) && (errno == EINVAL)) {
				errormsg("invalid argument to option -b: \"%s\"\n",argv[c]);
			} else {
				Numblocks = nb;
				optBset = 1;
			}
			debugmsg("Numblocks = %llu\n",Numblocks);
		} else if (!strcmp("--tcpbuffer",argv[c])) {
			TCPBufSize = calcint(argv,++c,TCPBufSize);
			debugmsg("TCPBufSize = %lu\n",TCPBufSize);
		} else if (!argcheck("-d",argv,&c,argc)) {
#ifdef HAVE_STRUCT_STAT_ST_BLKSIZE
			SetOutsize = 1;
			debugmsg("setting output size according to the blocksize of the device\n");
#else
			fatal("cannot determine blocksize of device (unsupported by OS)\n");
#endif
		} else if (!argcheck("-v",argv,&c,argc)) {
			/* has been parsed already */
		} else if (!argcheck("-u",argv,&c,argc)) {

			long p = strtol(argv[c],0,0);
			if ((p == 0) && (errno == EINVAL))
				errormsg("invalid argument to option -u: \"%s\"\n",argv[c]);
			else
				Pause = p;
			debugmsg("Pause = %d\n",Pause);
		} else if (!argcheck("-r",argv,&c,argc)) {
			MaxReadSpeed = calcint(argv,c,0);
			debugmsg("MaxReadSpeed = %lld\n",MaxReadSpeed);
		} else if (!argcheck("-R",argv,&c,argc)) {
			MaxWriteSpeed = calcint(argv,c,0);
			debugmsg("MaxWriteSpeed = %lld\n",MaxWriteSpeed);
		} else if (!argcheck("-n",argv,&c,argc)) {
			long nv = strtol(argv[c],0,0);
			if ((nv < 0) || ((nv == 0) && (errno == EINVAL)))
				fatal("invalid argument to option -n: \"%s\"\n",argv[c]);
			else
				NumVolumes = nv;
			if (NumVolumes < 0)
				fatal("argument for number of volumes must be > 0\n");
			debugmsg("NumVolumes = %d\n",NumVolumes);
		} else if (!argcheck("-i",argv,&c,argc)) {
			if (strcmp(argv[c],"-")) {
				Infile = argv[c];
				debugmsg("Infile = %s\n",Infile);
			} else {
				Infile = STDIN_FILENO;
				debugmsg("Infile is stdin\n");
			}
		} else if (!argcheck("-o",argv,&c,argc)) {
			dest_t *dest = malloc(sizeof(dest_t));
			if (strcmp(argv[c],"-")) {
				debugmsg("output file: %s\n",argv[c]);
				dest->arg = argv[c];
				dest->name = argv[c];
				dest->fd = -1;
				dest->mode = O_CREAT|O_WRONLY|optMode|Direct|LARGEFILE|OptSync;
			} else {
				if (numstdout++) 
					fatal("cannot output multiple times to stdout\n");
				debugmsg("output to stdout\n",argv[c]);
				dest->fd = dup(STDOUT_FILENO);
				err = dup2(STDERR_FILENO,STDOUT_FILENO);
				assert(err != -1);
				dest->arg = "<stdout>";
				dest->name = "<stdout>";
				dest->mode = 0;
			}
			optMode = O_EXCL;
			dest->port = 0;
			dest->result = 0;
			bzero(&dest->thread,sizeof(dest->thread));
			dest->next = Dest;
			Dest = dest;
			if (outfile == 0)
				outfile = argv[c];
			++numOut;
			++NumSenders;
#ifdef AF_INET6
		} else if (!strcmp("-0",argv[c])) {
			AddrFam = AF_UNSPEC;
		} else if (!strcmp("-4",argv[c])) {
			AddrFam = AF_INET;
		} else if (!strcmp("-6",argv[c])) {
			AddrFam = AF_INET6;
#endif
		} else if (!argcheck("-I",argv,&c,argc)) {
			initNetworkInput(argv[c]);
		} else if (!argcheck("-O",argv,&c,argc)) {
			dest_t *d = createNetworkOutput(argv[c]);
			if (d->fd == -1) {
				free(d);
			} else {
				d->next = Dest;
				Dest = d;
				++NumSenders;
			}
			++numOut;
		} else if (!argcheck("-T",argv,&c,argc)) {
			Tmpfile = malloc(strlen(argv[c]) + 1);
			if (!Tmpfile)
				fatal("out of memory\n");
			(void) strcpy(Tmpfile, argv[c]);
			Memmap = 1;
			debugmsg("Tmpfile = %s\n",Tmpfile);
		} else if (!strcmp("-t",argv[c])) {
			Memmap = 1;
			debugmsg("Memmap = 1\n");
		} else if (!argcheck("-l",argv,&c,argc)) {
			Log = open(argv[c],O_WRONLY|O_APPEND|O_TRUNC|O_CREAT|LARGEFILE,0666);
			if (-1 == Log) {
				Log = STDERR_FILENO;
				errormsg("error opening log file: %s\n",strerror(errno));
			}
			debugmsg("logFile set to %s\n",argv[c]);
		} else if (!strcmp("-f",argv[c])) {
			optMode &= ~O_EXCL;
			debugmsg("overwrite = 1\n");
		} else if (!strcmp("-q",argv[c])) {
			debugmsg("disabling display of status\n");
			Quiet = 1;
		} else if (!strcmp("-Q",argv[c])) {
			debugmsg("disabling logging of status\n");
			StatusLog = 0;
		} else if (!strcmp("-c",argv[c])) {
			debugmsg("enabling full synchronous I/O\n");
			OptSync = O_SYNC;
		} else if (!strcmp("-e",argv[c])) {
			debugmsg("will terminate on any kind of error\n");
			ErrorsFatal = 1;
		} else if (!argcheck("-a",argv,&c,argc)) {
			long at = strtol(argv[c],0,0) - 1;
			if ((at == 0) && (errno == EINVAL)) 
				errormsg("invalid argument to option -a: \"%s\"\n",argv[c]);
			else {
				Autoloader = 1;
				AutoloadTime = at;
			}
			debugmsg("Autoloader time = %d\n",AutoloadTime);
		} else if (!argcheck("-A",argv,&c,argc)) {
			Autoloader = 1;
			AutoloadCmd = argv[c];
			debugmsg("Autoloader command = \"%s\"\n", AutoloadCmd);
		} else if (!argcheck("-P",argv,&c,argc)) {
			if (1 != sscanf(argv[c],"%lf",&StartWrite))
				StartWrite = 0;
			StartWrite /= 100;
			if ((StartWrite > 1) || (StartWrite <= 0))
				fatal("error in argument -P: must be bigger than 0 and less or equal 100\n");
			debugmsg("StartWrite = %1.2lf\n",StartWrite);
		} else if (!argcheck("-p",argv,&c,argc)) {
			if (1 == sscanf(argv[c],"%lf",&StartRead))
				StartRead /= 100;
			else
				StartRead = 1.0;
			if ((StartRead >= 1) || (StartRead < 0))
				fatal("error in argument -p: must be bigger or equal to 0 and less than 100\n");
			debugmsg("StartRead = %1.2lf\n",StartRead);
		} else if (!strcmp("-L",argv[c])) {
#ifdef _POSIX_MEMLOCK_RANGE
			Memlock = 1;
			debugmsg("memory locking enabled\n");
#else
			warning("POSIX memory locking is unsupported on this system.\n");
#endif
		} else if (!argcheck("-W",argv,&c,argc)) {
			Timeout = strtol(argv[c],0,0);
			if (Timeout <= 0)
				fatal("invalid argument to option -W\n");
			if (Timeout <= AutoloadTime)
				fatal("timeout must be bigger than autoload time\n");
		} else if (!strcmp("--direct",argv[c])) {
#ifdef O_DIRECT
			debugmsg("using O_DIRECT to open file descriptors\n");
			Direct = O_DIRECT;
#else
			warningmsg("--direct is unsupported on this system\n");
#endif
		} else if (!strcmp("--help",argv[c]) || !strcmp("-h",argv[c])) {
			usage();
		} else if (!strcmp("--version",argv[c]) || !strcmp("-V",argv[c])) {
			version();
		} else if (!strcmp("--md5",argv[c]) || !strcmp("-H",argv[c])) {
#ifdef HAVE_MD5
			addHashAlgorithm("MD5");
#else
			fatal("hash calculation support has not been compiled in!\n");
#endif
		} else if (!strcmp("--hash",argv[c])) {
			++c;
			if (c == argc)
				fatal("missing argument to option --hash\n");
#if HAVE_LIBMHASH
			if (!strcmp(argv[c],"list")) {
				(void) fprintf(stderr,"valid hash functions are:\n");
				int algo = mhash_count();
				while (algo >= 0) {
					const char *algoname = (const char *) mhash_get_hash_name_static(algo);
					if (algoname)
						(void) fprintf(stderr,"\t%s\n",algoname);
					--algo;
				}
				exit(EXIT_SUCCESS);
			}
#elif defined HAVE_MD5
			if (!strcmp(argv[c],"list")) {
				(void) fprintf(stderr,"valid hash functions are:\n");
				(void) fprintf(stderr,"\tMD5\n");
				exit(EXIT_SUCCESS);
			}
#else
			fatal("hash calculation support has not been compiled in!\n");
#endif
			addHashAlgorithm(argv[c]);
		} else if (!argcheck("-D",argv,&c,argc)) {
			OutVolsize = calcint(argv,c,0);
			debugmsg("OutVolsize = %llu\n",OutVolsize);
		} else
			fatal("unknown option \"%s\"\n",argv[c]);
	}

	/* consistency check for options */
	if (AutoloadTime && Timeout && Timeout <= AutoloadTime)
		fatal("autoload time must be smaller than watchdog timeout\n");
	if (optBset&optSset&optMset) {
		if (Numblocks * Blocksize != Totalmem)
			fatal("inconsistent options: blocksize * number of blocks != totalsize!\n");
	} else if (((!optBset)&optSset&optMset) || (optMset&(!optBset)&(!optSset))) {
		if (Totalmem <= Blocksize)
			fatal("total memory must be larger than block size\n");
		Numblocks = Totalmem / Blocksize;
		infomsg("Numblocks = %llu, Blocksize = %llu, Totalmem = %llu\n",(unsigned long long)Numblocks,(unsigned long long)Blocksize,(unsigned long long)Totalmem);
	} else if (optBset&!optSset&optMset) {
		if (Blocksize == 0)
			fatal("blocksize must be greater than 0\n");
		if (Totalmem <= Blocksize)
			fatal("total memory must be larger than block size\n");
		Blocksize = Totalmem / Numblocks;
		infomsg("blocksize = %llu\n",(unsigned long long)Blocksize);
	}
	if ((StartRead < 1) && (StartWrite > 0))
		fatal("setting both low watermark and high watermark doesn't make any sense...\n");
	if ((NumSenders-Hashers > 0) && (Autoloader || OutVolsize))
		fatal("multi-volume support is unsupported with multiple outputs\n");
	if (Autoloader) {
		if ((!outfile) && (!Infile))
			fatal("Setting autoloader time without using a device doesn't make any sense!\n");
		if (outfile && Infile) {
			fatal("Which one is your autoloader? Input or output? Replace input or output with a pipe.\n");
		}
	}
	if (Infile && netPortIn)
		fatal("Setting both network input port and input file doesn't make sense!\n");
	if (outfile && netPortOut)
		fatal("Setting both network output and output file doesn't make sense!\n");

	/* multi volume input consistency checking */
	if ((NumVolumes != 1) && (!Infile))
		fatal("multi volume support for input needs an explicit given input device (option -i)\n");

	/* SPW: Volsize consistency checking */
	if (OutVolsize && !outfile)
		fatal("Setting OutVolsize without an output device doesn't make sense!\n");
	if ((OutVolsize != 0) && (OutVolsize < Blocksize))
		/* code assumes we can write at least one block */
		fatal("If non-zero, OutVolsize must be at least as large as the buffer blocksize (%llu)!\n",Blocksize);
	/* SPW END */

	/* check that we stay within system limits */
#ifdef __FreeBSD__
	{
		size_t semvmx_size = sizeof(mxnrsem);
		if (sysctlbyname("kern.ipc.semvmx", &mxnrsem, &semvmx_size, 0, 0) == -1)
			mxnrsem = -1;
	}
#else
	mxnrsem = sysconf(_SC_SEM_VALUE_MAX);
#endif
	if (-1 == mxnrsem) {
#ifdef SEM_MAX_VALUE
		mxnrsem = SEM_MAX_VALUE;
#else
		mxnrsem = LONG_MAX;
		warningmsg("unable to determine maximum value of semaphores\n");
#endif
	}
	if (Numblocks < 5)
		fatal("Minimum block count is 5.\n");
	if (Numblocks > mxnrsem) {
		fatal("cannot allocate more than %d blocks.\nThis is a system dependent limit, depending on the maximum semaphore value.\nPlease choose a bigger block size.\n",mxnrsem);
	}

	if ((Blocksize * (long long)Numblocks) > (long long)SSIZE_MAX)
		fatal("Cannot address so much memory (%lld*%d=%lld>%lld).\n",Blocksize,Numblocks,Blocksize*(long long)Numblocks,(long long)SSIZE_MAX);
	/* create buffer */
	Buffer = (char **) valloc(Numblocks * sizeof(char *));
	if (!Buffer)
		fatal("Could not allocate enough memory (%d requested): %s\n",Numblocks * sizeof(char *),strerror(errno));
	if (Memmap) {
		infomsg("mapping temporary file to memory with %llu blocks with %llu byte (%llu kB total)...\n",(unsigned long long) Numblocks,(unsigned long long) Blocksize,(unsigned long long) ((Numblocks*Blocksize) >> 10));
		if (!Tmpfile) {
			char tmplname[] = "mbuffer-XXXXXX";
			char *tmpdir = getenv("TMPDIR") ? getenv("TMPDIR") : "/var/tmp";
			char tfilename[sizeof(tmplname) + strlen(tmpdir) + 1];
			(void) strcpy(tfilename,tmpdir);
			(void) strcat(tfilename,"/");
			(void) strcat(tfilename,tmplname);
			Tmp = mkstemp(tfilename);
			Tmpfile = malloc(strlen(tfilename));
			if (!Tmpfile)
				fatal("out of memory: %s\n",strerror(errno));
			(void) strcpy(Tmpfile,tfilename);
			infomsg("tmpfile is %s\n",Tmpfile);
		} else {
			mode_t mode = O_RDWR | LARGEFILE;
			if (strncmp(Tmpfile,"/dev/",5))
				mode |= O_CREAT|O_EXCL;
			Tmp = open(Tmpfile,mode,0600);
		}
		if (-1 == Tmp)
			fatal("could not create temporary file (%s): %s\n",Tmpfile,strerror(errno));
		if (strncmp(Tmpfile,"/dev/",5))
			(void) unlink(Tmpfile);
		/* resize the file. Needed - at least under linux, who knows why? */
		if (-1 == lseek(Tmp,Numblocks * Blocksize - sizeof(int),SEEK_SET))
			fatal("could not resize temporary file: %s\n",strerror(errno));
		if (-1 == write(Tmp,&c,sizeof(int)))
			fatal("could not resize temporary file: %s\n",strerror(errno));
		Buffer[0] = mmap(0,Blocksize*Numblocks,PROT_READ|PROT_WRITE,MAP_SHARED,Tmp,0);
		if (MAP_FAILED == Buffer[0])
			fatal("could not map buffer-file to memory: %s\n",strerror(errno));
		debugmsg("temporary file mapped to address %p\n",Buffer[0]);
	} else {
		infomsg("allocating memory for %d blocks with %llu byte (%llu kB total)...\n",Numblocks,(unsigned long long) Blocksize,(unsigned long long) ((Numblocks*Blocksize) >> 10));
		Buffer[0] = (char *) valloc(Blocksize * Numblocks);
		if (Buffer[0] == 0)
			fatal("Could not allocate enough memory (%lld requested): %s\n",(unsigned long long)Blocksize * Numblocks,strerror(errno));
#ifdef MADV_DONTFORK
		if (-1 == madvise(Buffer[0],Blocksize * Numblocks, MADV_DONTFORK))
			warningmsg("unable to advise memory handling of buffer: %s\n",strerror(errno));
#endif
	}
	for (c = 1; c < Numblocks; c++) {
		Buffer[c] = Buffer[0] + Blocksize * c;
		*Buffer[c] = 0;	/* touch every block before locking */
	}

#ifdef _POSIX_MEMLOCK_RANGE
	if (Memlock) {
		uid_t uid;
#ifndef HAVE_SETEUID
#define seteuid setuid
#endif
		uid = geteuid();
		if (0 != seteuid(0))
			warningmsg("could not change to uid 0 to lock memory (is mbuffer setuid root?)\n");
		else if ((0 != mlock((char *)Buffer,Numblocks * sizeof(char *))) || (0 != mlock(Buffer[0],Blocksize * Numblocks)))
			warningmsg("could not lock buffer in memory: %s\n",strerror(errno));
		else
			infomsg("memory locked successfully\n");
		err = seteuid(uid);	/* don't give anyone a chance to attack this program, so giveup uid 0 after locking... */
		assert(err == 0);
	}
#endif

	debugmsg("creating semaphores...\n");
	if (0 != sem_init(&Buf2Dev,0,0))
		fatal("Error creating semaphore Buf2Dev: %s\n",strerror(errno));
	if (0 != sem_init(&Dev2Buf,0,Numblocks))
		fatal("Error creating semaphore Dev2Buf: %s\n",strerror(errno));

	debugmsg("opening input...\n");
	if (Infile) {
		int flags = O_RDONLY | LARGEFILE | Direct;
		In = open(Infile,flags);
		if (-1 == In) {
			if (errno == EINVAL) {
				flags &= ~LARGEFILE;
				In = open(Infile,flags);
			}
			if (-1 == In)
				fatal("could not open input file: %s\n",strerror(errno));
		}
	} else if (In == -1) {
		In = STDIN_FILENO;
	}
#ifdef __sun
	if (0 == directio(In,DIRECTIO_ON))
		infomsg("direct I/O hinting enabled for input\n");
	else
		infomsg("direct I/O hinting failed for input: %s\n",strerror(errno));
#endif
	if (numOut == 0) {
		dest_t *d = malloc(sizeof(dest_t));
		d->fd = dup(STDOUT_FILENO);
		err = dup2(STDERR_FILENO,STDOUT_FILENO);
		assert(err != -1);
		d->name = "<stdout>";
		d->arg = "<stdout>";
		d->port = 0;
		d->result = 0;
		bzero(&d->thread,sizeof(d->thread));
		d->next = Dest;
		Dest = d;
		++NumSenders;
	}
	openDestinationFiles(Dest);
	if (NumSenders == -1) {
		fatal("no output left - nothing to do\n");
	}

	debugmsg("checking if we have a controlling terminal...\n");
	sig.sa_handler = SIG_IGN;
	err = sigaction(SIGTTIN,&sig,0);
	assert(err == 0);
	fl = fcntl(STDERR_FILENO,F_GETFL);
	err = fcntl(STDERR_FILENO,F_SETFL,fl | O_NONBLOCK);
	assert(err == 0);
	if ((read(STDERR_FILENO,&c,1) == -1) && (errno != EAGAIN)) {
		int tty = open("/dev/tty",O_RDWR);
		if (-1 == tty) {
			Terminal = 0;
			if (Autoloader == 0)
				warningmsg("No controlling terminal and no autoloader command specified.\n");
		} else {
			Terminal = 1;
			err = dup2(tty,STDERR_FILENO);
			assert(err != -1);
		}
	}
	err = fcntl(STDERR_FILENO,F_SETFL,fl);
	assert(err == 0);
	if ((Terminal == 1) && (NumVolumes != 1)) {
		struct termios tset;
		if (-1 == tcgetattr(STDERR_FILENO,&tset)) {
			warningmsg("unable to get terminal attributes: %s\n",strerror(errno));
		} else {
			tset.c_lflag &= (~ICANON) & (~ECHO);
			tset.c_cc[VTIME] = 0;
			tset.c_cc[VMIN] = 1;
			if (-1 == tcsetattr(STDERR_FILENO,TCSANOW,&tset))
				warningmsg("unable to set terminal attributes: %s\n",strerror(errno));
		}
	}

	debugmsg("registering signals...\n");
	sig.sa_handler = sigHandler;
	err = sigemptyset(&sig.sa_mask);
	assert(err == 0);
	err = sigaddset(&sig.sa_mask,SIGINT);
	assert(err == 0);
	sig.sa_flags = SA_RESTART;
	if (0 != sigaction(SIGINT,&sig,0))
		warningmsg("error registering new SIGINT handler: %s\n",strerror(errno));
	err = sigemptyset(&sig.sa_mask);
	assert(err == 0);
	err = sigaddset(&sig.sa_mask,SIGHUP);
	assert(err == 0);
	if (0 != sigaction(SIGHUP,&sig,0))
		warningmsg("error registering new SIGHUP handler: %s\n",strerror(errno));

	debugmsg("starting threads...\n");
	(void) gettimeofday(&Starttime,0);
	err = sigfillset(&signalSet);
	assert(0 == err);
	(void) pthread_sigmask(SIG_BLOCK, &signalSet, NULL);

	/* select destination for output thread */
	dest = Dest;
	while (dest->fd == -1) {
		dest->name = 0;
		debugmsg("skipping destination %s\n",dest->arg);
		assert(dest->next);
		dest = dest->next;
	}

#ifdef HAVE_STRUCT_STAT_ST_BLKSIZE
	debugmsg("checking output device...\n");
	if (-1 == fstat(dest->fd,&st))
		errormsg("could not stat output: %s\n",strerror(errno));
	else if (S_ISBLK(st.st_mode) || S_ISCHR(st.st_mode)) {
		infomsg("blocksize is %d bytes on output device\n",st.st_blksize);
		if (Blocksize % st.st_blksize != 0) {
			warningmsg("Blocksize should be a multiple of the blocksize of the output device!\n"
				"This can cause problems with some device/OS combinations...\n"
				"Blocksize on output device is %d (transfer block size is %lld)\n", st.st_blksize, Blocksize);
			if (SetOutsize)
				fatal("unable to set output blocksize\n");
		} else {
			if (SetOutsize) {
				infomsg("setting output blocksize to %d\n",st.st_blksize);
				Outsize = st.st_blksize;
			}
		}
	} else
		infomsg("no device on output stream\n");
	debugmsg("checking input device...\n");
	if (-1 == fstat(In,&st))
		warningmsg("could not stat input: %s\n",strerror(errno));
	else if (S_ISBLK(st.st_mode) || S_ISCHR(st.st_mode)) {
		infomsg("blocksize is %d bytes on input device\n",st.st_blksize);
		if (Blocksize % st.st_blksize != 0) {
			warningmsg("Blocksize should be a multiple of the blocksize of the input device!\n"
				"Use option -s to adjust transfer block size if you get an out-of-memory error on input.\n"
				"Blocksize on input device is %d (transfer block size is %lld)\n", st.st_blksize, Blocksize);
		}
	} else
		infomsg("no device on input stream\n");
#else
	warningmsg("Could not stat output device (unsupported by system)!\n"
		   "This can result in incorrect written data when\n"
		   "using multiple volumes. Continue at your own risk!\n");
#endif
	if (((Verbose < 4) || (StatusLog == 0)) && (Quiet != 0))
		Status = 0;
	if (Status) {
		if (-1 == pipe(TermQ))
			fatal("could not create termination pipe: %s\n",strerror(errno));
	} else {
		TermQ[0] = -1;
		TermQ[1] = -1;
	}
	err = pthread_create(&dest->thread,0,&outputThread,dest);
	assert(0 == err);
	if (Timeout) {
		err = pthread_create(&Watchdog,0,&watchdogThread,(void*)0);
		assert(0 == err);
	}
	if (Status) {
		err = pthread_create(&Reader,0,&inputThread,0);
		assert(0 == err);
		(void) pthread_sigmask(SIG_UNBLOCK, &signalSet, NULL);
		statusThread();
		err = pthread_join(Reader,0);
		if (err != 0)
			errormsg("error joining reader: %s\n",strerror(errno));
	} else {
		(void) pthread_sigmask(SIG_UNBLOCK, &signalSet, NULL);
		(void) inputThread(0);
		debugmsg("waiting for output to finish...\n");
		if (TermQ[0] != -1) {
			err = read(TermQ[0],&null,1);
			assert(err == 1);
		}
	}
	if (Timeout) {
		err = pthread_cancel(Watchdog);
		assert(err == 0);
	}
	if (Dest) {
		dest_t *d = Dest;
		int ret;

		infomsg("waiting for senders...\n");
		if (Terminate)
			cancelAll();
		do {
			if (d->name) {
				void *status;
				if (d->arg) {
					debugmsg("joining sender for %s\n",d->arg);
				} else {
					debugmsg("joining hasher for %s\n",d->name);
				}
				ret = pthread_join(d->thread,&status);
				if (ret != 0)
					errormsg("error joining %s: %s\n",d->arg,d->name,strerror(errno));
				if (status == 0)
					++numthreads;
			}
			d = d->next;
		} while (d);
	}
	if (Status || Log != STDERR_FILENO)
		summary(Numout * Blocksize + Rest, numthreads);
	if (Memmap) {
		int ret = munmap(Buffer[0],Blocksize*Numblocks);
		assert(ret == 0);
	}
	if (Tmp != -1)
		(void) close(Tmp);
	if (Dest) {
		dest_t *d = Dest;
		do {
			dest_t *n = d->next;
			if (d->result) {
				if (d->arg) {
					warningmsg("error during output to %s: %s\n",d->arg,d->result);
				} else {
					(void) write(STDERR_FILENO,d->result,strlen(d->result));
					if (Log != STDERR_FILENO) 
						(void) write(Log,d->result,strlen(d->result));
				}
			}
			free(d);
			d = n;
		} while (d);
	}
	if (ErrorOccurred)
		exit(EXIT_FAILURE);
	exit(EXIT_SUCCESS);
}

/* vim:tw=0
 */
