/*
 *  Copyright (C) 2000-2009, Thomas Maier-Komor
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

#include "log.h"

#if !(defined(__sun) || defined(__linux) || defined(__GLIBC__))
#define NEED_IO_INTERLOCK
#endif

#include <assert.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__linux)
#include <linux/limits.h>
#elif defined(__bsd)
#include <sys/syslimits.h>
#endif
#include <sys/types.h>
#include <unistd.h>

int Verbose = 3, Log = STDERR_FILENO, ErrorOccurred = 0, ErrorsFatal = 0;
extern char *Prefix;
extern size_t PrefixLen;

#ifdef NEED_IO_INTERLOCK
#include <pthread.h>
pthread_mutex_t
	LogMut = PTHREAD_MUTEX_INITIALIZER;
#endif


#ifdef DEBUG
void logdebug(const char *msg, ...)
{
	va_list val;
	char buf[256];
	size_t b = PrefixLen;

	va_start(val,msg);
	(void) memcpy(buf,Prefix,b);
	b += vsnprintf(buf + b,sizeof(buf)-b,msg,val);
	assert(b < sizeof(buf));
#ifdef NEED_IO_INTERLOCK
	if (b <= PIPE_BUF) {
		(void) write(Log,buf,b);
	} else {
		int err;
		err = pthread_mutex_lock(&LogMut);
		assert(err == 0);
		(void) write(Log,buf,b);
		err = pthread_mutex_unlock(&LogMut);
		assert(err == 0);
	}
#else
	(void) write(Log,buf,b);
#endif
	va_end(val);
}
#define debugmsg if (Verbose >= 5) logdebug
#else
#define debugmsg
#endif


void infomsg(const char *msg, ...)
{
	if (Verbose >= 4) {
		va_list val;
		char buf[256], *b = buf + PrefixLen;
		size_t s;

		va_start(val,msg);
		(void) memcpy(buf,Prefix,PrefixLen);
		b += vsnprintf(b,sizeof(buf)-(b-buf),msg,val);
		s = b - buf;
		assert(s < sizeof(buf));
#ifdef NEED_IO_INTERLOCK
		if (s <= PIPE_BUF) {
			(void) write(Log,buf,s);
		} else {
			int err;
			err = pthread_mutex_lock(&LogMut);
			assert(err == 0);
			(void) write(Log,buf,s);
			err = pthread_mutex_unlock(&LogMut);
			assert(err == 0);
		}
#else
		(void) write(Log,buf,s);
#endif
		va_end(val);
	}
}


void warningmsg(const char *msg, ...)
{
	if (Verbose >= 3) {
		va_list val;
		char buf[256], *b = buf + PrefixLen;
		size_t s;

		va_start(val,msg);
		(void) memcpy(buf,Prefix,PrefixLen);
		(void) memcpy(b,"warning: ",9);
		b += 9;
		b += vsnprintf(b,sizeof(buf)-(b-buf),msg,val);
		s = b - buf;
		assert(s < sizeof(buf));
#ifdef NEED_IO_INTERLOCK
		if (s <= PIPE_BUF) {
			(void) write(Log,buf,s);
		} else {
			int err;
			err = pthread_mutex_lock(&LogMut);
			assert(err == 0);
			(void) write(Log,buf,s);
			err = pthread_mutex_unlock(&LogMut);
			assert(err == 0);
		}
#else
		(void) write(Log,buf,s);
#endif
		va_end(val);
	}
}


void errormsg(const char *msg, ...)
{
	ErrorOccurred = 1;
	if (Verbose >= 2) {
		va_list val;
		char buf[256], *b = buf + PrefixLen;
		size_t s;

		va_start(val,msg);
		(void) memcpy(buf,Prefix,PrefixLen);
		(void) memcpy(b,"error: ",7);
		b += 7;
		b += vsnprintf(b,sizeof(buf)-(b-buf),msg,val);
		s = b - buf;
		assert(s < sizeof(buf));
#ifdef NEED_IO_INTERLOCK
		if (s <= PIPE_BUF) {
			(void) write(Log,buf,s);
		} else {
			int err;
			err = pthread_mutex_lock(&LogMut);
			assert(err == 0);
			(void) write(Log,buf,s);
			err = pthread_mutex_unlock(&LogMut);
			assert(err == 0);
		}
#else
		(void) write(Log,buf,s);
#endif
		va_end(val);
	}
	if (ErrorsFatal) {
		close(Log);
		exit(EXIT_FAILURE);
	}
}


void fatal(const char *msg, ...)
{
	if (Verbose >= 1) {
		va_list val;
		char buf[256], *b = buf + PrefixLen;
		size_t s;

		va_start(val,msg);
		(void) memcpy(buf,Prefix,PrefixLen);
		(void) memcpy(b,"fatal: ",7);
		b += 7;
		b += vsnprintf(b,sizeof(buf)-(b-buf),msg,val);
		s = b - buf;
		assert(s < sizeof(buf));
#ifdef NEED_IO_INTERLOCK
		if (s <= PIPE_BUF) {
			(void) write(Log,buf,s);
		} else {
			int err;
			err = pthread_mutex_lock(&LogMut);
			assert(err == 0);
			(void) write(Log,buf,s);
			err = pthread_mutex_unlock(&LogMut);
			assert(err == 0);
		}
#else
		(void) write(Log,buf,s);
#endif
		va_end(val);
	}
	exit(EXIT_FAILURE);
}


