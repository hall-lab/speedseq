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

#include "config.h"

#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#elif defined __GNUC__
#define alloca __builtin_alloca
#elif defined _AIX
#define alloca __alloca
#else
#include <stddef.h>
void *alloca(size_t);
#endif


#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include "dest.h"
#include "network.h"
#include "log.h"

extern int In;
int32_t TCPBufSize = 1 << 20;
#if defined(PF_INET6) && defined(PF_UNSPEC)
int AddrFam = PF_UNSPEC;
#else
int AddrFam = PF_INET;
#endif


static void setTCPBufferSize(int sock, unsigned buffer)
{
	int err;
	int32_t osize, size;
	socklen_t bsize = sizeof(osize);

	assert(buffer == SO_RCVBUF || buffer == SO_SNDBUF);
	err = getsockopt(sock,SOL_SOCKET,buffer,&osize,&bsize);
	assert((err == 0) && (bsize == sizeof(osize)));
	if (osize < TCPBufSize) {
		size = TCPBufSize;
		do {
			err = setsockopt(sock,SOL_SOCKET,buffer,(void *)&size,sizeof(size));
			size >>= 1;
		} while ((-1 == err) && (errno == ENOMEM) && (size > osize));
		if (err == -1) {
			warningmsg("unable to set socket buffer size: %s\n",strerror(errno));
			return;
		}
	}
	bsize = sizeof(size);
	err = getsockopt(sock,SOL_SOCKET,buffer,&size,&bsize);
	assert(err != -1);
	if (buffer == SO_RCVBUF) 
		infomsg("set TCP receive buffer size to %d\n",size);
	else
		infomsg("set TCP send buffer size to %d\n",size);
}


#ifdef HAVE_GETADDRINFO

void initNetworkInput(const char *addr)
{
	char *host, *port;
	struct addrinfo hint, *pinfo = 0, *x, *cinfo = 0;
	int err, sock = -1, l;

	debugmsg("initNetworkInput(\"%s\")\n",addr);
	l = strlen(addr) + 1;
	host = alloca(l);
	memcpy(host,addr,l);
	port = strrchr(host,':');
	if (port == 0) {
		port = host;
		host = 0;
	} else if (port == host) {
		port = host + 1;
		host = 0;
	} else {
		*port = 0;
		++port;
		bzero(&hint,sizeof(hint));
		hint.ai_family = AddrFam;
		hint.ai_protocol = IPPROTO_TCP;
		hint.ai_socktype = SOCK_STREAM;
#ifdef __FreeBSD__
		hint.ai_flags = AI_ADDRCONFIG;
#else
		hint.ai_flags = AI_ADDRCONFIG | AI_V4MAPPED;
#endif
		err = getaddrinfo(host,0,&hint,&cinfo);
		if (err != 0) 
			fatal("unable to resolve address information for expected host '%s': %s\n",host,gai_strerror(err));
	}
	bzero(&hint,sizeof(hint));
	hint.ai_family = AddrFam;
	hint.ai_protocol = IPPROTO_TCP;
	hint.ai_socktype = SOCK_STREAM;
	hint.ai_flags = AI_PASSIVE | AI_ADDRCONFIG;
	err = getaddrinfo(0,port,&hint,&pinfo);
	if (err != 0)
		fatal("unable to get address information for port/service '%s': %s\n",port,gai_strerror(err));
	assert(pinfo);
	for (x = pinfo; x; x = x->ai_next) {
		int reuse_addr = 1;
		debugmsg("creating socket for address familiy %d\n",x->ai_family);
		sock = socket(x->ai_family, SOCK_STREAM, 0);
		if (sock == -1) {
			warningmsg("unable to create socket for input: %s\n",strerror(errno));
			continue;
		}
		if (-1 == setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &reuse_addr, sizeof(reuse_addr)))
			warningmsg("cannot set socket to reuse address: %s\n",strerror(errno));
		if (0 == bind(sock, x->ai_addr, x->ai_addrlen)) {
			debugmsg("successfully bound socket - address length %d\n",x->ai_addrlen);
			break;
		}
		warningmsg("could not bind to socket for network input: %s\n",strerror(errno));
		(void) close(sock);
	}
	if (x == 0)
		fatal("Unable to initialize network input.\n");
	infomsg("listening on socket...\n");
	if (0 > listen(sock,0))		/* accept only 1 incoming connection */
		fatal("could not listen on socket for network input: %s\n",strerror(errno));
	for (;;) {
		char chost[NI_MAXHOST], serv[NI_MAXSERV];
		struct sockaddr_in6 caddr;
		struct addrinfo *c;
		socklen_t len = sizeof(caddr);
		int err;

		debugmsg("waiting for incoming connection\n");
		In = accept(sock, (struct sockaddr *) &caddr, &len);
		if (0 > In)
			fatal("Unable to accept connection for network input: %s\n",strerror(errno));
		err = getnameinfo((struct sockaddr *) &caddr,len,chost,sizeof(chost),serv,sizeof(serv),NI_NUMERICHOST|NI_NUMERICSERV|NI_NOFQDN);
		if (0 != err) {
			fatal("unable to get name information for hostname of incoming connection: %s\n",gai_strerror(err));
		}
		infomsg("incoming connection from %s:%s\n",chost,serv);
		if (host == 0)
			break;
		for (c = cinfo; c; c = c->ai_next) {
			char xhost[NI_MAXHOST];
			if (0 == getnameinfo((struct sockaddr *)c->ai_addr,c->ai_addrlen,xhost,sizeof(xhost),0,0,NI_NUMERICHOST|NI_NOFQDN)) {
				debugmsg("checking against host '%s'\n",xhost);
				if (0 == strcmp(xhost,chost))
					break;
			}
		}
		if (c)
			break;
		warningmsg("rejected connection from %s\n",chost);
		if (-1 == close(In))
			warningmsg("error closing rejected input: %s\n",strerror(errno));
	}
	freeaddrinfo(pinfo);
	if (cinfo)
		freeaddrinfo(cinfo);
	debugmsg("input connection accepted\n");
	setTCPBufferSize(In,SO_RCVBUF);
	(void) close(sock);
}


dest_t *createNetworkOutput(const char *addr)
{
	char *host, *port;
	struct addrinfo hint, *ret = 0, *x;
	int err, fd = -1;
	dest_t *d;

	assert(addr);
	host = strdup(addr);
	assert(host);
	port = strrchr(host,':');
	if (port == 0) {
		fatal("syntax error - target must be given in the form <host>:<port>\n");
	}
	*port++ = 0;
	bzero(&hint,sizeof(hint));
	hint.ai_family = AddrFam;
	hint.ai_protocol = IPPROTO_TCP;
	hint.ai_socktype = SOCK_STREAM;
	hint.ai_flags = AI_ADDRCONFIG;
	debugmsg("getting address info for %s\n",addr);
	err = getaddrinfo(host,port,&hint,&ret);
	if (err != 0)
		fatal("unable to resolve address information for '%s': %s\n",addr,gai_strerror(err));
	for (x = ret; x; x = x->ai_next) {
		fd = socket(x->ai_family, SOCK_STREAM, 0);
		if (fd == -1) {
			errormsg("unable to create socket: %s\n",strerror(errno));
			continue;
		}
		if (0 == connect(fd, x->ai_addr, x->ai_addrlen)) {
			debugmsg("successfully connected to %s\n",addr);
			break;
		}
		(void) close(fd);
		fd = -1;
		warningmsg("error connecting to %s: %s\n",addr,strerror(errno));
	}
	if ((x == 0) || (fd == -1))
		errormsg("unable to connect to %s\n",addr);
	freeaddrinfo(ret);
	if (fd != -1)
		setTCPBufferSize(fd,SO_SNDBUF);
	d = (dest_t *) malloc(sizeof(dest_t));
	d->arg = addr;
	d->name = host;
	d->port = port;
	d->fd = fd;
	bzero(&d->thread,sizeof(d->thread));
	d->result = 0;
	d->next = 0;
	return d;
}


#else	/* HAVE_GETADDRINFO */


static void openNetworkInput(const char *host, unsigned short port)
{
	struct sockaddr_in saddr;
	struct hostent *h = 0, *r = 0;
	const int reuse_addr = 1;
	int sock;

	debugmsg("openNetworkInput(\"%s\",%hu)\n",host,port);
	sock = socket(AF_INET, SOCK_STREAM, 6);
	if (0 > sock)
		fatal("could not create socket for network input: %s\n",strerror(errno));
	if (-1 == setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &reuse_addr, sizeof(reuse_addr)))
		warningmsg("cannot set socket to reuse address: %s\n",strerror(errno));
	setTCPBufferSize(sock,SO_RCVBUF);
	if (host[0]) {
		debugmsg("resolving hostname '%s' of input...\n",host);
		if (0 == (h = gethostbyname(host)))
#ifdef HAVE_HSTRERROR
			fatal("could not resolve server hostname: %s\n",hstrerror(h_errno));
#else
			fatal("could not resolve server hostname: error code %d\n",h_errno);
#endif
	}
	bzero((void *) &saddr, sizeof(saddr));
	saddr.sin_family = AF_INET;
	saddr.sin_addr.s_addr = htonl(INADDR_ANY);
	saddr.sin_port = htons(port);
	debugmsg("binding socket to port %d...\n",port);
	if (0 > bind(sock, (struct sockaddr *) &saddr, sizeof(saddr)))
		fatal("could not bind to socket for network input: %s\n",strerror(errno));
	debugmsg("listening on socket...\n");
	if (0 > listen(sock,1))		/* accept only 1 incoming connection */
		fatal("could not listen on socket for network input: %s\n",strerror(errno));
	for (;;) {
		struct sockaddr_in caddr;
		socklen_t clen = sizeof(caddr);
		char **p;
		debugmsg("waiting to accept connection...\n");
		In = accept(sock, (struct sockaddr *)&caddr, &clen);
		if (0 > In)
			fatal("could not accept connection for network input: %s\n",strerror(errno));
		if (host[0] == 0) {
			infomsg("accepted connection from %s\n",inet_ntoa(caddr.sin_addr));
			(void) close(sock);
			return;
		}
		for (p = h->h_addr_list; *p; ++p) {
			if (0 == memcmp(&caddr.sin_addr,*p,h->h_length)) {
				infomsg("accepted connection from %s\n",inet_ntoa(caddr.sin_addr));
				(void) close(sock);
				return;
			}
		}
		r = gethostbyaddr((char *)&caddr.sin_addr,sizeof(caddr.sin_addr.s_addr),AF_INET);
		if (r)
			warningmsg("rejected connection from %s (%s)\n",r->h_name,inet_ntoa(caddr.sin_addr));
		else
			warningmsg("rejected connection from %s\n",inet_ntoa(caddr.sin_addr));
		if (-1 == close(In))
			warningmsg("error closing rejected input: %s\n",strerror(errno));
	}
}


void initNetworkInput(const char *addr)
{
	char *host, *portstr;
	unsigned pnr;
	size_t l;

	debugmsg("initNetworkInput(\"%s\")\n",addr);
	l = strlen(addr) + 1;
	host = alloca(l);
	memcpy(host,addr,l);
	portstr = strrchr(host,':');
	if (portstr == 0) {
		portstr = host;
		host = "";
	} else if (portstr == host) {
		portstr = host + 1;
		host = "";
		*portstr = 0;
	} else {
		*portstr = 0;
		++portstr;
	}
	if (1 != sscanf(portstr,"%u",&pnr))
		fatal("invalid port string '%s' - port must be given by its number, not service name\n", portstr);
	openNetworkInput(host,pnr);
}


static void openNetworkOutput(dest_t *dest)
{
	struct sockaddr_in saddr;
	struct hostent *h = 0;
	int out;
	unsigned short pnr;

	debugmsg("creating socket for output to %s:%d...\n",dest->name,dest->port);
	if (1 != sscanf(dest->port,"%hu",&pnr))
		fatal("port must be given by its number, not service name\n");
	out = socket(PF_INET, SOCK_STREAM, 0);
	if (0 > out) {
		errormsg("could not create socket for network output: %s\n",strerror(errno));
		return;
	}
	setTCPBufferSize(out,SO_SNDBUF);
	bzero((void *) &saddr, sizeof(saddr));
	saddr.sin_port = htons(pnr);
	infomsg("resolving host %s...\n",dest->name);
	if (0 == (h = gethostbyname(dest->name))) {
#ifdef HAVE_HSTRERROR
		dest->result = hstrerror(h_errno);
		errormsg("could not resolve hostname %s: %s\n",dest->name,dest->result);
#else
		dest->result = "unable to resolve hostname";
		errormsg("could not resolve hostname %s: error code %d\n",dest->name,h_errno);
#endif
		dest->fd = -1;
		(void) close(out);
		return;
	}
	saddr.sin_family = h->h_addrtype;
	assert(h->h_length <= sizeof(saddr.sin_addr));
	(void) memcpy(&saddr.sin_addr,h->h_addr_list[0],h->h_length);
	infomsg("connecting to server at %s...\n",inet_ntoa(saddr.sin_addr));
	if (0 > connect(out, (struct sockaddr *) &saddr, sizeof(saddr))) {
		dest->result = strerror(errno);
		errormsg("could not connect to %s:%d: %s\n",dest->name,dest->port,dest->result);
		(void) close(out);
		out = -1;
	}
	dest->fd = out;
}


dest_t *createNetworkOutput(const char *addr)
{
	char *host, *portstr;
	dest_t *d = (dest_t *) malloc(sizeof(dest_t));

	debugmsg("createNetworkOutput(\"%s\")\n",addr);
	host = strdup(addr);
	portstr = strrchr(host,':');
	if ((portstr == 0) || (portstr == host))
		fatal("argument '%s' doesn't match <host>:<port> format\n",addr);
	*portstr++ = 0;
	bzero(d, sizeof(dest_t));
	d->fd = -1;
	d->arg = addr;
	d->name = host;
	d->port = portstr;
	openNetworkOutput(d);
	return d;
}


#endif /* HAVE_GETADDRINFO */


/* vim:tw=0
 */
