#include "membuffer.h"
#include <stdlib.h>
#include <sys/types.h>

int membuffer_init(membuffer_t b, mb_memtype_t type, size_t size, int nbuffers, 
		   void (*cin)(void *dest, void *src, void *arg),
		   void (*cout)(void *dest, void *src, void *arg))
{
	int i;

	b->type = type;
	b->nbuffers = nbuffers;
	b->buffers = malloc(sizeof(struct mb_elem) * nbuffers);
	b->copyin = cin;
	b->copyout = cout;
	pthread_mutex_init(&b->m, NULL);
	pthread_cond_init(&b->c, NULL);
	
	RingBufferInit(&b->requests, nbuffers + 1, &b->m, &b->c);
	for (i = 0; i < nbuffers; i++) {
		b->buffers[i].status = MB_STAT_EMPTY;
		/* XXX This malloc should take into account the type of memory
 		 * we're buffering in. */
		b->buffers[i].bufp = aligned_alloc(32, size); 
		b->buffers[i].buflen = size;

		b->buffers[i].request = MB_REQ_NONE;
		b->buffers[i].extbufp = NULL;
		b->buffers[i].extbuflen = 0;
	}

	/* XXX Start our worker threads */
	b->enabled = 1;
	return 0;
}

int membuffer_start_copyin(membuffer_t b, int bufnum, void *src, void *arg)
{
	int ret = 0;
	pthread_mutex_lock(&b->m);
	/* Right now we're purely synchronous and single-threaded. We'll
         * need a threading/synchronization approach at some point. */
	if (b->buffers[bufnum].status != MB_STAT_EMPTY) {
		ret = -1;
		goto out;
	}

	b->buffers[bufnum].status = MB_STAT_BUSY;
	b->buffers[bufnum].request = MB_REQ_COPYIN;
	b->buffers[bufnum].extbufp = src;
	b->buffers[bufnum].extarg = arg;
	
	pthread_mutex_unlock(&b->m);
	RingBufferAppend(&b->requests, bufnum);
	pthread_mutex_lock(&b->m);

    out:
	pthread_mutex_unlock(&b->m);
	return ret;
}

void * membuffer_get_block(membuffer_t b, int bufnum) {
	pthread_mutex_lock(&b->m);

	while (b->buffers[bufnum].status != MB_STAT_FULL)
		pthread_cond_wait(&b->c, &b->m);

	pthread_mutex_unlock(&b->m);

	return b->buffers[bufnum].bufp;
}

void * membuffer_worker_thread(pthread_t t, void *arg)
{
	membuffer_t *b = (membuffer_t *)arg;
	int req;
	do {
		void *out;
		int ret;
		ret = (int)(long)RingBufferRemove(&b->requests, &out);
		if (!ret) continue;

		req = (int)(long)out;
		if (req >= 0) {
			pthread_mutex_lock(&b->m);
			/* Process the request here */
			switch (b->buffers[req].request) {
			    case MB_REQ_COPYIN:
				break;
			    case MB_REQ_COPYOUT:
				break;
			    case MB_REQ_NONE:
			    default:
				break;
			} 
			pthread_mutex_unlock(&b->m);
		}
	} while (req != -1);
	return;
}
int membuffer_dirty_block(membuffer_t b, int bufnum)
{
	/* XXX */
	b->bufstatus[bufnum] = MB_STAT_DIRTY;
	return 0;
}

int membuffer_release_block(membuffer_t b, int bufnum)
{
	/* XXX */
	b->bufstatus[bufnum] = MB_STAT_EMPTY;
	return 0;
}

int membuffer_deinit(membuffer_t b)
{
	int i;
	/* XXX */
	for (i = 0; i < b->nbuffers; i++) {
		free(b->bufdata[i]);
	}
	return 0;
}
