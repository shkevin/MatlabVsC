/* A replicating multiple producer/multiple consumer ring buffer.
 *
 * - A global head pointer points to the main head of the FIFO
 * - A global tail pointer points to the NEXT FREE ITEM in FIFO

 * - head == tail means the FIFO is empty - writing to the tail will fill in 
     the head
 * - tail + 1 == head means that the FIFO is full (yes, we waste one entry)

 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "ringbuffer.h"

int RingbufferAppend(ringbuffer_t *rb, void *item)
{
    ringbuffer_cell_t *next;

    pthread_mutex_lock(rb->m);
    next = rb->tail + 1;
    
    if (next >= rb->end)
        next = rb->buf;

    while (next == rb->head) { 
	pthread_cond_wait(rb->c, rb->m);
    }

    rb->tail->value = item;
    rb->tail = next;

    pthread_cond_broadcast(rb->c);
    pthread_mutex_unlock(rb->m);
    return 0;
}

int RingbufferRemove(ringbuffer_t *rb, void **item)
{
    ringbuffer_cell_t *head;
    ringbuffer_cell_t *curr;

    pthread_mutex_lock(rb->m);

    /* Check that the ring's not empty */
    while (rb->head == rb->tail) {
	pthread_cond_wait(rb->c, rb->m);
    }

    /* get the item we want and advance the per-thread head. Beause
     * this is all per-thread, there are no atomicity problems
     * here. */
    *item = rb->head->value;
    curr = rb->head + 1;
    if (curr >= rb->end)
	rb->head = rb->buf;
    else
	rb->head = curr;

    pthread_cond_broadcast(rb->c);
    pthread_mutex_unlock(rb->m);
    return 0;
}

int RingbufferDestroy(ringbuffer_t *rb)
{
    if (rb->buf) free(rb->buf);
    return 0;
}

/* Initialize a ringbuffer, giving the length in entries */
int RingbufferInit(ringbuffer_t *rb, int len,
		   pthread_mutex_t *m, pthread_cond_t *c)
{
    rb->m = m;
    rb->c = c;
    rb->head = rb->tail = rb->buf 
	= malloc(len * sizeof(ringbuffer_cell_t));
    rb->end = rb->buf+len;
    return 0;
}
