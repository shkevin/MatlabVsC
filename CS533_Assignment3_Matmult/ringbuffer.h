#ifndef _RINGBUFFER_H_
#define _RINGBUFFER_H_

typedef struct ringbuffer_cell {
    void *value;
} ringbuffer_cell_t;

typedef struct ringbuffer {
    pthread_mutex_t *m;
    pthread_cond_t *c;
    ringbuffer_cell_t *head, *tail;
    ringbuffer_cell_t *buf, *end;
} ringbuffer_t;

int RingbufferInit(ringbuffer_t *rb, int len, 
		   pthread_mutex_t *m, pthread_cond_t *c);
int RingbufferAppend(ringbuffer_t *rb, void *item);
int RingbufferRemove(ringbuffer_t *rb, void **item);

#endif /* _RINGBUFFER_H_ */
