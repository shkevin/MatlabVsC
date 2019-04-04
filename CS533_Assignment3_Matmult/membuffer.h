#ifndef _MEMBUFFER_H_
#define _MEMBUFFER_H_

#include <sys/types.h>

typedef enum memtype {MB_TYPE_DCACHE = 0, MB_TYPE_MCDRAM, MB_TYPE_DRAM} mb_memtype_t;
typedef enum bufstat {MB_STAT_EMPTY = 0, MB_STAT_FULL, MB_STAT_DIRTY} mb_bufstat_t;
typedef enum reqtype {MB_REQ_NONE = 0, MB_REQ_COPYIN, MB_REQ_COPYOUT} mb_reqtype_t;

struct mb_elem {
	mb_memstat_t status;
	void *bufp;
	size_t buflen;

	mb_reqtype_t request;
	void *extbufp;
	void *extarg;
};
typedef struct mb_elem mb_elem_t;

struct membuffer {
	mb_memtype_t type;
	int enabled;

	struct mb_elem *buffers;
	int nbuffers;

	void (*copyin)(void *dest, void *src, void *arg);
	void (*copyout)(void *dest, void *src, void *arg);

	pthread_mutex_t m;
	pthread_cond_t c;

	ringbuffer_t requests;
};
typedef struct membuffer membuffer_t;


int membuffer_init(membuffer_t b, mb_memtype_t type, size_t size, int nbuffers, 
		   void (*copyin)(void *dest, void *src, void *arg),
		   void (*copyout)(void *dest, void *src, void *arg));
int membuffer_start_copyin(membuffer_t b, int bufnum, void *src, void *arg);
int membuffer_start_copyout(membuffer_t b, int bufnum, void *src, void *arg);
void *membuffer_get_block(membuffer_t b, int bufnum);
int membuffer_release_block(membuffer_t b, int bufnum);
int membuffer_dirty_block(membuffer_t b, int bufnum);
int membuffer_deinit(membuffer_t b);

#endif /* _MEMBUFFER_H_ */
