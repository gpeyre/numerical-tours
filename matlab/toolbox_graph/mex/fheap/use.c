#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "fib.h"

#define TESTCASE	1

#define COUNT	100000
#define DIF	1000
#define MAXEXT	10
#define VERBOSE	1

int cmp(void *, void *);

int
cmp(void *x, void *y)
{
	int a, b;
	a = (int)x;
	b = (int)y;

	if (a < b)
		return -1;
	if (a == b)
		return 0;
	return 1;
}

int
main(void)
{
#if !TESTCASE
	struct fibheap_el *w;
#else
	int e, j, k;
#endif
	struct fibheap *a;
	int i, x;

	a = fh_makeheap();
	fh_setcmp(a, cmp);

	srandom(time(NULL));
#if TESTCASE
#if VERBOSE
	printf("inserting: ");
#endif
	e = 0;
	for (i = 0; i < COUNT; i++) {
#if VERBOSE
		if (i)
			printf(", ");
#endif
		fh_insert(a, (void *)(x = random()/10));
#if VERBOSE
		printf("%d", x);
#endif
		if (i - e > DIF) {
			k = random() % MAXEXT;
			for (j = 0; j < k; j++, e++)
				printf("throwing: %d\n", (int)fh_extractmin(a));
		}
	}

#if VERBOSE
	printf("\nremaining: %d\n", COUNT - e);
	printf("extracting: ");
#endif
	for (i = 0; i < COUNT - e; i++) {
#if VERBOSE
		if (i)
			printf(", ");
		printf("%d", (int)fh_extractmin(a));
#else
		fh_extractmin(a);
#endif
	}
#if VERBOSE
	printf("\n");
#endif
	if ((int)fh_extractmin(a) == 0)
		printf("heap empty!\n");
	else {
		printf("heap not empty! ERROR!\n");
		exit(1);
	}
#else
	w = fh_insert(a, (void *)6);
	printf("adding: %d\n", 6);
	fh_insert(a, (void *)9);
	printf("adding: %d\n", 9);
	fh_insert(a, (void *)1);
	printf("adding: %d\n", 1);
	for(i = 0; i < 5; i++) {
		x = random()/10000;
		printf("adding: %d\n", x);
		fh_insert(a, (void *)x);
	}
	fh_insert(a, (void *)4);
	printf("adding: %d\n", 4);
	fh_insert(a, (void *)8);
	printf("adding: %d\n", 8);
	printf("first: %d\n", (int)fh_extractmin(a));
	printf("deleting: %d\n", (int)fh_delete(a, w));
	printf("first: %d\n", (int)fh_extractmin(a));
	printf("first: %d\n", (int)fh_extractmin(a));
	printf("first: %d\n", (int)fh_extractmin(a));
	printf("first: %d\n", (int)fh_extractmin(a));
	printf("first: %d\n", (int)fh_extractmin(a));
	for(i = 0; i < 5; i++) {
		x = random()/10000;
		printf("adding: %d\n", x);
		fh_insert(a, (void *)x);
	}
	printf("first: %d\n", (int)fh_extractmin(a));
	printf("first: %d\n", (int)fh_extractmin(a));
	printf("first: %d\n", (int)fh_extractmin(a));
	printf("first: %d\n", (int)fh_extractmin(a));
	printf("first: %d\n", (int)fh_extractmin(a));
#endif

	fh_deleteheap(a);

	return 0;
}
