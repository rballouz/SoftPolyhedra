/*Define vertex indices*/
#define X 0
#define Y 1
#define Z 2

/*Boolean type*/
typedef enum {FALSE,TRUE} bool;

/*Define Flags*/
#define ONHULL      TRUE
#define REMOVED     TRUE
#define VISIBLE     TRUE
#define PROCESSED   TRUE

/*MACROS*/
#define EXIT_FAILURE 1

/*creates new p type inside if statement, and checks for proper memory allocation*/
#define NEW(p, type) \
        if ((p=(type*) malloc (sizeof(type))) == NULL) {\
        printf ("NEW: Out of Memory!\n");\
        exit(EXIT_FAILURE);\
    }

/*link p to end of head (a circular linked list), if head is NULL set head to p*/
#define ADD(head, p) if (head) {\
        p->next = head;\
        p->prev = head->prev;\
        head->prev = p;\
        p->prev->next = p;\
    }\
    else {\
        head=p;\
        head->next = head->prev = p;\
    }

#define FREE(p) if (p) {free ((char *) p); p = NULL; }

/*Delete p from head*/
#define DELETE(head, p) if (head) {\
    if (head == head->next) \
        head=NULL; \
    else if (p==head) \
        head=head->next; \
    p->next->prev = p->prev; \
    p->prev->next = p->next; \
    FREE(p);\
}

/*Swapping*/
#define SWAP(t,x,y) { t = x; x = y; y = t; }