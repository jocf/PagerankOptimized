/*
 * header file for SOFT3410 assignment four - "pagerank"
 *
 * DO NOT MODIFY THIS HEADER FILE
 */

#ifndef __PAGERANK_H
#define __PAGERANK_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#define EPSILON 5E-3
#define NAME_SIZE 21
#define BUFFER_SIZE 101

/* forward type definitions */
typedef struct page page;
typedef struct node node;
typedef struct list list;

/* data structure to store page information */
struct page
{
  char name[21]; /* page name */
  int index;     /* index of the page */
  int noutlinks; /* number of outlinks from this page */
  list* inlinks; /* pointer to linked list of pages with inlinks to this page */
};

/* singly linked list node to store each page */
struct node
{
  page* page; /* pointer to page data structure */
  node* next; /* pointer to next page in list */
};

/* singly linked list to store all pages */
struct list
{
  node* head; /* pointer to the head of the list */
  node* tail; /* pointer to the tail of the list */
  int length; /* length of the entire list */
};


/* ========== FUNCTION PROTOTYPES ========== */

static void die(list* plist);
static void page_destroy(page* p);
static page* page_create(char* name, int index);
static list* page_list_create(void);
static void page_list_destroy(list* plist);
static node* page_list_find(list* plist, char* name);
static void _read_edges(char* buffer, list* plist, int* nedges);
static void _read_page_list(char* buffer, list** plist, int npages);
static void read_input(list** plist, int* ncores, int* npages, int* nedges, double* damping_factor);

/* ========================================== */

/**
 * Outputs 'error' to stdout
 * cleans up the given page list and exits with non-zero status
 */
static void die(list* plist)
{
  printf("error\n");
  page_list_destroy(plist);
  exit(1);
}

/**
 * Creates a page struct given a name and index
 *     > Returns NULL when name is too long or malloc fails
 */
static page* page_create(char* name, int index)
{
  if (strlen(name) >= NAME_SIZE) /* page name too long */
    return NULL;

  page *p = (page *) malloc(sizeof(page));
  if (p == NULL) /* failed to allocate memory */
    return NULL;

  strcpy(p->name, name);
  p->index = index;
  p->noutlinks = 0;
  p->inlinks = NULL;

  return p;
}

/**
 * Free the memory of used by a dynamically allocated specified page
 */
static void page_destroy(page* p)
{
  if (p == NULL) /* empty page */
    return;

  list* inlinks = p->inlinks;
  if (inlinks != NULL) /* empty inlinks */
  {
    node *next;
    node *curr;

    curr = inlinks->head;
    while (curr != NULL)
    {
      next = curr->next;
      free(curr);
      curr = next;
    }
    free(inlinks);
  }
  free(p);
}

/**
 * Create an empty page linked list
 *     > Returns the created list, otherwise NULL if malloc fails
 */
static list* page_list_create(void)
{
  list *plist = (list *) malloc(sizeof(list));

  if (plist == NULL) /* failed to allocate memory */
    return NULL;

  plist->length = 0;
  plist->head = NULL;
  plist->tail = NULL;

  return plist;
}

/* frees a page linked list and all of the pages it contains.
 * assumes pl was allocated by malloc
 */
static void page_list_destroy(list* plist)
{
  node* next;
  node* current;

  if (plist == NULL) /* null page list */
    return;

  current = plist->head;
  while (current != NULL)
  {
    next = current->next;
    page_destroy(current->page);
    free(current);
    current = next;
  }

  free(plist);
}

/* adds a page p to the start of page linked list pl.
 * returns a pointer to the new front of the linked list
 * returns NULL if malloc fails
 */
static node* page_list_add_front(list* plist, page* p)
{
  node* new_node = (node *) malloc(sizeof(node));

  if (new_node == NULL) /* failed to allocate memory */
    return NULL;

  /* assign the page */
  new_node->page = p;

  if (plist->head == NULL) /* empty list */
  {
    plist->head = new_node;
    plist->tail = new_node;

    plist->head->next = NULL;
    plist->tail->next = NULL;
  }

  else /* insert into the front */
  {
    new_node->next = plist->head;
    plist->head = new_node;
  }

  plist->length++;

  /* front of the list */
  return plist->head;
}

/* adds a page p to the end of a page linked list pl
 * returns the new end of the linked list (where next = NULL)
 * returns NULL if malloc fails
 */
static node* page_list_add_end(list* plist, page* p)
{
  if (plist == NULL) /* null list */
    return NULL;

  node* new_node = (node *) malloc(sizeof(node));

  if (new_node == NULL) /* failed to allocate memory */
    return NULL;
  else /* assign the page */
    new_node->page = p;

  if (plist->head == NULL) /* empty list */
  {
    plist->head = new_node;
    plist->tail = new_node;

    plist->head->next = NULL;
    plist->tail->next = NULL;
  }

  else /* insert into the end */
  {
    plist->tail->next = new_node;
    plist->tail = new_node;
    new_node->next = NULL;
  }

  plist->length++;

  /* end of the list */
  return plist->tail;
}

/* given a page linked list and a name, performs a forward iteration through
 * the linked list and finds the item containing a page with the same name.
 *
 * returns NULL if no matching page can be found
 */
static node* page_list_find(list* plist, char* name)
{
  if (plist == NULL) /* null list */
    return NULL;

  node* curr = plist->head;
  while (curr != NULL)
  {
    if (strcmp(curr->page->name, name) == 0)
      return curr;
    curr = curr->next;
  }
  return NULL;
}

/* helper function to read the list of pages into a page_list
 *
 * buffer is used by fgets to read input (assumed to be of sufficient size)
 *
 * ppl is a pointer to a pointer to a page_list, so that this function can
 * modify the pointer value to point to the filled page list
 *
 * die() is called if there are any input errors
 */
static void _read_page_list(char* buffer, list** plist, int npages)
{
  page* p;
  char name[NAME_SIZE];

  if (plist == NULL) /* null list */
    die(NULL);

  *plist = page_list_create();

  for (int i = 0; i < npages; i++)
  {
    if (fgets(buffer, BUFFER_SIZE, stdin) == NULL
        || sscanf(buffer, "%20s\n", name) != 1)
      die(*plist);

    if ((p = page_create(name, i)) == NULL)
      die(*plist);

    if (page_list_add_end(*plist, p) == NULL)
      die(*plist);
  }
}

/* helper function to read the list of edges into pages in a page_list
 *
 * buffer is used by fgets to read input (assumed to be of sufficient size)
 *
 * die() is called if there are any input errors
 */
static void _read_edges(char* buffer, list* plist, int* nedges)
{
  node* node1;
  node* node2;

  char name1[NAME_SIZE];
  char name2[NAME_SIZE];
  char excess[NAME_SIZE];

  if (fgets(buffer, BUFFER_SIZE, stdin) == NULL
      || sscanf(buffer, "%d %s\n", nedges, excess) != 1)
    die(plist);

  for (int i = 0; i < *nedges; i++)
  {
    if (fgets(buffer, BUFFER_SIZE, stdin) == NULL
        || sscanf(buffer, "%20s %20s %s\n", name1, name2, excess) != 2)
      die(plist);

    node1 = page_list_find(plist, name1);
    node2 = page_list_find(plist, name2);

    if (node1 == NULL || node2 == NULL) /* undefined pages */
      die(plist);

    /* Make sure we can add links into a page */
    if (node2->page->inlinks == NULL &&
        (node2->page->inlinks = page_list_create()) == NULL)
      die(plist);

    /* Add the link to the page */
    if (page_list_add_front(node2->page->inlinks, node1->page) == NULL)
      die(plist);

    node1->page->noutlinks++;
  }
}

/* function to read input in the correct format and handle input errors */
static void read_input(list** plist, int* ncores, int* npages,
    int* nedges, double* dampener)
{
  char buffer[BUFFER_SIZE];

  /* check for invalid input */

  if (fgets(buffer, BUFFER_SIZE, stdin) == NULL
      || sscanf(buffer, "%d\n", ncores) != 1
      || *ncores == 0)
    die(*plist);

  if (fgets(buffer, BUFFER_SIZE, stdin) == NULL
      || sscanf(buffer, "%lf\n", dampener) != 1
      || *dampener < 0 || fabs(*dampener) > 1)
    die(*plist);

  if (fgets(buffer, BUFFER_SIZE, stdin) == NULL
      || sscanf(buffer, "%d\n", npages) != 1
      || *npages == 0)
    die(*plist);

  /* read in the list of pages and edges once we have the main parameters */

  _read_page_list(buffer, plist, *npages);
  _read_edges(buffer, *plist, nedges);
}

#endif
