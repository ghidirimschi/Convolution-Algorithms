/* based on listTypeAndFunctions.c by Gerard Renardel, 31 January 2018 */

#ifndef LIBLIST_H
#define LIBLIST_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct ListNode* List;

struct ListNode {
  int item;
  List next;
};

List newEmptyList();

int isEmptyList(List li);
void listEmptyError();
List addItem(int n, List li);
int firstItem(List li);
List removeFirstNode(List li);
void freeList(List li);
void visit(List li);
void visitList(List li);
void visitListRec(List li);
void listTooShort();
int itemAtPos(List li, int p);
List addItemAtPos(List li, int n, int p);
List addItemAtPosIt(List li, int n, int p);
List removeItem(List li, int n);
List removeItemIt(List li, int n);

#endif
