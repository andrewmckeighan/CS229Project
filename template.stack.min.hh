//@author Andrew McKeighan
#ifndef Stack_G_H
#define Stack_G_H

#include <iostream>
#include <cstdlib>
#include <typeinfo>
#include <fstream>
#include <sstream>
using namespace std;

#include <stdexcept>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Encoded.hh"
#include "wrapper.hh"
#include "prototype.hh"

template <class E>
// This one says that the template parameter E is used as a type in the class.
class Stack
{
  private:
  struct Node // a struct type for Node
  {
     E  *element; // a reference to type E
     Node *next;
     Node(E *elem) : element(elem) {} // inline constructor for reference initialization
  };
  Node  *top;
  // Node  *top = NULL; // not allowed
  
  public:
   Stack();
   ~Stack();
   void  push(E *obj); // a reference argument
   E*    pop();        // a reference return type
   E*    peek();       // a reference return type
   bool  isEmpty();
};

template <class E> // constructor
Stack<E>:: Stack()
{
  top = NULL;
}

template <class E> // desstructor
Stack<E>:: ~Stack()
{
  while ( top != NULL )
   {
     Node *tmp = top;
     top = top->next;
     delete tmp;
   }
}

template <class E>
void Stack<E>:: push(E *obj)
{
  Node *one = new Node(obj);
  one->next = top;
  top = one;
}


template <class E>  
E* Stack<E>:: pop() 
{
  if ( top == NULL )
    fatal("The stack is empty");
  Node *tmp = top;
// E element = top->element; 
  E* element = top->element; 
  top = top->next;
  delete tmp;
  return element;
}

template <class E>
E* Stack<E>:: peek()
{
  if ( top == NULL )
    fatal("The stack is empty");
  E* element = top->element;

  return element;
}

template <class E>
bool Stack<E>:: isEmpty()
{
  return top == NULL;
}


template <class T>
T* findMin(Stack<T> *stack)
{
  Stack<T> tmp;
  if ( stack->isEmpty() )
    fatal("The stack is empty");
  T* ref = stack->pop(); // an error occurs if T ref is used. Why?
  tmp.push(ref);
  T*  min = ref;//does this become a pointer?
  while ( ! stack->isEmpty() )
   {
     T* element = stack->pop(); // the scope of element is limited to one iteration. 
     tmp.push(element);
     if ( element->operator<=(*min) )
        min = element;
   }

  while ( ! tmp.isEmpty() )   // restores stack
     stack->push( tmp.pop() );

  return min;
}

template <class T>
T* findMax(Stack<T> *stack)
{
  Stack<T> tmp;
  if ( stack->isEmpty() )
    fatal("The stack is empty");
  T* ref = stack->pop(); // an error occurs if T ref is used. Why?
  tmp.push(ref);
  T*  max = ref;
  while ( ! stack->isEmpty() )
   {
     T* element = stack->pop(); // the scope of element is limited to one iteration. 
     tmp.push(element);
     if ( max->operator<=(*element) )
        max = element;
   }

  while ( ! tmp.isEmpty() )   // restores stack
     stack->push( tmp.pop() );

  return max;
}


#endif