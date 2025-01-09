#ifndef CMTOBJECT_H
#define CMTOBJECT_H 1
/*
  * This object is meant to be a base class for any other multithreaded object
  * It's only point is to manage this nb_threads static member variable
  * That is, all classes inheriting from MTObject will share this variable
*/

// As it is a static variable we have to declare it outside of the class
// Also, I think it is better to initialise it at 1 to avoid unitialisation issue
//unsigned short CMTObject::nb_threads = 1;

class CMTObject
{
public:
  static unsigned short nb_threads;
};


#endif //MTOBJECT_H
