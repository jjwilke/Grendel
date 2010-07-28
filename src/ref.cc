#include <src/ref.h>
#include <iostream>

using namespace std;
using namespace smartptr;

#define concheck cout << "Constructor: " << __FILE__ << " " << __LINE__ << endl

void
smartptr::intrusive_ptr_add_ref(Countable* c)
{
    if (c == NULL) return;

    nref_t count = c->incref();
}

void
smartptr::intrusive_ptr_add_ref(const Countable* c)
{
    Countable* nonconst = const_cast<Countable*>(c);
    intrusive_ptr_add_ref(nonconst);
}

void
smartptr::intrusive_ptr_release(Countable* c)
{
    if (c == NULL) return;
        

    nref_t count = c->decref();

    if (count == 0)
    {
        //cerr << "Deleting " << c << endl;
        delete c;
    }
    else if (count < 0)
    {
        cerr << "Pointer reference to " << c << " already zero.  Cannot release." << endl;
        abort();
    }
    else
    {
        //cerr << "Decrementing " << c << endl;
    }
}

void
smartptr::intrusive_ptr_release(const Countable* c)
{
    Countable* nonconst = const_cast<Countable*>(c);
    intrusive_ptr_release(nonconst);
}

Countable::Countable()
    : refcount_(0)
{
    //cerr << "Allocating " << this << endl;
}

nref_t
Countable::incref()
{
    refcount_++;
    return refcount_;
}

nref_t
Countable::decref()
{
    refcount_--;
    return refcount_;
}

nref_t
Countable::nref() const 
{
    return refcount_;
}

