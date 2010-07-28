
#ifndef gigide_ref_h
#define gigide_ref_h

#include <boost/intrusive_ptr.hpp>

namespace smartptr {

typedef unsigned int nref_t;

class Countable {
    
    private:
        nref_t refcount_;

    public:
        Countable();

        nref_t incref();

        nref_t decref();

        nref_t nref() const;

        virtual ~Countable(){}

};

void
intrusive_ptr_add_ref(const Countable* c);

void
intrusive_ptr_add_ref(Countable* c);

void
intrusive_ptr_release(const Countable* c);

void
intrusive_ptr_release(Countable* c);

template <class T>
class RefPtr : public boost::intrusive_ptr<T> {

    protected:
        typedef boost::intrusive_ptr<T> parent;

    public:
        RefPtr() : parent(NULL)
        {
        }

        RefPtr(T* p) : parent(p)
        {
        }

        RefPtr(const parent& p) : parent(p)
        {
        }

        RefPtr(const RefPtr<T>& p) : parent(p)
        {
        }

        template <class Y>
        RefPtr(const RefPtr<Y>& p) : parent(p)
        {
        }

        bool null(){return parent::get() == NULL;}

        bool nonnull(){return parent::get() != NULL;}

};

}

#endif
