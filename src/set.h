#ifndef gigide_set_h
#define gigide_set_h

#include <vector>
#include <src/smartptr/ref.h>

namespace gigide {

template <class T>
class SetType {

    public:
        typedef T type;
};

template <class T>
class SetType<const T>
{
    public:
        typedef T type;
};

template <class T>
class Set {
    
    public:
        typedef typename SetType<T>::type val_t;
        typedef typename std::vector<val_t>::const_iterator const_iterator;
        typedef typename std::vector<val_t>::iterator iterator;

    private:
        //this is potentially empty and only exists in cases
        //of type conversion
        std::vector<val_t> items_;

    public:
        Set()
        {
        }

        Set(const std::vector<T>& v)
        {
            //we need to recreate the vector as the new elements
            //inefficient, but otherwise no type safety
            typename std::vector<T>::const_iterator it;
            for (it = v.begin(); it != v.end(); ++it)
                items_.push_back(*it); //do the type conversion
        }

        Set(const Set<val_t>& set) 
        {
            //typename Set<T>::const_iterator it;
            typename std::vector<T>::const_iterator it;
            for (it = set.begin(); it != set.end(); ++it)
                items_.push_back(*it);
        }

        template <class Y>
        Set(const std::vector<Y>& v) 
        {
            //we need to recreate the vector as the new elements
            //inefficient, but otherwise no type safety
            typename std::vector<Y>::const_iterator it;
            for (it = v.begin(); it != v.end(); ++it)
                items_.push_back(*it); //do the type conversion
        }

        template <class Y>
        Set(const Set<Y>& set) 
        {
            //we need to recreate the vector as the new elements
            //inefficient, but otherwise no type safety
            typename Set<Y>::const_iterator it;
            for (it = set.begin(); it != set.end(); ++it)
                items_.push_back(*it);
        }

        template <class Y>
        void
        append(const Y& element)
        {
            items_.push_back(element);
        }

        iterator begin() {return items_.begin();}
        iterator end() {return items_.end();}
        const_iterator begin() const {return items_.begin();}
        const_iterator end() const {return items_.end();}

        int size() const {return items_.size();}

        const T& operator[](int n) const {return items_[n];}

        
};

}

#endif
