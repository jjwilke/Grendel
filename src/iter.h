#ifndef gigide_iter_h
#define gigide_iter_h

#include <vector>

template <class T, class U>
class Iterator : public std::vector<T>::const_iterator {

    private:
        typedef typename std::vector<T>::const_iterator parent;

    public:
        Iterator()
        {
        }
        
        Iterator(const parent& iter)
            : parent(iter)
        {
        }

        const U operator*()
        {
            return parent::operator*();
        }
    
};

#endif
