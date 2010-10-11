
#ifndef gigide_permutation_h
#define gigide_permutation_h

#include <vector>
#include <algorithm>

#include <src/deftypes.h>

namespace gigide {

template <class Item>
class PermutationGenerator : public smartptr::Countable
{
    private:
        bool finished_;
        int* index_used_;
        int* indices_;
        std::vector<Item> items_;
        int running_index_;
        int length_;
        int n_;
        void nextindex(int idx);

    public:
        PermutationGenerator(std::vector<Item>& items, int n);
        ~PermutationGenerator();
        bool nextPermutation();
        void generatePermutations(
            std::vector< std::vector<Item> >& all_perms
        );
};

template <class Item>
PermutationGenerator<Item>::PermutationGenerator(std::vector<Item>& items, int n)
{
    //first thing's first! The permutation generator expects sorted items
    sort(items.begin(), items.end());

    running_index_ = n - 1;
    length_ = items.size();

    this->items_ = items;

    index_used_ = new int[length_];
    indices_ = new int[n];
    //set the initial indices - these will get reset later
    for (int i=0; i < n; i++)
    {
        indices_[i] = 0;
    }
    //zero out the rest of the index used array
    for (int i=0; i < length_; i++)
        index_used_[i] = 0;

    //we can't just start the indices at 0,1,2, etc
    //in case of repeats, we may have to shift the indices along
    for (int idx = 0; idx <= running_index_; idx++)
    {
        while ( (indices_[idx] + 1 < length_) && 
                (!index_used_[indices_[idx] + 1]) &&
                (items_[indices_[idx]] == items_[indices_[idx] + 1]) ||
                (index_used_[indices_[idx]]) )
        {
            indices_[idx]++;
        }
        //and now that we have found the final value, claim it!
        index_used_[indices_[idx]] = 1;
    }
    finished_ = false;

    n_ = n;
}

template <class Item>
PermutationGenerator<Item>::~PermutationGenerator()
{
    delete[] index_used_;
    delete[] indices_;
}

template <class Item>
bool
PermutationGenerator<Item>::nextPermutation()
{
    nextindex(running_index_);
    return !finished_;
}

template <class Item>
void
PermutationGenerator<Item>::nextindex(int idx)
{
    if (idx < 0) 
    {
        finished_ = true; 
        return;
    }
    

    int startindex = indices_[idx];
    bool zero_return = false;
    while( 
           index_used_[indices_[idx]] || //if we can shift to other values
           ((indices_[idx] + 1 < length_) && (items_[indices_[idx]] == items_[indices_[idx]+1]) && (!index_used_[indices_[idx] + 1]))
         )
    {
        indices_[idx]++;
        //if we have reached the end of the indices, go back to the beginning
        if (indices_[idx] == length_)
        {
            indices_[idx] = 0;
            index_used_[startindex] = 0;
            zero_return = true;
            nextindex(idx - 1);
        }
    }
    //lay claim to new value
    index_used_[indices_[idx]] = 1;
    //and release claim on old value, assuming we haven't loop back on ourselves
    //if we have looped back, then we have already release the old index
    if (!zero_return) index_used_[startindex] = 0;
    //
    return;
}

template <class Item>
void
PermutationGenerator<Item>::generatePermutations(
    std::vector< std::vector<Item> >& all_perms
)
{
    //we always add a first permutation
    std::vector<Item> first_perm;
    for (int i=0; i < n_; i++)
        first_perm.push_back(items_[indices_[i]]);
    all_perms.push_back(first_perm);

    //now get the rest of the permutations
    while (nextPermutation())
    {
        std::vector<Item> next_perm;
        for (int i=0; i < n_; i++)
            next_perm.push_back(items_[indices_[i]]);
        all_perms.push_back(next_perm);
    }
}


}

#endif
