#ifndef gigide_keyword_hpp
#define gigide_keyword_hpp

namespace gigide {

class KeywordValue;
class KeywordSet;
class KeywordIterator;

typedef boost::intrusive_ptr<KeywordValue> KeywordValuePtr;
typedef boost::intrusive_ptr<KeywordSet> KeywordSetPtr;
typedef boost::intrusive_ptr<KeywordIterator> KeywordIteratorPtr;

typedef boost::intrusive_ptr<const KeywordValue> ConstKeywordValuePtr;
typedef boost::intrusive_ptr<const KeywordSet> ConstKeywordSetPtr;
typedef boost::intrusive_ptr<const KeywordIterator> ConstKeywordIteratorPtr;

}

#endif

