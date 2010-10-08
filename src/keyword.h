
#ifndef intdif_keyword_h
#define intdif_keyword_h

#include <iostream>
#include <sstream>
#include <map>
#include <src/pyregexp.h>
#include <src/printstream.h>
#include <src/defines.h>
#include <src/deftypes.h>

namespace gigide {

/** 
    @class KeywordValue

    Encapsulates a keyword value.  This only knows its value, not its name, i.e.
	suppose I specify debug = 1.  This keeps track of the value 1, not the name debug. 
*/
class KeywordValue : public smartptr::Countable {
    
    public:
		//typedef enum { DOUBLE, STRING, INTEGER, BOOLEAN, VECTOR_DOUBLE } datatype;
    
    private:
		/** 
            The text associated with a given keyword.  Based on the method called later,
		 	the text will be cast as integer, double, etc. 
        */
        std::string valtext_;

        /** 
            Pop off next value in a list 
        */
        template <class T> 
        T pop();

        /**
            Pop the next few items as a vector of type T
            @param n The number of elements in the vector
            @return The vector of elements
        */
        template <class T> std::vector<T> popVector(int n);

    public:
		/** 
            Constructor
			@param val 	The keyword value to store
		*/
        KeywordValue(const std::string& val);

		/** 
            Casts the option text as a double
			@return The keyword value as a double
		*/
        double getValueDouble() const;

		/** 
            Sends back just the string representation of the value
			@return The keyword value as a string
		*/
        std::string getValueString() const;

		/** 
            Casts the option text as an integer
			@return The keyword value as an integer 
		*/
        int getValueInteger() const;

		/** 
            Casts the option text as a boolean
			@return The keyword value as a boolean
		*/
        bool getValueBoolean() const;

		/** 
            Casts the option text as a vector of doubles
			@return The keyword value as a vector of doubles
		*/
        std::vector<double> getValueVectorDouble() const;

		/** 
            Casts the option text as a vector of strings
			@return The keyword value as a vector of strings
		*/
        void getValueVectorString(std::vector<std::string>& vec) const;

		/** 
            Casts the option text as a vector of doubles
			@return The keyword value as a vector of doubles
		*/
        std::vector<int> getValueVectorInteger() const;

        /**
            @return The number of elements remaining
        */
        int count() const;

        /**
            @return The next item in the list as an integer
        */
        int popInteger();

        /**
            @return The next item in the list as a double
        */
        double popDouble();

        /**
            @return The next item in the list as a string
        */
        std::string popString();

        /**
            @param n The number of elements to pop
            @return Vector of strings of length n
        */
        std::vector<std::string> popVectorString(int n);

        /**
            @param n The number of elements to pop
            @return Vector of integers of length n
        */
        std::vector<int> popVectorInteger(int n);

        /**
            @param n The number of elements to pop
            @return Vector of doubles of length n
        */
        std::vector<double> popVectorDouble(int n);

        /**
        */
		~KeywordValue();
};

/** 
    Encapsulates a set of keyword-value pairs. 

    @class KeywordSet
*/
class KeywordSet : public smartptr::Countable {

    public:
		/** 
            Constructor
			@param inputtest The text with the user specified options
			@param defaults  The text defining default values that user may or may not have overridden
		*/
        KeywordSet(const std::string& inputtext, const std::string& defaults);

		/** 
            Constructor
			@param inputtest The text with the user specified options
		*/
        KeywordSet(const std::string& inputtext);

		/** 
            Returns a keyword value object with the given name
			@param key The keyword name
			@return The keyword value
		*/
        static KeywordValuePtr getKeyword(const std::string& key);

        /**
            Prints the universal (static) keyword set

            @param os The output stream to print to
        */
        static void printKeywords(std::ostream& os = std::cout);

        /**
            @param os The output stream to print to
        */
        void print(std::ostream& os = std::cout) const;

        /**
            @return The number of keywords in the set
        */
        int size() const;

        /**
            @param key The name of the keyword to return
            @return The keyword value object 
        */
        KeywordValuePtr get(const std::string& key) const;

        /**
            @param i The keyword number 
            @return The name of the keyword
            @throw GigideException If index is invalid
        */
        std::string getKeyname(int i) const;

        /**
        */
        ~KeywordSet(){}

    private:
		/** 
            The std::map storing all keyword, value pairs 
        */
        std::map<std::string, KeywordValuePtr > keymap_;
        
        /** 
            Orders the keynames depending on order input in file 
        */
        std::map<int, std::string> keynames_;

        static KeywordSetPtr keyset_;

		/** 
            Given a text input, adds all keyword value pairs found therein
			@param text The text specifying keywords
            @param check 
		*/
        void addKeywords(const std::string& text, bool check = false);

};

/**
    Iterates a a KeywordSet

    @class KeywordIterator
*/
class KeywordIterator : public smartptr::Countable {

    private:
        /**
            The keyword set to be iterated
        */
        KeywordSetPtr keyset_;

        /**
            The current iteration index
        */
        int ival_;
    
    public:
        /**
            The text of the section to create keywords for
        */
        KeywordIterator(const std::string& section);

        /**
            @return The current keyword in the iteration
        */
        KeywordValuePtr getKeyword() const;

        /**
            @return The name of the current keyword
        */
        std::string getName() const;

        /**
            Move to the next keyword.
        */
        void next();

        /**
            @return Whether the iteration is finished
        */
        bool finished() const;

        /**
            Start the iterations going
        */
        void start();

};


/**
    @class Encapsulates repetitive gigide keyword function calls that require processing of the data
*/
class GigideKeyword {

    public:

        static std::vector<double> getDisplacementSizes(int ncoord);

        static std::string getBondUnits();

        static std::string getEnergyUnits();

};

}

#endif
