#include "gigide.h"

using namespace std;
using namespace gigide;
using namespace smartptr;

KeywordValue::KeywordValue(const std::string& val)
    : valtext_(val)
{
}

double
KeywordValue::getValueDouble() const
{
	double newval;
	stringstream sstr(valtext_);
	sstr >> newval;
	return newval;
}

string 
KeywordValue::getValueString() const
{
	return valtext_;
}

size_t
KeywordValue::getValueMemory() const
{
    double prefactor;
    size_t inc;

    std::string incstr;
    stringstream sstr(valtext_);

    sstr >> prefactor;
    sstr >> incstr;

    if      (incstr == "b")
        inc = 1;
    else if (incstr == "kb")
        inc = 1 << 10;
    else if (incstr == "mb")
        inc = 1 << 20;
    else if (incstr == "gb")
        inc = 1 << 30;
    else if (incstr == "mw")
        inc = 8 * (1 << 20);
    else
    {
        except("Invalid memory storage specifier.  I need something like B, KB, MB, GB, or MW");
    }

    double dblmem = prefactor * inc;
    size_t mem = (size_t) dblmem;
    return mem;
}

int
KeywordValue::count() const
{
    stringstream sstr(valtext_);
    int cnt = 0;
    string nextval;
    while (!sstr.eof())
    {
        sstr >> nextval;
        ++cnt;
    }
    return cnt;
}

int
KeywordValue::getValueInteger() const
{
	int newval;
	stringstream sstr(valtext_);
	sstr >> newval;
	return newval;
}

bool
KeywordValue::getValueBoolean() const
{
    if (valtext_ == "true")
        return true;
    else if (valtext_ == "false")
        return false;
    else
    {
        except(stream_printf("Invalid boolean value %s", valtext_.c_str()));
    }

	bool newval;
	stringstream sstr(valtext_);
	sstr >> newval;
	return newval;
}

void
KeywordValue::getValueVectorString(vector<string>& vec) const
{
	findmatch(vec, "([a-zA-Z\\d]+)", valtext_, FindAll);
}

vector<int>
KeywordValue::getValueVectorInteger() const
{
	vector<int> vec;
	stringstream sstr(valtext_);
	while (!sstr.eof())
	{
		int nextval;
		sstr >> nextval;
		vec.push_back(nextval);
	}
	return vec;
}

vector<double>
KeywordValue::getValueVectorDouble() const
{
	vector<double> dblvec;
	stringstream sstr(valtext_);
	while (!sstr.eof())
	{
		double nextval;
		sstr >> nextval;
		dblvec.push_back(nextval);
	}
	return dblvec;
}

template <
    class T
> T 
KeywordValue::pop()
{
    T val;
    stringstream sstr(valtext_);
    sstr >> val;
    sstr.ignore(256, ' ');

    stringbuf scratch;
    sstr.get(scratch);
    valtext_ = scratch.str();

    return val;
}

template <
    class T
> vector<T>
KeywordValue::popVector(int n)
{
    vector<T> vec;
    for (int i=0; i < n; ++i)
        vec.push_back(pop<T>());
    return vec;
}

int
KeywordValue::popInteger()
{
    return pop<int>();
}

double
KeywordValue::popDouble()
{
    return pop<double>();
}

string
KeywordValue::popString()
{
    return pop<string>();
}

vector<int>
KeywordValue::popVectorInteger(int n)
{
    return popVector<int>(n);
}

vector<double>
KeywordValue::popVectorDouble(int n)
{
    return popVector<double>(n);
}

vector<string>
KeywordValue::popVectorString(int n)
{
    return popVector<string>(n);
}

KeywordValue::~KeywordValue()
{
	//do not delete anything... when the data is returned it must remain
}

KeywordSetPtr KeywordSet::keyset_;

KeywordSet::KeywordSet(const std::string& inputtext, const std::string& defaults)
{
	//first add the defaults
	addKeywords(defaults);
	//now, go through and set all the options explicitly set
	addKeywords(inputtext, true);

    keyset_ = this;
}

KeywordSet::KeywordSet(const std::string& inputtext)
{
	addKeywords(inputtext, false); //do not validate
}

KeywordValuePtr
KeywordSet::getKeyword(const std::string& key)
{
    return keyset_->get(key);
}

KeywordValuePtr
KeywordSet::get(const std::string& key) const
{
    map<string, KeywordValuePtr>::const_iterator it;
    it = keymap_.find(key);

	if (it == keymap_.end())
	{
        stringstream sstr;
		sstr << "No value for keyword " << key << " has been given!" << endl;
        print(sstr);
        except(sstr.str());
	}

	return it->second;
}

void
KeywordSet::print(ostream& os) const
{
    map<string, KeywordValuePtr>::const_iterator it;
    for (it = keymap_.begin(); it != keymap_.end(); ++it)
    {
        os << stream_printf("%30s = %s", it->first.data(), it->second->getValueString().data()) << endl;
    }
}

void
KeywordSet::printKeywords(ostream& os)
{
    keyset_->print(os);
}

void
KeywordSet::addKeywords(const std::string& text, bool check)
{
	stringstream defaultstr(text);
	string key_regexp = "\\s*([a-zA-Z \\d_]*[a-zA-Z\\d])\\s*[=]\\s*([a-zA-Z\\-\\d][+a-zA-Z \\-\\d.]*[a-zA-Z\\d]*)\\s*";
    string blankline_regexp = "\\s+";
	char* line = new char[MAX_LINE_SIZE];

    int keynum = 0;
	while (!defaultstr.eof())
	{
		defaultstr.getline(line, MAX_LINE_SIZE);
        string nextline = line; 

        if (!nextline.size())
            continue;

		if (has_regexp_match(key_regexp, nextline, LowerCase))
		{
			vector<string> keyoption_pair; findmatch(keyoption_pair, key_regexp, line, LowerCase | StripWhitespace); //2 groups
			string key = keyoption_pair[0];
			string option = keyoption_pair[1];
			KeywordValuePtr keyval = new KeywordValue(option);
            //this should have overwritten the default value... grab it to make sure it is nonnull
            KeywordValuePtr test = keymap_[key];
            if (check && !test)
            {
                keymap_.erase(keymap_.find(key));
                stringstream sstr;
                sstr << "Keyword " << key << " is not a valid keyword name. Valid keywords are:" << endl;
                sstr << endl << endl;
                print_keys<KeywordValuePtr >(keymap_, sstr);
                string msg = sstr.str();
                except(msg);
            }
			keymap_[key] = keyval;
            keynames_[keynum] = key;
            ++keynum;
		}
        else if (has_regexp_match(blankline_regexp, nextline))
        {
            //cout << "blank line" << endl;
        }
        else
        {
            stringstream sstr;
            sstr << "Input section not properly formatted. Please revise the following:" << endl;
            sstr << line << endl;
            except(sstr.str());
        }
	}

    delete[] line;
}

string
KeywordSet::getKeyname(int i) const
{
    map<int, string>::const_iterator it;
    it = keynames_.find(i);
    if (it == keynames_.end())
        except(stream_printf("Element %d invalid for keyword vector of size %d", i, keymap_.size()));

    return it->second;
}

int
KeywordSet::size() const
{
    return keymap_.size();
}

KeywordIterator::KeywordIterator(
    const std::string& section
)
{
    keyset_ = new KeywordSet(section);
    ival_ = 0;
}

KeywordValuePtr
KeywordIterator::getKeyword() const
{
    return keyset_->get( getName() );
}

string
KeywordIterator::getName() const
{
    return keyset_->getKeyname(ival_);
}

void
KeywordIterator::start()
{
    ival_ = 0;
}

bool
KeywordIterator::finished() const
{
    return ival_ == keyset_->size();
}

void
KeywordIterator::next()
{
    ++ival_;
}

vector<double>
GigideKeyword::getDisplacementSizes(int ncoord)
{
    vector<double> disp_sizes = KeywordSet::getKeyword("dispsizes")->getValueVectorDouble();
    if (disp_sizes.size() == 1) //only a single value was given, replicate it for all displacements
    {
        for (int i=1; i < ncoord; i++)
            disp_sizes.push_back( disp_sizes[0] );
    }
    return disp_sizes;
}

string
GigideKeyword::getBondUnits()
{
    return KeywordSet::getKeyword("bond units")->getValueString();
}

string
GigideKeyword::getEnergyUnits()
{
    return KeywordSet::getKeyword("energy units")->getValueString();
}



