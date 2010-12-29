#ifndef gigide_input_hpp
#define gigide_input_hpp

namespace gigide {

class InputFile;
class GigideInputFile;

typedef boost::intrusive_ptr<InputFile> InputFilePtr;
typedef boost::intrusive_ptr<GigideInputFile> GigideInputFilePtr;

typedef boost::intrusive_ptr<const InputFile> ConstInputFilePtr;
typedef boost::intrusive_ptr<const GigideInputFile> ConstGigideInputFilePtr;

}

#endif

