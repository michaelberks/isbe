#ifndef __ncm_image_grade_h__
#define __ncm_image_grade_h__

//:
// \file
// \brief Structure for holding the image grade of an nailfold capillaroscopy 
//        image
// \author Phil Tresadern
//
// \verbatim
//  Modifications
// \endverbatim

#include <vcl_string.h>

// Forward declarations
class ncm_annotation;

class ncm_image_grade
{
//  INTERFACE

public:
  //  No member variables here, please

  //: Quasi-subjective grading of subject's condition
  enum ImageGrade { GradeUndefined = 0,
                    GradeNormal,
                    GradeEarly,
                    GradeActive,
                    GradeLate,
                    GradeExtreme, // Ungradeable
                    GradePoorQuality, // Ungradeable
                    GradeNonspecific,
    
                    GradeFirst = GradeUndefined,
                    GradeLast = GradeNonspecific };

  //: Default constructor
  ncm_image_grade(ncm_annotation* parent_annotation = NULL);

  //: Destructor
  ~ncm_image_grade();

  //: Image grade as a string
  void set_value(ImageGrade image_grade);
  void set_value(int image_grade);
  void set_from_string(vcl_string image_grade);
  ImageGrade value() const;
  vcl_string as_string() const;


  // IMPLEMENTATION

private:
  //  Members and functions visible only to objects of this class

  //: Classification of the image by the observer
  ImageGrade image_grade_;

  //: Pointer to the owner of this grade
  ncm_annotation* parent_annotation_;
};

#endif // __ncm_image_grade_h__